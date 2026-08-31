#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# GDSC-MGX Pipeline (BioBakery Implementation)
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#≠≠≠≠≠ To Do ≠≠≠≠≠#
# Method for extracting uniprot IDs for groups of interest and then searching in Centrifuger output
# Need some way of memory management for temporary humann bowtie tsvs once queried (delete?)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SET GLOBAL SCOPE PYTHON VARIABLES (EXECUTED BEFORE SNAKEMAKE)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
import sys
import pandas as pd

#----- Set path to config file
configfile: "config.yaml"

#----- read in sample data and extract sample names
samples_df = pd.read_csv(config["sample_csv"])

#----- sample_id has to exist before it can be used as the index, otherwise
# pandas raises a bare KeyError that does not say which file is at fault
if "sample_id" not in samples_df.columns:
    raise ValueError(
        "sample sheet {} has no 'sample_id' column. Found: {}".format(
            config["sample_csv"], ", ".join(map(str, samples_df.columns)))
    )

samples_df = samples_df.set_index("sample_id", drop=False)
sample_list = list(samples_df['sample_id'])

#----- Set directories
projDir = config["proj_dir"] + '/'
outDir = config["outDir"] + '/'
raw = projDir + "raw_data/"

#----- Host read filtering
kneaddata = outDir + "kneaddata/"

#----- MetaPhlAn, per sample outputs then the merged tables built from them
mpaDir     = outDir + "metaphlan_taxonomy/"
sams       = mpaDir + "sams/"
profiles   = mpaDir + "profiled_metagenomes/"
mpaMerged  = mpaDir + "merged_abundances/"
bt2logs    = mpaDir + "bowtie2_logs/"

#----- MMUPHin batch correction. run_mmuphin.R emits the corrected table, the
# before/after PCoA and MMUPHin's own diagnostic pdf off a single path prefix,
# so they land together here rather than being split across results and plots.
mmuphinDir = outDir + "mmuphin/"

#----- HUMAnN, same shape as MetaPhlAn above
humannDir    = outDir + "humann_function/"
humannBt2    = humannDir + "bowtie2_aligned/"
humannProf   = humannDir + "profiled_functions/"
humannMerged = humannDir + "merged_abundances/"

#----- Figures and the tables written alongside them
plots = outDir + "plots/"

#----- Create necessary directories
tmp = config["tmp_dir"]
if not os.path.isdir(tmp):
    os.mkdir(tmp)

#----- Skip host read filtering and profile the sample sheet fastqs directly
skip_kneaddata = config.get("skip_kneaddata", False)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# VALIDATE THE SAMPLE SHEET
#
# The sheet needs sample_id and fastq_1 always, plus fastq_2 when layout is
# paired. Everything else in it is sample metadata, free for the user to name,
# and available to downstream rules as a batch or covariate variable.
#
# Problems are collected and reported together rather than one per run, since
# fixing them one at a time across cluster submissions is slow.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#----- Columns the pipeline reserves. Structural problems stop here, because
# nothing else is worth checking until the sheet has the right shape.
read_cols = ["fastq_1"] + (["fastq_2"] if config["layout"] == "paired" else [])
required_cols = ["sample_id"] + read_cols

absent = [c for c in required_cols if c not in samples_df.columns]
if absent:
    raise ValueError(
        "sample sheet {} is missing column(s) {}. layout is '{}', so it needs {}.".format(
            config["sample_csv"], ", ".join(absent), config["layout"],
            ", ".join(required_cols))
    )

#----- Value level problems are collected so they can be fixed in one pass
ids = samples_df["sample_id"].astype(str).str.strip()
problems = []

if ids.eq("").any() or ids.str.lower().eq("nan").any():
    problems.append("sample_id is blank for at least one row")
if ids.duplicated().any():
    problems.append("duplicate sample_id(s): " +
                    ", ".join(sorted(set(ids[ids.duplicated()]))))

# sample_id becomes part of output filenames, so it cannot carry separators
odd = sorted(ids[ids.str.contains(r"[\s/\\]")])
if odd:
    problems.append("sample_id contains whitespace or a slash, which breaks "
                    "output paths: " + ", ".join(odd))

for _c in read_cols:
    blank = sorted(ids[samples_df[_c].astype(str).str.strip().isin(["", "nan"])])
    if blank:
        problems.append("{} is blank for sample(s): {}".format(_c, ", ".join(blank)))

if problems:
    raise ValueError("problems in sample sheet {}:\n  - {}".format(
        config["sample_csv"], "\n  - ".join(problems)))

#----- Anything not reserved is metadata, usable downstream as a batch or
# covariate variable
metadata_cols = [c for c in samples_df.columns if c not in required_cols]
if metadata_cols:
    sys.stderr.write("Sample sheet metadata columns available: {}\n".format(
        ", ".join(metadata_cols)))

#----- MMUPHin batch correction is opt in. It needs a batch column in the sample
# sheet and not every project has one, so a config without run_mmuphin keeps the
# rule out of the DAG entirely.
run_mmuphin = config.get("run_mmuphin", False)
mmuphin_ranks = config.get("mmuphin_ranks", ["species"])

#----- Covariates are a list in config.yaml, but a comma separated string is an
# easy thing to write by accident, so take either. They reach run_mmuphin.R as
# separate trailing arguments, which is why an empty list can simply disappear
# from the command line without breaking the argument positions.
mmuphin_covariate_list = config.get("mmuphin_covariates") or []
if isinstance(mmuphin_covariate_list, str):
    mmuphin_covariate_list = mmuphin_covariate_list.split(",")
mmuphin_covariate_list = [c.strip() for c in mmuphin_covariate_list if str(c).strip()]
mmuphin_covariates = " ".join(mmuphin_covariate_list)

#----- Fail here rather than several hours into a run, inside an R script, on a
# compute node. run_mmuphin.R checks this too, but by then the taxonomy and
# function rules have already burned their walltime.
if run_mmuphin:
    mmuphin_cols = [config["mmuphin_batch"]] + mmuphin_covariate_list
    absent_mmuphin = [c for c in mmuphin_cols if c not in samples_df.columns]
    if absent_mmuphin:
        raise ValueError(
            "run_mmuphin is on but sample sheet {} has no column(s) {}. "
            "Available metadata columns: {}".format(
                config["sample_csv"], ", ".join(absent_mmuphin),
                ", ".join(metadata_cols) or "none")
        )
    if not mmuphin_covariate_list:
        sys.stderr.write(
            "WARNING: run_mmuphin is on with no mmuphin_covariates. Correcting on "
            "'{}' alone removes any biological signal that tracks the batch.\n".format(
                config["mmuphin_batch"]))

#----- Report every unreadable read file at once, not one per failed job
unreadable = ["  {} {}: {}".format(_s, _c, samples_df.loc[_s, _c])
              for _s in sample_list for _c in read_cols
              if not os.access(str(samples_df.loc[_s, _c]), os.R_OK)]
if unreadable:
    sys.stderr.write(
        "WARNING: {} file(s) in {} are missing or unreadable. Jobs needing them "
        "will fail, repeatedly if restart-times is set:\n{}\n".format(
            len(unreadable), config["sample_csv"], "\n".join(unreadable)))

#----- metaphlan will not start if its bowtie2out file already exists
# (metaphlan.py: "BowTie2 output file detected ... Exiting"). Snakemake deletes
# outputs of failed jobs, but a job killed outright, by a timeout or a node
# failure, can leave one behind. Every retry then dies the same way, which with
# restart-times set means a pile of identical failures.
# Only flag the blocking case: a bowtie2out with no profile beside it, meaning
# rule taxonomy still has to run for that sample.
stale_bt2 = [_s for _s in sample_list
             if os.path.exists(bt2logs + str(_s) + ".bowtie2out.txt")
             and not os.path.exists(profiles + str(_s) + "-profiled_metagenome.txt")]
if stale_bt2:
    sys.stderr.write(
        "WARNING: {} leftover bowtie2out file(s) with no matching profile. "
        "metaphlan will refuse to re-run these samples until they are removed:\n".format(
            len(stale_bt2)))
    for _s in stale_bt2:
        sys.stderr.write("  {}\n".format(bt2logs + str(_s) + ".bowtie2out.txt"))
    sys.stderr.write("  remove with: rm {}*.bowtie2out.txt\n".format(bt2logs))

#----- Where the profiling rules get their reads
# With skip_kneaddata on, rule filter never runs and the fastq named in the
# sample sheet is profiled as-is. Otherwise the kneaddata output is used.
# Both rule taxonomy and rule function read this, so they stay in step.
def profiling_reads(wildcards):
    if skip_kneaddata:
        return samples_df.loc[wildcards.sample, "fastq_1"]
    return kneaddata + wildcards.sample + ".fastq.gz"

#----- Paired data loses its R2 in skip mode, so say so rather than hide it
if skip_kneaddata and config["layout"] == "paired":
    sys.stderr.write(
        "WARNING: skip_kneaddata is on and layout is 'paired'. Only fastq_1 is "
        "profiled, fastq_2 is not used. humann takes a single input file, so "
        "there is nowhere to put the second mate without concatenating.\n"
    )

#----- Define all outputs
all_outputs = [
        #----- Taxonomy outputs
        expand(profiles + "{sample}-profiled_metagenome.txt", sample = sample_list),
        expand(sams + "{sample}.sam.bz2", sample = sample_list),
        expand(bt2logs + "{sample}.bowtie2out.txt", sample = sample_list),
        mpaMerged + "merged_abundance_table.txt",
        mpaMerged + "merged_genus_abundance_table.txt",
        mpaMerged + "merged_species_abundance_table.txt",
        #----- Function outputs
        expand(humannProf + "{sample}-humann.log", sample = sample_list),
        expand(humannProf + "{sample}_genefamilies.tsv", sample = sample_list),
        expand(humannProf + "{sample}_pathcoverage.tsv", sample = sample_list),
        expand(humannProf + "{sample}_pathabundance.tsv", sample = sample_list),
        expand(humannBt2 + "{sample}_bowtie2_aligned.tsv", sample = sample_list),
        humannMerged + "merged_genefamilies.txt",
        humannMerged + "merged_pathabundances.txt",
        humannMerged + "merged_pathcoverages.txt",
        humannMerged + "merged_genefamilies_cpm.txt",
        humannMerged + "merged_pathabundances_cpm.txt",
        humannMerged + "merged_pathcoverages_cpm.txt",
        #----- Functional results
        humannMerged + "merged_genefamilies_cpm_named.txt",
        humannMerged + "merged_pathabundances_cpm_named.txt",
        humannMerged + "merged_pathcoverages_cpm_named.txt",
        humannMerged + "merged_taxonomic_profiles.tsv"
]

#----- The filtered fastqs and the diagnostics plot are kneaddata's own
# products, so they are only asked for when kneaddata actually runs.
if not skip_kneaddata:
    all_outputs.append(expand(kneaddata + "{sample}.fastq.gz", sample = sample_list))
    all_outputs.append(plots + "KneadData_Diagnostics.png")

#----- Only the corrected table is asked for. See rule batch_correct for why the
# figure is not a declared output.
if run_mmuphin:
    all_outputs.append(expand(mmuphinDir + "{rank}_mmuphin_corrected_pct.txt",
                              rank = mmuphin_ranks))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SNAKEMAKE RULES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#    
rule all:
    input:
        all_outputs
    output:
        plots + "Species_Shannon_Diversity_Per_Sample.png",
        plots + "Genus_Shannon_Diversity_Per_Sample.png"
    params:
        sample_csv = config["sample_csv"]
    conda: 
        "env_config/r-plotting.yaml"
    resources: cpus="10", maxtime="7:00:00", mem_mb="60gb"
    shell: """
    
        #----- Run Alpha diversity
        # merged-abundance and plots dirs are passed explicitly, they hang off out_dir
        # and are not relative to where snakemake runs from
        Rscript scripts/plot_alpha.R {params.sample_csv} {mpaMerged} {plots}

    """


#----- Filter host reads with Kneaddata
rule filter:
    output:
        kneaddata + "{sample}.fastq.gz"
    params:
        sample = lambda wildcards: wildcards.sample,
        R1 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_1"],
        R2 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_2"] if config["layout"]=="paired" else "None",
        layout = config["layout"],
        kneaddata = config["kneaddata_path"],
        bowtie_ref = config["bowtie_ref"]
    conda:
        "env_config/kneaddata.yaml",
    resources: cpus="10", maxtime="7:00:00", mem_mb="60gb",
    threads: 16
    shell: """
        echo {params.sample}

        TRIM_JAR=$(ls "$CONDA_PREFIX"/share/trimmomatic*/trimmomatic.jar 2>/dev/null | head -n 1 || true)
        if [ -z "$TRIM_JAR" ]; then
            echo "ERROR: no trimmomatic.jar under $CONDA_PREFIX/share" >&2
            echo "       Run snakemake with --use-conda so env_config/kneaddata.yaml is active." >&2
            exit 1
        fi
        TRIM_DIR="{tmp}/trimmomatic_{params.sample}"
        mkdir -p "$TRIM_DIR"
        ln -sf "$TRIM_JAR" "$TRIM_DIR/trimmomatic.jar"
        echo "Using trimmomatic jar: $TRIM_JAR"

        if [ {params.layout} == "single" ]
        then
            #== Run Kneaddata
            python --version
                {params.kneaddata} \
                    --unpaired {params.R1} \
                    -db {params.bowtie_ref} \
                    -t {threads} \
                    --bypass-trf \
                    --remove-intermediate-output \
                    --output {kneaddata}/{params.sample} \
                    --cat-final-output \
                    --output-prefix {params.sample} \
                    --trimmomatic "$TRIM_DIR"

            #----- gzip kneaddata outputs
            pigz {kneaddata}{params.sample}/{params.sample}.fastq
            mv {kneaddata}{params.sample}/{params.sample}.fastq.gz {kneaddata}
            mv {kneaddata}{params.sample}/{params.sample}*contam*.fastq {tmp}

        else

            #== Run Kneaddata
            python --version
                {params.kneaddata} \
                    --input1 {params.R1} \
                    --input2 {params.R2} \
                    -db {params.bowtie_ref} \
                    -t {threads} \
                    --bypass-trf \
                    --remove-intermediate-output \
                    --output {kneaddata}/{params.sample} \
                    --cat-final-output \
                    --output-prefix {params.sample} \
                    --trimmomatic "$TRIM_DIR"

            #----- gzip kneaddata outputs
            pigz {kneaddata}{params.sample}/{params.sample}.fastq
            mv {kneaddata}{params.sample}/{params.sample}.fastq.gz {kneaddata}
            mv {kneaddata}{params.sample}/{params.sample}*contam*.fastq {tmp}
        fi

        #----- Drop the jar symlink directory
        rm -rf "$TRIM_DIR"

"""

#----- Rule to plot kneadData diagnostics
rule plot_filtering:
    input:
        expand(kneaddata + "{sample}.fastq.gz", sample = sample_list)    
    output:
        plots + "KneadData_Diagnostics.png"
    conda: 
        "env_config/r-plotting.yaml"
    resources: cpus="40", maxtime="20:00:00", mem_mb="60gb"
    shell: """
    
        #----- Run plotting script
        # kneaddata and plots dirs are passed explicitly, they hang off
        # out_dir and are not relative to where snakemake runs from
        Rscript scripts/plot_kneadData.R {kneaddata} {plots}

    """

#----- Rule to assign taxonomy
rule taxonomy:
    input:
        profiling_reads
    output: 
        profile = profiles + "{sample}-profiled_metagenome.txt",
        sam = sams + "{sample}.sam.bz2",
        bt2out = bt2logs + "{sample}.bowtie2out.txt"
    params:
        metaphlan_path = config["metaphlan_path"],
        metaphlan_db = config["metaphlan_db"],
        metaphlan_analysis_type = config["metaphlan_analysis_type"],
        mpIndex = config["mpIndex"]
    conda:
        "env_config/metaphlan.yaml",
    resources: cpus="40", maxtime="20:00:00", mem_mb="60gb",
    threads: 40
    shell: """
        {params.metaphlan_path} \
            {input[0]} \
            --input_type fastq \
            -o {output.profile} \
            --samout {output.sam} \
            --bowtie2out {output.bt2out} \
            --bowtie2db {params.metaphlan_db} \
            -x {params.mpIndex} \
            -t {params.metaphlan_analysis_type} \
            --offline \
            --nproc {threads}
      
"""

#----- Rule to merge abundance tables
rule merge_abundances:
    input: 
        expand(profiles + "{sample}-profiled_metagenome.txt", sample = sample_list)
    output:
        mpaMerged + "merged_abundance_table.txt"
    conda: 
        "env_config/metaphlan.yaml"
    resources: cpus="40", maxtime="20:00:00", mem_mb="60gb"
    shell: """
    
        #----- Merge metaphlan tables
        merge_metaphlan_tables.py {input} > {output}
    
    """

#----- Rule to split the merged table into per-rank abundance tables
# Ranks come from the wildcard, so adding e.g. "phylum" is a one word change:
# add it to the constraint below and to all_outputs.
rule extract_rank:
    input:
        table = mpaMerged + "merged_abundance_table.txt",
        samples = config["sample_csv"]
    output:
        mpaMerged + "merged_{rank}_abundance_table.txt"
    wildcard_constraints:
        rank = "kingdom|phylum|class|order|family|genus|species|sgb"
    resources: cpus="1", maxtime="1:00:00", mem_mb="8gb"
    shell: """

        #----- Pull one rank out of the merged table, keeping sample names
        python3 scripts/stratify_abundance.py \
            --input {input.table} \
            --output {output} \
            --rank {wildcards.rank} \
            --samples {input.samples}

    """

#----- Rule to correct batch effects with MMUPHin
rule batch_correct:
    input:
        abd  = mpaMerged + "merged_{rank}_abundance_table.txt",
        meta = config["sample_csv"]
    output:
        mmuphinDir + "{rank}_mmuphin_corrected_pct.txt"
    params:
        batch      = config["mmuphin_batch"],
        covariates = mmuphin_covariates,
        # a lambda rather than a plain string, so the rank wildcard is resolved
        # by snakemake rather than relying on params formatting
        prefix     = lambda wildcards: mmuphinDir + wildcards.rank
    wildcard_constraints:
        rank = "kingdom|phylum|class|order|family|genus|species|sgb"
    conda:
        "env_config/mmuphin.yaml"
    resources: cpus="1", maxtime="4:00:00", mem_mb="16gb"
    shell: """

        #----- adjust_batch is single threaded, hence one cpu.
        # covariates comes last and is variadic, so an empty list expands to
        # nothing and the four required arguments keep their positions. No
        # quoting needed.
        Rscript scripts/run_mmuphin.R \
            {input.abd} {input.meta} \
            {params.batch} {params.prefix} \
            {params.covariates}

    """

#----- Rule to assign function
rule function:
    input: 
        profiling_reads,
        profiles + "{sample}-profiled_metagenome.txt"  
    output: 
        humannProf + "{sample}-humann.log",
        humannProf + "{sample}_genefamilies.tsv",
        humannProf + "{sample}_pathcoverage.tsv",
        humannProf + "{sample}_pathabundance.tsv",
        humannBt2 + "{sample}_bowtie2_aligned.tsv"
    params:
        humann_path = config["humann_path"],
        sample = lambda wildcards: wildcards.sample,
        mode = config["mode"],
        nt_db = config["nt_db"],
        aa_db = config["aa_db"]
    conda:
        "env_config/humann3.yaml",
    resources: cpus="40", maxtime="25:00:00", mem_mb="60gb",
    shell: """
		{params.humann_path} -i {input[0]} \
			-o {humannProf} \
            --output-basename {params.sample} \
			--threads 16 --o-log {output[0]} \
			--taxonomic-profile {input[1]} \
			--search-mode {params.mode} \
			--nucleotide-database {params.nt_db} \
			--protein-database {params.aa_db} &&
        
        #----- Move the bowtie2 alignments out of humann's temp folder.
        mkdir -p {humannBt2} &&
        mv {humannProf}{params.sample}_humann_temp/{params.sample}_bowtie2_aligned.tsv {humannBt2}{params.sample}_bowtie2_aligned.tsv &&

        #----- Remove temp folders
        rm -r {humannProf}{params.sample}_humann_temp
   """

#----- Rule to join functional tables
rule merge_function:
    input:
        expand(humannProf + "{sample}_genefamilies.tsv", sample = sample_list),
        expand(humannProf + "{sample}_pathabundance.tsv", sample = sample_list),
        expand(humannProf + "{sample}_pathcoverage.tsv", sample = sample_list)
    output:
        gfs = humannMerged + "merged_genefamilies.txt",
        pas = humannMerged + "merged_pathabundances.txt",
        pcs = humannMerged + "merged_pathcoverages.txt"
    conda:
        "env_config/humann3.yaml"
    resources: cpus="40", maxtime="25:00:00", mem_mb="60gb",
    shell: """
    
        #----- Merge humann tables
        humann_join_tables --input {humannProf} --output {output.gfs} --file_name genefamilies &&
        humann_join_tables --input {humannProf} --output {output.pas} --file_name pathabundance &&
        humann_join_tables --input {humannProf} --output {output.pcs} --file_name pathcoverage

    """

#----- Rule to normalize tables
rule normalize_function:
    input:
        merged_gfs = humannMerged + "merged_genefamilies.txt",
        merged_pas = humannMerged + "merged_pathabundances.txt",
        merged_pcs = humannMerged + "merged_pathcoverages.txt"
    output:
        cpm_gfs = humannMerged + "merged_genefamilies_cpm.txt",
        cpm_pas = humannMerged + "merged_pathabundances_cpm.txt",
        cpm_pcs = humannMerged + "merged_pathcoverages_cpm.txt"
    params:
        normMethod = config["normMethod"]
    conda: 
        "env_config/humann3.yaml"
    resources: cpus="40", maxtime="25:00:00", mem_mb="60gb",
    shell: """
    
        #----- Run normalization
        humann_renorm_table --input {input.merged_gfs} --units {params.normMethod} --output {output.cpm_gfs} &&
        humann_renorm_table --input {input.merged_pas} --units {params.normMethod} --output {output.cpm_pas} &&
        humann_renorm_table --input {input.merged_pcs} --units {params.normMethod} --output {output.cpm_pcs}
    
    """
    
#----- Rule to rename tables
rule rename_function:
    input:
        cpm_gfs = humannMerged + "merged_genefamilies_cpm.txt",
        cpm_pas = humannMerged + "merged_pathabundances_cpm.txt",
        cpm_pcs = humannMerged + "merged_pathcoverages_cpm.txt"
    output:
        rn_gfs = humannMerged + "merged_genefamilies_cpm_named.txt",
        rn_pas = humannMerged + "merged_pathabundances_cpm_named.txt",
        rn_pcs = humannMerged + "merged_pathcoverages_cpm_named.txt"
    params:
        names = config["mode"],
        utility_mapping = config["utility_mapping"]
    conda: 
        "env_config/humann3.yaml"
    resources: cpus="40", maxtime="25:00:00", mem_mb="60gb",
    shell: """

        #----- Point humann at utility mapping folder
        humann_config --update database_folders utility_mapping {params.utility_mapping} &&

        #----- Run renaming
        humann_rename_table --input {input.cpm_gfs} --names {params.names} --output {output.rn_gfs} &&
        humann_rename_table --input {input.cpm_pas} --names {params.names} --output {output.rn_pas} &&
        humann_rename_table --input {input.cpm_pcs} --names {params.names} --output {output.rn_pcs}
    
    """

#----- Rule to infer taxonomy from function
rule infer_taxonomy_from_function:
    input:
        rn_gfs = humannMerged + "merged_genefamilies_cpm_named.txt",
    output:
        taxProfiles = humannMerged + "merged_taxonomic_profiles.tsv"
    params:
        taxons = config["taxons"],
        utility_mapping = config["utility_mapping"]
    conda: 
        "env_config/humann3.yaml"
    resources: cpus="40", maxtime="25:00:00", mem_mb="60gb",
    shell: """

        #----- Update utility mapping path
        humann_config --update database_folders utility_mapping {params.utility_mapping} &&

        #----- Infer taxonomy
        humann_infer_taxonomy -i {input.rn_gfs} -d {params.taxons} > {output.taxProfiles}
    
    """