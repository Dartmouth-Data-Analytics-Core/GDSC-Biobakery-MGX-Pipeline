# Dartmouth GDSC Biobakery Pipeline
<img src="img/cqb_logo.jpg" alt="CQB Logo" width="200" align="right"/>

![Version](https://img.shields.io/badge/version-1.0-blue)

The GDSC Biobakery (metagenomics) pipeline provides preprocessing, taxonomic profiling and functional profiling of shotgun metagenomic sequencing data, implemented through [Snakemake](https://snakemake.readthedocs.io/en/stable/) for use on the [Dartmouth Discovery HPC](https://rc.dartmouth.edu/discoveryhpc/). The pipeline is built on [The Biobakery](https://github.com/biobakery) toolset and supports both single- and paired-end libraries. Optional batch effect correction across sequencing runs, plates or studies is available through [MMUPHin](https://huttenhower.sph.harvard.edu/mmuphin/). Software dependencies are installed per rule by Snakemake from the conda environment files in [`env_config/`](env_config/). 

## Documentation

- [Summary](#summary)
- [Quick start](#quick-start)
- [Installation](#installation)
- [Configuration](#configuration)
- [Optional Features](#optional-features)
- [Parameters](params.md)
- [Understanding the Outputs](understanding_outputs.md)
- [Contact](#contact)

## Summary

**Currently the pipeline performs the following:**

- (**Optional**) Host read filtering and adapter trimming using [KneadData v0.12.0](https://huttenhower.sph.harvard.edu/kneaddata/) and [Trimmomatic v 0.39](https://github.com/usadellab/trimmomatic$0)
- Taxonomic profiling using [MetaPhlAn v4.0.6](https://huttenhower.sph.harvard.edu/metaphlan/) against the CHOCOPhlAn SGB database
- Merging of per sample profiles into a single abundance table
- Splitting of the merged table into per-taxonomic rank abundance tables
- Functional profiling using [HUMAnN v3.7](https://huttenhower.sph.harvard.edu/humann/) against ChocoPhlAn and UniRef
- CPM normalization and renaming of gene family and pathway tables
- Inference of taxonomic profiles from the functional output
- Alpha diversity calculation and plotting using [vegan](https://cran.r-project.org/package=vegan)
- Read filtering diagnostics from the KneadData logs
- **(optional)** Batch effect correction using [MMUPHin v1.16.0](https://huttenhower.sph.harvard.edu/mmuphin/)


## Quick Start
1. Populate [`sample_fastq_list_single.csv`](sample_fastq_list_single.csv) or [`sample_fastq_list_paired.csv`](sample_fastq_list_paired.csv) with your sample information
2. Set the paths and options in [`config.yaml`](config.yaml) (see [Parameters](params.md))

3. Submit [`job_script.sh`](job_script.sh) to the SLURM scheduler

## Installation

Clone the repository:

```shell
git clone https://github.com/Dartmouth-Data-Analytics-Core/GDSC-MGX-Pipieline
cd GDSC-MGX-Pipieline
```

## Configuration

**1. Sample sheet**

Populate [`sample_fastq_list_single.csv`](sample_fastq_list_single.csv) or [`sample_fastq_list_paired.csv`](sample_fastq_list_paired.csv) with your sample information. This is a comma-separated file with the following columns:

| Column | Description |
|--------|-------------|
| `sample_id` | Short sample identifier used to name all output files |
| `fastq_1` | Path to the R1 FASTQ file |
| `fastq_2` | Path to the R2 FASTQ file. Required only when `layout` is `paired` |
| `batch` | Optional. Technical variable for MMUPHin, such as sequencing run or plate |
| `group` | Optional. Biological variable for MMUPHin, such as condition or treatment |

>[!IMPORTANT] Any column that is not `sample_id`, `fastq_1` or `fastq_2` is treated as sample metadata and is available to MMUPHin. The columns can be named anything as long as the names in `config.yaml` match.

>[!IMPORTANT] `sample_id` becomes part of every output filename, so it cannot contain whitespace or slashes. The Snakefile checks this along with duplicate and blank sample IDs, and reports every problem it finds in one pass.

An example of a bare-bones sample sheet for single-end data is shown below. This configuration implied MMUPHin will not be ran. 

```
sample_id,fastq_1
S1743,/path/to/CFKC00252.fastq.gz
S1916,/path/to/CFKC00264.fastq.gz
```



**2. Pipeline parameters**

All settings live in [`config.yaml`](config.yaml). At minimum, set the following before the first run:

| Parameter | Description |
|-----------|-------------|
| `sample_csv` | Path to your sample sheet |
| `layout` | Either `single` or `paired` |
| `proj_dir` | Path to the working directory |
| `outDir` | Path to where results are written |
| `tmp_dir` | Path to a temporary directory |

The reference database paths under *Tooling and databases* point at the shared references on `/dartfs/rc/nosnapshots/` and do not normally need changing. See [Parameters](params.md) for a full description of every parameter.

**3. Job submission script**

[`job_script.sh`](job_script.sh) activates the shared Snakemake environment and runs the workflow against [`cluster_profile/`](cluster_profile/).

>[!IMPORTANT]
> Change the `--mail-user` line in [`job_script.sh`](job_script.sh) to your own address before submitting.

Conda environments are built on first use. A prefix for where to install them can be set in the Snakemake call through the `--conda-prefix` argument. Building only occurs once, and concurrent runs can point at that path for ready usage.


**4. Submitting the job**

```shell
sbatch job_script.sh
```

Snakemake stays running as the parent job and submits one SLURM job per rule instance, up to 10 at a time. The parent log is written to `MGX_<jobid>.log`, and per rule logs go to `slurm_logs/log_<rule>_<jobid>.log`.

>[!NOTE]
> MetaPhlAn refuses to start if its bowtie2out file already exists, so a job killed by a timeout or node failure leaves a file behind that makes every retry fail the same way. The Snakefile detects this at startup and prints the affected samples along with the command to clear them.

For an explanation of all output files, see [Understanding the Outputs](understanding_outputs.md).

## Optional Features

### Skipping host read filtering

>[!IMPORTANT]
> Only enable this when your FASTQs have already been quality trimmed and depleted of host reads. Set `skip_kneaddata: true` in your config to profile them as they are.

```yaml
skip_kneaddata: true
```

When enabled, `rule filter` never runs and MetaPhlAn and HUMAnN read the FASTQs named in the sample sheet directly. For paired data this uses `fastq_1` only, since HUMAnN takes a single input file, and the pipeline prints a warning saying so at startup.

### MMUPHin batch correction

[MMUPHin](https://huttenhower.sph.harvard.edu/mmuphin/) removes technical batch effects from abundance tables while holding onto the biological signal you tell it to keep. `adjust_batch` fits a location and scale model per feature per batch, shrinks the estimates with an empirical Bayes step, and returns an adjusted table on the same scale as the input.

>[!IMPORTANT]
> MMUPHin requires a `batch` column in the sample sheet with at least two levels and at least two samples per level. Set `run_mmuphin: true` in your config to enable it. With `run_mmuphin: false` the rule is not in the DAG at all.

An example sample sheet header is shown below configured to match the config.

```
sample_id,batch,group,fastq_1

```

Appropriate entries in **`config.yaml`** informs the pipeline that the batch variable is under the column name `batch` and the covariate you want included is under the column `group`. 

>[!IMPORTANT]
> Any number of covariates can be included, so long as they have an associated column in the sample sheet and are listed as such in the config as a comma separated list, with each list element in quotes. 

```yaml
run_mmuphin: true
mmuphin_batch: "batch"
mmuphin_covariates: ["group"]
mmuphin_ranks: ["species"]
```

| Parameter | Description |
|-----------|-------------|
| `run_mmuphin` | Boolean. `false` keeps `rule batch_correct` out of the workflow entirely |
| `mmuphin_batch` | Sample sheet column holding the technical variable. Exactly one, MMUPHin does not accept more |
| `mmuphin_covariates` | List of sample sheet columns holding the biological variables to preserve. An empty list `[]` corrects on batch alone |
| `mmuphin_ranks` | List of taxonomic ranks to correct. Valid values are `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species`, `sgb` |

The rule runs on the per-rank tables from `rule extract_rank`, not on `merged_abundance_table.txt`. MetaPhlAn writes one row per taxonomic level with the full lineage, so the merged table stacks all eight ranks into one file. It sums to roughly 100 per rank, and a phylum row is the same reads as its child species counted a second time. `adjust_batch` assumes each sample column is a composition summing to 1 and that features are independent, and the stacked table breaks both. The pipeline checks the column sums and stops with an explanation if it is handed the wrong table.

Each rank is corrected on its own, so genus values will not sum to their species-level children after correction. This is expected.

**Batch and group must mean different things.** `batch` is the technical variable being removed, such as sequencing run, plate or extraction date. The covariates are the biological variables being preserved, such as condition, treatment or timepoint. Each covariate needs its own column, since MMUPHin builds a design matrix column per covariate level. Every sample needs a value in both, and blank cells are rejected with the list of samples that are missing one.

>[!WARNING]
> Setting `mmuphin_covariates: []` corrects on batch alone. With nothing declared to hold onto, any biological difference that tracks the batch is removed along with the batch effect. The pipeline warns at startup when this is set, but it cannot detect whether it did any damage.

**Cases that stop the run.** Three designs are checked before anything is fit:

| Case | Behavior |
|------|----------|
| One batch level | Input written through unchanged, job exits successfully. There is nothing to correct |
| A batch with a single sample | Stops. `adjust_batch` cannot estimate parameters from one sample |
| A covariate perfectly confounded with the batch | Stops. The correction would remove the biology along with the batch effect |

A confounded design means every batch contains only one level of a covariate, so the batch effect and that biological effect are the same contrast. MMUPHin does not warn about this, it just returns a table with the signal gone, so the pipeline stops instead:

```
Error: 'group' is perfectly confounded with 'batch': every batch holds only one level of it.
     CF NonCF
  A   9     0
  B   0     9
```

There is no way to correct this design. It has to be handled when samples are assigned to batches.

**Alpha diversity stays on the raw table.** Batch correction matters for beta diversity, ordination and differential abundance. Feed the corrected table into those analyses, not into per sample richness or evenness.

## Contact

**Contact and questions:** Please address questions to *DataAnalyticsCore@groups.dartmouth.edu* or submit an issue in the GitHub repository.

**This pipeline was created with funds from the COBRE grant 1P20GM130454. If you use the pipeline in your own work, please acknowledge the pipeline by citing the grant number in your manuscript.**
