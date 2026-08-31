# Understanding the Outputs

Everything the pipeline produces is written under the `outDir` set in [`config.yaml`](config.yaml).

- [Directory layout](#directory-layout)
- [Host read filtering](#host-read-filtering)
- [Taxonomic profiling](#taxonomic-profiling)
- [Functional profiling](#functional-profiling)
- [Batch correction](#batch-correction)
- [Plots](#plots)
- [Which table to use](#which-table-to-use)

## Directory layout

```
outDir/
├── kneaddata/
├── metaphlan_taxonomy/
│   ├── profiled_metagenomes/
│   ├── merged_abundances/
│   ├── sams/
│   └── bowtie2_logs/
├── humann_function/
│   ├── profiled_functions/
│   ├── bowtie2_aligned/
│   └── merged_abundances/
├── mmuphin/
└── plots/
```

## Host read filtering

`kneaddata/`

| File | Description |
|------|-------------|
| `<sample>.fastq.gz` | Reads remaining after adapter trimming and removal of host reads. These feed MetaPhlAn and HUMAnN |

Not produced when `skip_kneaddata: True`, in which case the FASTQs named in the sample sheet are profiled directly.

## Taxonomic profiling

`metaphlan_taxonomy/profiled_metagenomes/`

| File | Description |
|------|-------------|
| `<sample>-profiled_metagenome.txt` | Per sample MetaPhlAn profile. Holds every taxonomic rank, and the estimated read counts when `metaphlan_analysis_type` is `rel_ab_w_read_stats` |

`metaphlan_taxonomy/merged_abundances/`

| File | Description |
|------|-------------|
| `merged_abundance_table.txt` | All samples merged into one table. Stacks all eight ranks, so it sums to roughly 100 per rank |
| `merged_species_abundance_table.txt` | Species rows only, with real sample names on the columns |
| `merged_genus_abundance_table.txt` | Genus rows only |
| `named_merged_species_abundance_table.txt` | The species table after the sample name check in `plot_alpha.R` |
| `named_merged_genus_abundance_table.txt` | The genus table after the sample name check |

`metaphlan_taxonomy/sams/` holds the compressed alignments (`<sample>.sam.bz2`), and `metaphlan_taxonomy/bowtie2_logs/` holds the intermediate bowtie2 output (`<sample>.bowtie2out.txt`).

>[!NOTE]
> `merged_abundance_table.txt` carries every rank at once. A phylum row is the same reads as its child species counted again, so it is not a matrix you can hand to a statistical model. Use the per-rank tables for anything downstream.

>[!IMPORTANT]
> A `.bowtie2out.txt` left behind by a killed job will make MetaPhlAn refuse to re-run that sample on every retry. The Snakefile detects this at startup and prints the command to clear them.

**Table format.** Clade names in column 1, one column per sample, values as percentages from 0 to 100. The per-rank tables carry only the leaf name, so a species row reads `s__Escherichia_coli` rather than the full lineage.

## Functional profiling

`humann_function/profiled_functions/`

| File | Description |
|------|-------------|
| `<sample>_genefamilies.tsv` | UniRef gene family abundances in RPK |
| `<sample>_pathabundance.tsv` | MetaCyc pathway abundances |
| `<sample>_pathcoverage.tsv` | Pathway presence and absence calls |
| `<sample>-humann.log` | HUMAnN run log, including the fraction of reads aligned at each stage |

`humann_function/bowtie2_aligned/` holds the nucleotide alignment output (`<sample>_bowtie2_aligned.tsv`).

`humann_function/merged_abundances/`

| File | Description |
|------|-------------|
| `merged_genefamilies.txt` | All samples merged, RPK |
| `merged_pathabundances.txt` | All samples merged |
| `merged_pathcoverages.txt` | All samples merged |
| `merged_genefamilies_cpm.txt` | Normalized using `normMethod` |
| `merged_pathabundances_cpm.txt` | Normalized |
| `merged_pathcoverages_cpm.txt` | Normalized |
| `merged_genefamilies_cpm_named.txt` | Normalized, with human readable names attached |
| `merged_pathabundances_cpm_named.txt` | Normalized and named |
| `merged_pathcoverages_cpm_named.txt` | Normalized and named |
| `merged_taxonomic_profiles.tsv` | Taxonomy inferred from the functional output |

The `_named` tables are the ones to work from. Gene family IDs alone are not readable, and the naming step attaches the descriptions from `utility_mapping`.

## Batch correction

`mmuphin/`, produced only when `run_mmuphin: True`. One set of files per rank listed in `mmuphin_ranks`.

| File | Description |
|------|-------------|
| `<rank>_mmuphin_corrected_pct.txt` | The corrected table. Same shape and same 0 to 100 percent scale as the input |
| `<rank>_mmuphin_pcoa.png` | Bray-Curtis PCoA before and after correction, coloured by batch and shaped by the first covariate |
| `<rank>_mmuphin_diagnostic.pdf` | MMUPHin's own diagnostic plot of the empirical Bayes shrinkage step |

Only the corrected table is a declared Snakemake output. The other two are written alongside it, except when the batch column has a single level, in which case the input is written through unchanged and no figure is produced.

**Reading the PCoA.** The two panels are the same samples before and after correction, and the batch R2 from a PERMANOVA on Bray-Curtis distances is printed in the subtitle. What you want to see is points separating by colour on the left and mixing on the right, with the batch R2 dropping. If the covariate R2 drops with it, the correction is taking biology out along with the batch effect, which usually means batch and group are close to confounded.

**Feature counts.** `adjust_batch` does not drop features. It models only the features that clear its own prevalence filter and carries the rest through, renormalizing each sample afterwards, so the corrected table has exactly the rows the input had. The number actually modelled is reported in the job log:

```
Adjusting for (after filtering) 18 features
```

On a sparse table this can be a small fraction of the total. Features below the filter still move slightly, because each sample is renormalized around the ones that were adjusted.

## Plots

`plots/`

| File | Description |
|------|-------------|
| `KneadData_Diagnostics.png` | Reads surviving each filtering step, per sample |
| `Species_Shannon_Diversity_Per_Sample.png` | Shannon index per sample at species level |
| `Genus_Shannon_Diversity_Per_Sample.png` | Shannon index per sample at genus level |
| `Species_alpha_diversity_for_plotting.csv` | Shannon, Simpson and observed richness per sample, long format |
| `Genus_alpha_diversity_for_plotting.csv` | Same at genus level |

The CSVs hold all three alpha diversity metrics, not just the one plotted, so they are the place to go for Simpson or observed richness.

## Which table to use

| Analysis | Table |
|----------|-------|
| Alpha diversity | `merged_species_abundance_table.txt`, uncorrected |
| Beta diversity, ordination | `species_mmuphin_corrected_pct.txt` when batch correction is on, otherwise the uncorrected species table |
| Differential abundance | Same as beta diversity |
| Functional analysis | `merged_genefamilies_cpm_named.txt` or `merged_pathabundances_cpm_named.txt` |

>[!IMPORTANT]
> Alpha diversity is computed on the raw table by design, and `plot_alpha.R` is deliberately left pointing at the uncorrected table. Batch correction changes per sample composition, which changes richness and evenness in ways that are hard to interpret. Batch correction matters for between sample comparisons, not within sample ones.
