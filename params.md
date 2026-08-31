# Parameters

All pipeline settings live in [`config.yaml`](config.yaml). This page describes every parameter in that file.

- [General settings](#general-settings)
- [Optional rules](#optional-rules)
- [MetaPhlAn](#metaphlan)
- [HUMAnN](#humann)
- [MMUPHin](#mmuphin)
- [Tooling and databases](#tooling-and-databases)
- [Unused parameters](#unused-parameters)

## General settings

| Parameter | Type | Description |
|-----------|------|-------------|
| `sample_csv` | string | Path to the sample sheet. Needs `sample_id` and `fastq_1` always, plus `fastq_2` when `layout` is `paired`. Any other column is treated as sample metadata |
| `layout` | string | Either `single` or `paired`. Controls whether `fastq_2` is required and whether KneadData runs in paired mode |
| `proj_dir` | string | Path to the working directory |
| `outDir` | string | Path to where results are written. Everything the pipeline produces hangs off this directory, so it does not have to sit inside `proj_dir` |
| `tmp_dir` | string | Path to a temporary directory. Created at startup if it does not exist |

## Optional rules

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `skip_kneaddata` | boolean | `False` | Skip host read filtering and profile the FASTQs named in the sample sheet directly |
| `run_mmuphin` | boolean | `False` | Run MMUPHin batch correction. When absent or `False`, `rule batch_correct` is not in the DAG |

`skip_kneaddata: True` with `layout: paired` uses `fastq_1` only, because HUMAnN takes a single input file and there is nowhere to put the second mate without concatenating. The pipeline prints a warning at startup when this combination is set.

## MetaPhlAn

| Parameter | Type | Description |
|-----------|------|-------------|
| `metaphlan_analysis_type` | string | Either `rel_ab_w_read_stats` for estimated counts alongside relative abundance, or `rel_ab` for relative abundance only. See `metaphlan -h` for other accepted values |
| `metaphlan_db` | string | Path to the directory holding the MetaPhlAn bowtie2 database |
| `mpIndex` | string | Name of the database index inside `metaphlan_db`. Currently `mpa_vOct22_CHOCOPhlAnSGB_202212` |

>[!NOTE]
> `merge_metaphlan_tables.py` reads only the clade name and relative abundance columns, so the estimated read counts produced by `rel_ab_w_read_stats` do not reach the merged tables. They have to be read from the per sample profiles in `metaphlan_taxonomy/profiled_metagenomes/`.

## HUMAnN

| Parameter | Type | Description |
|-----------|------|-------------|
| `mode` | string | Protein database to search. Either `uniref90` or `uniref50` |
| `normMethod` | string | Normalization method passed to `humann_renorm_table`. Default is `cpm` for counts per million. `relab` for relative abundance is also accepted |
| `taxons` | string | Taxonomic mode passed to `humann_infer_taxonomy`. Currently `uniref90-tol-lca` |
| `nt_db` | string | Path to the ChocoPhlAn nucleotide database |
| `aa_db` | string | Path to the UniRef protein database |
| `utility_mapping` | string | Path to the HUMAnN utility mapping files, used when renaming gene family and pathway tables |

`mode` and `aa_db` have to agree. Setting `mode: uniref50` while `aa_db` points at a UniRef90 database will not error, it will just search the wrong database.

## MMUPHin

| Parameter | Type | Description |
|-----------|------|-------------|
| `mmuphin_batch` | string | Sample sheet column holding the technical variable being corrected for. Exactly one. MMUPHin does not accept more |
| `mmuphin_covariates` | list | Sample sheet columns holding the biological variables to preserve, for example `["group"]` or `["group", "sex"]`. An empty list `[]` corrects on batch alone |
| `mmuphin_ranks` | list | Taxonomic ranks to correct, for example `["species"]` or `["species", "genus"]`. Accepted values are `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species`, `sgb` |

Column names are case-sensitive and must match the sample sheet exactly. The pipeline checks them at startup and fails with the list of available metadata columns if one is missing:

```
run_mmuphin is on but sample sheet sample_fastq_list_single.csv has no column(s) run_id.
Available metadata columns: batch, group
```

Each rank listed in `mmuphin_ranks` needs its per-rank table built, which `rule extract_rank` does automatically.

**Covariates and degrees of freedom.** Each covariate level costs a design matrix column. On small studies, two or three factor covariates can use up the available degrees of freedom, and MMUPHin will stop with `Covariates are confounded!` even when the covariates are not redundant. Start with one covariate and add more only if the sample size supports it.

**Choosing a rank.** Species is the usual target. Genus is cheap to add by listing both, since the same rule runs once per rank. Each rank is corrected independently, so the results will not reconcile across ranks.

## Tooling and databases

| Parameter | Type | Description |
|-----------|------|-------------|
| `kneaddata_path` | string | Command or path used to invoke KneadData. Default `kneaddata` resolves it from the conda environment |
| `metaphlan_path` | string | Command or path used to invoke MetaPhlAn |
| `humann_path` | string | Command or path used to invoke HUMAnN |
| `bowtie_ref` | string | Path to the host bowtie2 index used by KneadData for read filtering. Switch this to the mouse index for mouse data. A commented mouse path is included in the config |


