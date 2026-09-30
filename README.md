# scTAPAS

**scTAPAS** (single-cell **T**ranscriptome **A**llelic **P**rediction and **A**ssociation **S**tudies) is a set of [Nextflow](https://www.nextflow.io/) pipelines for imputing genotypes (and preparing MHC-region variants for downstream HLA imputation) from donor-wise hg38-aligned single-cell RNA-seq BAM files using [QUILT2](https://github.com/rwdavies/QUILT).

It supports:

- **Genotype imputation** from low-coverage single-cell BAMs against a 1000 Genomes reference panel
- **MHC region extraction** for HLA imputation (e.g., Michigan Imputation Server)
- **Genetic ancestry inference** by projection onto 1000 Genomes principal components
- **cis-eQTL mapping** with [tensorQTL](https://github.com/broadinstitute/tensorqtl)
- **TCR V-gene usage QTL mapping** against imputed HLA (standalone pipeline)

> **Scope note:** scTAPAS imputes SNVs. It does **not** directly impute HLA alleles.

---

## Pipeline entry points

The main entry pipelines differ by the genomic interval they impute:

| Entry pipeline    | Scope                                      | Chunked? | Compatible downstream stages (optional) |
|-------------------|--------------------------------------------|----------|-----------------------------------------|
| `scTAPAS_MHC`     | MHC region on chr6 (default 28–34 Mb)      | No       | MHC extraction                          |
| `scTAPAS_REGION`  | Arbitrary `chrom:region_from-region_to`    | No       | —                                       |
| `scTAPAS_CHROM`   | One whole chromosome                       | Yes      | MHC (if chromosome 6) / ancestry / eQTL         |
| `scTAPAS_CHROMS`  | User-supplied list of chromosomes          | Yes      | MHC (if chromosome 6) / ancestry / eQTL         |
| `scTAPAS_AUTO`    | All autosomes (chr1–chr22)                 | Yes      | MHC / ancestry / eQTL                   |

---

## Installation

Clone the repository (typically on HPC):

```bash
git clone git@github.com:yang-luo-lab/scTAPAS.git
cd scTAPAS
```

### Required software

- **Nextflow** (>= 24.04.2)
- **QUILT2** in a dedicated conda environment (`quilt2_conda_env`)
- A **scTAPAS analysis environment**, either:
  - a conda environment (`conda_env`) containing at least:
    - `bcftools`
    - `plink2`
    - optionally, to run Ancestry inference, cis-eQTL mapping, or TCR V-gene usage QTLs: full environment defined by sctapas_env.yml
  - or HPC modules (`BCFtools`, `PLINK/2`) for the main pipelines and `prep_mhc_imp` (see `nextflow.config`)

### Example conda setup

```bash
# QUILT2 environment
conda create -c conda-forge -c bioconda -n quilt2 'r-quilt>=2.0.3'

# scTAPAS analysis environment
conda env create -f sctapas_env.yml
```

---

## Reference data

For each chromosome you intend to impute, you need:

1. **1000 Genomes WGS phased reference VCF** (hg38), e.g.  
   `CCDG_14151_B01_GRM_WGS_2020-08-05_chr{i}.filtered.shapeit2-duohmm-phased.vcf.gz`
2. **QUILT2 genetic map** (hg38) from the QUILT repository (`maps/hg38`)
3. **Reference genome FASTA** (hg38), with chromosome naming matching BAMs (`6` vs `chr6`)
4. For ancestry inference:
   - **1000 Genomes PLINK2 files** (`.pgen/.pvar/.psam`)
   - relatedness-exclusion sample list (e.g. `deg2_hg38.kin0.related`)

---

## Quick start

scTAPAS is designed for HPC (default executor: **SLURM**).

1. Start from an example params file in `examples/`
2. Edit paths and pipeline options
3. Launch with `nextflow run ...`

Example (`scTAPAS_CHROM`):

```bash
nextflow run workflows/sctapas_quilt2.nf \
  -entry scTAPAS_CHROM \
  -params-file examples/params.scTAPAS_CHROM.yaml \
  -resume
```

Notes:

- `-entry` selects one of: `scTAPAS_MHC`, `scTAPAS_REGION`, `scTAPAS_CHROM`, `scTAPAS_CHROMS`, `scTAPAS_AUTO`
- `-params-file` defines scTAPAS parameters
- Nextflow flags (`-entry`, `-resume`, `-profile`, `-C`, etc.) are command-line flags
- scTAPAS parameters belong in the YAML params file

---

## Common parameters

See `nextflow.config` for the full commented parameter set.

| Parameter             | Description |
|-----------------------|-------------|
| `outdir`              | Parent directory for all outputs |
| `bam_pattern`         | Glob for per-sample BAMs; may include `CHROM` if BAMs are split by chromosome |
| `genetic_map_pattern` | Path pattern to QUILT2 genetic map (`CHROM` placeholder for multi-chrom runs) |
| `thouG_vcf_pattern`   | Path pattern to 1000G phased VCF (`CHROM` placeholder for multi-chrom runs) |
| `ref_fasta`           | hg38 reference FASTA (QUILT2 + variant normalization) |
| `num_chrs`            | `true` if BAM contigs are numeric (`5`) instead of `chr5` (default `false`) |
| `conda_env`           | Path to scTAPAS analysis conda env (if unset, module-based setup is used) |
| `quilt2_conda_env`    | Path to QUILT2 conda env |
| `sample_name_awk`     | AWK snippet to derive sample IDs from BAM paths |

> **Note on `CHROM`:** reference file placeholders (`genetic_map_pattern`, `thouG_vcf_pattern`) are always resolved as `chrN` (`chr6`, `chr22`, etc.), independent of `num_chrs`.  
> `num_chrs` affects only BAM contig naming.

### Optional downstream stages (`scTAPAS_CHROM`, `scTAPAS_CHROMS`, `scTAPAS_AUTO`)

| Flag             | Enables                           | Additional required params |
|------------------|-----------------------------------|----------------------------|
| `prep_mhc_imp`   | MHC extraction for HLA imputation | `mhc_ref_panel_vars` |
| `infer_ancestry` | Ancestry inference                | `onek_gen_patt`, `onek_rel_sams` (optional: `pcs_to_use`, `mhc_from`, `mhc_to`) |
| `run_tensorqtl`  | cis-eQTL mapping                  | `identifiers`, `pseudobulk_pattern`, `covariates_pattern` (optional: `qval_lambda`, `n_pcs_tqtl`, `insert_position_pcs_tqtl`) |

---

## Important note: HLA imputation workflow

scTAPAS does not directly output HLA alleles.  
For HLA imputation from extracted MHC variants, one option is the Michigan Imputation Server.

To reproduce the manuscript setup, use:

1. Reference panel: **Four-digit Multi-ethnic HLA reference panel v2**
2. Array build: **GRCh38**
3. Phasing engine: **EAGLE2**

---

## Pipelines

### `scTAPAS_MHC`

Imputes the chr6 MHC region (default `mhc_from=28000000`, `mhc_to=34000000`) as a single QUILT2 region.  
If `prep_mhc_imp=true`, extracts `mhc_ref_panel_vars` for downstream HLA imputation.

```bash
nextflow run workflows/sctapas_quilt2.nf \
  -entry scTAPAS_MHC \
  -params-file examples/params.scTAPAS_MHC.yaml \
  -resume
```

### `scTAPAS_REGION`

Imputes one arbitrary region as a single QUILT2 interval.

> We recommend not exceeding **~6 Mbp** in a single `scTAPAS_REGION` run.

```bash
nextflow run workflows/sctapas_quilt2.nf \
  -entry scTAPAS_REGION \
  -params-file examples/params.scTAPAS_REGION.yaml \
  -resume
```

### `scTAPAS_CHROM`

Imputes one whole chromosome. QUILT2 chunks are computed (`quilt_chunk_map`), imputed independently, then merged and QC-filtered into one VCF.

```bash
nextflow run workflows/sctapas_quilt2.nf \
  -entry scTAPAS_CHROM \
  -params-file examples/params.scTAPAS_CHROM.yaml \
  -resume
```

### `scTAPAS_CHROMS`

Same as `scTAPAS_CHROM`, iterating over a YAML list `chroms`.

```bash
nextflow run workflows/sctapas_quilt2.nf \
  -entry scTAPAS_CHROMS \
  -params-file examples/params.scTAPAS_CHROMS.yaml \
  -resume
```

### `scTAPAS_AUTO`

Runs chromosome mode across all autosomes (chr1–chr22), no explicit `chroms` list needed.

```bash
nextflow run workflows/sctapas_quilt2.nf \
  -entry scTAPAS_AUTO \
  -params-file examples/params.scTAPAS_AUTO.yaml \
  -resume
```

---

## Outputs

All outputs are written under `outdir`:

| Path                                               | Produced by                   | Contents |
|----------------------------------------------------|-------------------------------|----------|
| `quilt2_output_raw/`                               | all pipelines                 | Raw per-chunk QUILT2 VCFs |
| `chunks_<chrom>.csv`                               | chunked pipelines             | QUILT2 chunk boundaries |
| `sctapas_output/*_info0.8_maf0.05_hwe1e-6.vcf.gz`  | all pipelines                 | Final QC-filtered imputed genotypes (`.tbi` index included) |
| `sctapas_mhc/*_mhc_ref_vars.vcf.gz`                | `prep_mhc_imp=true`           | MHC variants for HLA imputation |
| `sctapas_ancestry/sctapas_pop_assign.tsv` + `*.png`| `infer_ancestry=true`         | Sample-level ancestry assignment and PCA plots |
| `sctapas_tensorqtl/<identifier>/`                  | `run_tensorqtl=true`          | tensorQTL cis outputs |

QC retains variants with:

- `INFO_SCORE > 0.8`
- `0.05 < EAF < 0.95`
- `HWE > 1e-6`

Nextflow creates a `work/` directory. After finalizing analysis, you may remove it.  
Use `-resume` to reuse intermediate results from `work/`.

---

## Standalone pipelines

These run independently from an existing scTAPAS-imputed VCF (or imputed HLA input), without re-running imputation.

### Ancestry inference

```bash
nextflow run workflows/sctapas_ancestry.nf -params-file <params.yaml> -resume
```

Requires: `vcf`, `onek_gen_patt`, `onek_rel_sams`  
Optional: `pcs_to_use`, `mhc_from`, `mhc_to`

### cis-eQTL mapping

```bash
nextflow run workflows/sctapas_tensorqtl.nf -params-file <params.yaml> -resume
```

Requires: `vcf`, `identifiers`, `pseudobulk_pattern`, `covariates_pattern`

### MHC extraction

```bash
nextflow run workflows/sctapas_hla.nf -params-file <params.yaml> -resume
```

Requires: `vcf`, `mhc_ref_panel_vars`

### TCR V-gene usage QTLs

```bash
nextflow run workflows/sctapas_tcrqtl.nf -params-file <params.yaml> -resume
```

Requires: `hla_vcf`, `tcr_path`, `tcr_covs`, `anndata`, `subset_options`,  
`sample_id_col_name`, `clone_id_col_name`, plus TCR column-name parameters in `nextflow.config`.

Output: `sctapas_tcrqtl/plink_out/*_HLA_multi_assoc.*`

---

## Resource configuration

QUILT2 is memory-intensive. `QUILT2` retries with increasing memory; tune:

- `quilt2_cpus` (default `8`)
- `quilt2_memory_base_gb` (default `160`)
- `quilt2_memory_step_gb` (default `24` per retry)
- `quilt2_max_retries` (default `3`)
- `default_cpus` (default `4`)
- `default_queue` (default `short`)

Default executor is **SLURM**.  
Use `-profile local` to run on a local machine.

---

## Manuscript

`manuscript/` contains figure-generation and analysis scripts  
(e.g. `figure2.py`–`figure4.py`, `quilt2_eval.py`, `tcrqtl_eval.py`) used for manuscript reproducibility.  
These are not required for routine pipeline operation.

---

## Limitations and assumptions

- Currently designed for **hg38** workflows
- Input BAMs should be donor-wise and consistently indexed
- Contig naming between BAMs and references must be coherent (`chrN` vs `N` behavior described above)
- Compute demand can be high for large cohorts / whole-chromosome runs

---

## Citation

If you use scTAPAS, please cite the associated manuscript once available.

> TODO: add manuscript DOI / citation text when published.

---

## License

MIT license.
