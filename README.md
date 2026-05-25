# ldsc_continuous_covariate

Snakemake workflow for **continuous stratified LD Score Regression (S-LDSC)** using a sample-specific continuous annotation and a continuous covariate annotation.

The intended model is:

```text
baseline LDSC annotations + covariate continuous annotation + sample continuous annotation
```

For iPKM analyses, the covariate can be a cell-type mean iPKM annotation, and the sample annotation can be the iPKM profile of a specific sample or age.

## Repository structure

```text
ldsc_continuous_covariate/
├── workflow/
│   └── Snakefile
├── config/
│   ├── config.yaml
│   └── samples.tsv
├── scripts/
│   └── make_continuous_annot.py
├── README.md
├── environment.yml
└── .gitignore
```

## Workflow overview

```text
sample continuous BED
        │
        ├──> continuous .annot.gz
        │
        ├──> LD scores
        │
        ▼
covariate continuous BED
        │
        ├──> continuous .annot.gz
        │
        ├──> LD scores
        │
        ▼
LDSC-ready .sumstats.gz
        │
        ▼
ldsc.py --h2
  baseline + covariate + sample
        │
        ▼
results/{trait}/{celltype}.{sample}.results
        │
        ▼
results/sample_annotation_results.tsv
```

## Input BED format

Both sample-specific and covariate BED files must contain four columns:

```text
chr    start    end    value
```

Example:

```text
1    713739    714831    11.84087621
1    714876    715568    5.124349145
1    762078    763310    6.197571645
```

For each SNP, the value assigned in the LDSC annotation is the value of the interval overlapping that SNP. SNPs outside intervals receive value `0`.

If a SNP overlaps multiple intervals, the aggregation rule is controlled by:

```yaml
overlap_strategy: "mean"
```

Supported values are `mean`, `max`, `sum`, and `first`.

## Sample table

Edit:

```text
config/samples.tsv
```

Required columns:

```text
sample_id    celltype    sample_bed    covariate_bed
```

Example:

```text
sample_id	celltype	sample_bed	covariate_bed
sample_001	cortical	/path/to/beds/cortical/sample_001.bed	/path/to/beds/covariates/cortical_mean.bed
sample_002	Astro	/path/to/beds/Astro/sample_002.bed	/path/to/beds/covariates/Astro_mean.bed
```

Each row defines one sample-specific LDSC analysis unit. The same covariate BED can be reused across multiple samples from the same cell type.

## Configuration

All paths in `config/config.yaml` are interpreted relative to the repository root when running:

```bash
snakemake --snakefile workflow/Snakefile --cores 4
```

Edit:

```text
config/config.yaml
```

Main fields:

```yaml
samples_table: "config/samples.tsv"
plink_prefix_template: "/path/to/plink_reference/{chr}"
ldsc_py: "/path/to/ldsc/ldsc.py"
baseline_prefix: "/path/to/baseline/baseline."
weights_prefix: "/path/to/weights/"
frq_prefix: "/path/to/frequency/"
sumstats_dir: "/path/to/sumstats"
traits: []
```

If `traits: []`, the workflow automatically uses all files in `sumstats_dir` matching `*.sumstats.gz`.

## LDSC model

For each trait, cell type, and sample, the workflow runs:

```bash
ldsc.py \
  --h2 trait.sumstats.gz \
  --ref-ld-chr baseline.,covariate.,sample. \
  --w-ld-chr weights. \
  --frqfile-chr frequency. \
  --overlap-annot \
  --print-coefficients \
  --out results/trait/celltype.sample
```

The key term is the sample-specific annotation, interpreted conditional on both the standard LDSC baseline annotations and the continuous covariate annotation.

For downstream plots across sample age, the most useful statistic is usually `Coefficient_z-score` for the row corresponding to the sample-specific annotation.

The workflow writes these rows to:

```text
results/sample_annotation_results.tsv
```

## Install environment

Create the mamba/conda environment:

```bash
mamba env create -f environment.yml
mamba activate ldsc_continuous_covariate
```

Install LDSC separately or point the config to an existing installation:

```yaml
ldsc_py: "/path/to/ldsc/ldsc.py"
```

## Dry run

```bash
snakemake --snakefile workflow/Snakefile -n -p
```

## Run

```bash
snakemake --snakefile workflow/Snakefile --cores 4
```

Increase `--cores` depending on the available machine.

## Main outputs

```text
annotations/covariates/{celltype}/{chr}.annot.gz
annotations/samples/{sample}/{chr}.annot.gz

ldscores/covariates/{celltype}/{chr}.l2.ldscore.gz
ldscores/samples/{sample}/{chr}.l2.ldscore.gz

results/{trait}/{celltype}.{sample}.results
results/all_results.tsv
results/sample_annotation_results.tsv
```
