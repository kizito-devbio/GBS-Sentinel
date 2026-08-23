# Parameter Reference

All parameters are defined in `nextflow.config` and can be overridden on the command line.

## Input Parameters

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `curated_dir` | path | `false` | One of curated/raw | Directory containing assembled genome FASTA files (`.fa`, `.fna`, `.fasta`) |
| `raw_dir` | path | `false` | One of curated/raw | Directory containing paired-end FASTQ files named `sample_1.fastq` and `sample_2.fastq` |

## Output Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `outdir` | path | `./results` | Root output directory for all pipeline results |

## Reference Data

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `blast_db_dir` | path | `${projectDir}/references/streptococcus_genus_db` | Local cache for Streptococcus BLAST database |
| `background_ref_url` | URL | NCBI GCF_000009585.1 | Reference genome URL for background selection |

## Analysis Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mlst_scheme` | string | `sagalactiae` | PubMLST scheme for S. agalactiae MLST |
| `min_n50` | integer | `10000` | Minimum N50 (bp) for assemblies to pass quality filter |
| `min_length` | integer | `1800000` | Minimum expected genome length (reserved) |
| `max_length` | integer | `2400000` | Maximum expected genome length (reserved) |
| `max_contigs` | integer | `200` | Maximum contig count (reserved) |
| `skip_tbl2asn` | boolean | `true` | Skip NCBI tbl2asn submission step |

## Resource Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `max_cpus` | integer | auto (CPUs − 1) | Maximum CPU cores for compute-intensive processes |
| `max_memory` | string | `6 GB` | Maximum memory for high-compute processes |

## Utility

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `help` | boolean | `false` | Print parameter reference and exit |

## Examples

```bash
# Curated genomes with custom MLST scheme
nextflow run pipeline.nf -profile docker \
  --curated_dir genomes/ --mlst_scheme sagalactiae --outdir results

# Raw reads with relaxed N50 filter
nextflow run pipeline.nf -profile docker \
  --raw_dir reads/ --min_n50 5000 --outdir results

# HPC with Singularity
nextflow run pipeline.nf -profile singularity,cluster \
  --curated_dir genomes/ --outdir /scratch/user/gbs_results
```

## Deprecated Parameters

The following parameters were removed when human decontamination was frozen:

| Parameter | Former purpose |
|-----------|----------------|
| `human_ref_dir` | Local GRCh38 reference cache |
| `human_accession` | NCBI accession for human reference |

See `modules/frozen/README.md` for details.
