# GBS-Surveillance-Nextflow

[![Nextflow](https://img.shields.io/badge/Nextflow-DSL2-blue)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/Docker-kizitodevbio%2Fstrepto--pipeline-blue)](https://hub.docker.com/r/kizitodevbio/strepto-pipeline)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**A Modular Nextflow Workflow for End-to-End Genomic Surveillance, Capsular Serotyping, and Resistome Monitoring of Group B Streptococcus Toward Vaccine Preparedness in low-resource settings.**

Foundation: [GBS-Genomics-Pipeline](https://github.com/kizito-devbio/GBS-Genomics-Pipeline)  
Extension: CPS serotyping (Kaptive) + AMRFinderPlus + integrated surveillance reporting.

## Scientific purpose

*Streptococcus agalactiae* (Group B Streptococcus, GBS) is a leading cause of neonatal sepsis, meningitis and stillbirth. Maternal GBS vaccines are advancing toward deployment. This pipeline provides a single, reproducible, containerised workflow that:

- accepts raw paired-end reads **or** pre-assembled genomes,
- produces sequence type (ST), capsular serotype, resistome and virulence profiles,
- reconstructs core-genome phylogeny,
- delivers integrated surveillance outputs and figures,

suitable for workstations, laptops (where resources permit), HPC and cloud environments.

## Key scientific distinction

| Concept | Source | Meaning |
|---------|--------|---------|
| **ST** | MLST (`mlst`, scheme `sagalactiae`) | Sequence type / clonal lineage |
| **Serotype** | CPS typing (Kaptive) | Capsular polysaccharide type (Ia, Ib, II–IX) – vaccine antigen |
| **AMR / Resistome** | AMRFinderPlus | Antimicrobial resistance determinants |
| **Virulence** | Abricate / VFDB (preserved) | Virulence factors |

These are independent biological dimensions.

## Workflow

```
INPUT
  ├── Raw paired-end FASTQ → fastp → SPAdes → QUAST ─┐
  └── Assembled FASTA ───────────────────────────────┘
                         │
                         ▼
                  UNIFIED GENOMES
                         │
        ┌────────────────┼────────────────┐
        ▼                ▼                ▼
      MLST          SEROTYPE_CPS       Prokka
        │                │                │
        ▼                ▼       ┌────────┼────────┐
       ST            Serotype    ▼        ▼        ▼
                           AMRFinderPlus  Virulence  Panaroo
                                │         │         │
                                ▼         ▼         ▼
                               AMR       VF     Core genome
                                                    │
                                                    ▼
                                                 IQ-TREE
                                                    │
                                                    ▼
                                         INTEGRATED SURVEILLANCE
                                              REPORT + FIGURES
```

## Quick start

### Prerequisites
- Nextflow ≥ 22.10
- Docker (recommended) or Singularity/Apptainer or Conda

### Installation
```bash
git clone https://github.com/<your-username>/GBS-Surveillance-Nextflow.git
cd GBS-Surveillance-Nextflow
```

### Run with curated genomes
```bash
nextflow run pipeline.nf \
  -profile docker \
  --curated_dir curated_data \
  --outdir results
```

### Run with raw paired-end reads
```bash
nextflow run pipeline.nf \
  -profile docker \
  --raw_dir /path/to/fastq \
  --outdir results
```

### Resume
```bash
nextflow run pipeline.nf -profile docker --curated_dir curated_data --outdir results -resume
```

## Modules

| Module | Status | Tool |
|--------|--------|------|
| QC | Preserved | fastp + FastQC |
| Assembly | Preserved | SPAdes |
| Quality assessment | Preserved | QUAST |
| Prokka | Preserved | Prokka |
| MLST | Preserved | mlst (sagalactiae) |
| **CPS Serotyping** | **New** | Kaptive (BLAST-based cps locus) |
| **AMRFinderPlus** | **New** | NCBI AMRFinderPlus |
| Virulence | Preserved | Abricate / VFDB |
| Core genome | Preserved | Panaroo |
| Phylogeny | Preserved | IQ-TREE2 |
| Integration / figures | Preserved (to be extended) | existing |

## Serotyping method

**Kaptive** (Priority-1).

- Already present and verified in the foundation Docker image.
- Performs BLAST-based interrogation of the cps locus against curated GBS capsule reference sequences.
- Supports serotypes Ia, Ib, II–IX.
- Output: `results/Serotyping/<sample>_serotype.tsv`
- Ambiguous / no-match results reported explicitly (never forced).

## AMRFinderPlus

- Primary standardised resistome module.
- Existing Abricate/CARD path retained for compatibility.
- Organism flag `Streptococcus_agalactiae` used when supported.
- Database initialised at Docker build time (see `docker/Dockerfile`).
- Requires Docker image **v1.2.0** (AMRFinderPlus added to the foundation image).

## Output structure

```
results/
├── QC/
├── Assembly/
├── Annotation/
├── MLST/
├── Serotyping/          # NEW
├── AMR/                 # Abricate + AMRFinderPlus
├── Virulence/
├── CoreGenome/
├── Phylogeny/
├── Figures/
├── Reports/
└── Logs/
```

## Docker

- Foundation image: `kizitodevbio/strepto-pipeline:v1.1.1`
- Surveillance image: build from `docker/Dockerfile` and tag as `v1.2.0`
- Original `v1.1.1` is left intact.

```bash
cd docker
docker build -t kizitodevbio/strepto-pipeline:v1.2.0 .
```

## Citation

Cite the foundation pipeline and the upstream tools (Nextflow, fastp, SPAdes, QUAST, Prokka, mlst, Kaptive, AMRFinderPlus, Abricate, Panaroo, IQ-TREE, etc.).

## License

MIT License — see [LICENSE](LICENSE).

## Maintainer

Kizito Ibeojo Sylvester-Ali  
GitHub: https://github.com/kizito-devbio
