# Module Reference

Each module is a Nextflow DSL2 process file in `modules/`. Frozen (disabled) modules are in `modules/frozen/`.

## Active Modules

### `download_blast_db.nf` — DOWNLOAD_BLAST_DB

Downloads and indexes a Streptococcus agalactiae reference genome for taxonomic BLAST confirmation.

| | |
|---|---|
| **Input** | None (reference downloaded or cached) |
| **Output** | `blast_db/` directory with formatted nucleotide database |
| **Tools** | curl, makeblastdb |
| **Published to** | `references/streptococcus_genus_db/`, `Logs/` |

---

### `qc.nf` — QC

Trims and quality-assesses raw paired-end reads.

| | |
|---|---|
| **Input** | `tuple(sample, [R1, R2])` FASTQ |
| **Output** | Trimmed gzipped FASTQ, FastQC reports |
| **Tools** | fastp, FastQC |
| **Published to** | `QC/<sample>/`, `Logs/` |

---

### `assembly.nf` — ASSEMBLY

De novo genome assembly from trimmed reads.

| | |
|---|---|
| **Input** | `tuple(sample, [R1, R2])` trimmed FASTQ |
| **Output** | `<sample>_assembly.fasta` |
| **Tools** | SPAdes 4.2.0 |
| **Published to** | `Assembly/`, `Logs/` |

---

### `quality_assess.nf` — QUALITY_ASSESS

Assembly quality metrics via QUAST.

| | |
|---|---|
| **Input** | `tuple(sample, assembly.fasta)` |
| **Output** | QUAST TSV, N50/total/contig metric files |
| **Tools** | QUAST |
| **Published to** | `Assembly/`, `Logs/` |

---

### `taxonomy.nf` — TAXONOMY

Confirms species identity via BLASTn against Streptococcus reference.

| | |
|---|---|
| **Input** | `tuple(sample, genome.fasta)`, BLAST database |
| **Output** | `<sample>_taxonomy.txt` |
| **Tools** | BLAST+ |
| **Published to** | `Annotation/`, `Logs/` |

---

### `background_selection.nf` — BACKGROUND_SELECTION

Collects sample genomes and downloads a reference for contextual analysis.

| | |
|---|---|
| **Input** | Collected genome FASTA files |
| **Output** | `background_genomes/*.fna` |
| **Tools** | wget |
| **Published to** | `Reports/`, `Logs/` |

---

### `functional_annotation.nf` — FUNCTIONAL_ANNOTATION

Gene prediction (Prokka) and AMR screening (Abricate/CARD).

| | |
|---|---|
| **Input** | `tuple(sample, genome.fasta)` |
| **Output** | Prokka directory, `<sample>_amr.tsv` |
| **Tools** | Prokka 1.14.6, Abricate |
| **Published to** | `Annotation/<sample>/`, `AMR/`, `Logs/` |

---

### `virulence_factor.nf` — VIRULENCE_FACTOR

Virulence factor detection via Abricate/VFDB.

| | |
|---|---|
| **Input** | `tuple(sample, prokka.fna)` |
| **Output** | `<sample>_vf.tsv` |
| **Tools** | Abricate (VFDB) |
| **Published to** | `Virulence/`, `Logs/` |

---

### `mlst_extraction.nf` — MLST_EXTRACTION

Multi-Locus Sequence Typing assignment.

| | |
|---|---|
| **Input** | `tuple(sample, prokka.fna)` |
| **Output** | `<sample>_mlst.tsv` |
| **Tools** | mlst (PubMLST) |
| **Published to** | `MLST/`, `Logs/` |

---

### `core_genome.nf` — CORE_GENOME

Pangenome analysis and core gene alignment.

| | |
|---|---|
| **Input** | Collected Prokka GFF files (≥2 required) |
| **Output** | `core_alignment.aln`, Panaroo results |
| **Tools** | Panaroo 1.5.0 |
| **Published to** | `CoreGenome/`, `Logs/` |

---

### `phylogeny.nf` — PHYLOGENY

Maximum-likelihood phylogenetic tree construction.

| | |
|---|---|
| **Input** | Core gene alignment |
| **Output** | `gbs_phylogeny_tree.nwk` |
| **Tools** | IQ-TREE 2 |
| **Published to** | `Phylogeny/`, `Logs/` |

---

### `integration_visualization.nf` — INTEGRATION_VISUALIZATION

Final stage: generates publication figures from all analysis outputs.

| | |
|---|---|
| **Input** | Tree, AMR TSVs, VF TSVs, MLST TSVs |
| **Output** | Figures (PNG/SVG/PDF), data tables, skip reports |
| **Tools** | `bin/generate_figures.py` (matplotlib, seaborn, ete3) |
| **Published to** | `Figures/`, `Reports/` |

---

## Frozen Modules (disabled)

Located in `modules/frozen/`. See [frozen/README.md](../modules/frozen/README.md).

| Module | Reason frozen |
|--------|---------------|
| `decontam.nf.disabled` | Human decontamination no longer required |
| `download_human_ref.nf.disabled` | GRCh38 reference download no longer required |

## Unused Module

### `alignment.nf` — ALIGNMENT

MAFFT multiple sequence alignment. Not currently wired into the main workflow (Panaroo produces aligned core genes directly). Retained for potential future use.
