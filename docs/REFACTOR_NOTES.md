# REFACTOR_NOTES — GBS-Surveillance-Nextflow

## Preserved (unchanged)

- Dual-input architecture (raw FASTQ + curated FASTA → unified genome channel)
- modules/qc.nf (fastp + FastQC)
- modules/assembly.nf (SPAdes)
- modules/quality_assess.nf (QUAST + N50 filter)
- modules/download_blast_db.nf
- modules/taxonomy.nf
- modules/background_selection.nf
- modules/functional_annotation.nf (Prokka + Abricate/CARD)
- modules/virulence_factor.nf (Abricate/VFDB)
- modules/mlst_extraction.nf (mlst, scheme sagalactiae)
- modules/core_genome.nf (Panaroo)
- modules/phylogeny.nf (IQ-TREE2)
- modules/integration_visualization.nf
- modules/frozen/
- bin/
- curated_data/ (real GBS genomes A909, NEM316, NGBS128)
- test_data/
- nextflow.config (profiles, resources, resume, error handling)
- LICENSE, CONTRIBUTING.md, CODE_OF_CONDUCT.md, CITATION.cff
- Existing Docker base and tool installations (except addition of AMRFinderPlus)

## Modified

- pipeline.nf
  - Updated header/branding to GBS-Surveillance-Nextflow
  - Added includes for SEROTYPE_CPS and AMRFINDERPLUS
  - Wired SEROTYPE_CPS and AMRFINDERPLUS in parallel after genome channel
  - Updated help text and startup banner
  - All original processes and channel logic preserved

- README.md
  - Rewritten for surveillance / vaccine-preparedness purpose

- docker/Dockerfile
  - Appended section 46: AMRFinderPlus binary + database init + build-time smoke test
  - Target tag: v1.2.0 (does not alter historical v1.1.1)

## Added

- modules/serotyping/serotype_cps.nf  (Kaptive-based CPS serotyping)
- modules/amrfinderplus/amrfinderplus.nf  (AMRFinderPlus resistome)
- docs/REFACTOR_NOTES.md  (this file)

## Removed

- none
