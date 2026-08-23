# Frozen Modules — Human Genome Decontamination

These modules were **removed from the active workflow** as of the current release.

## Why removed

Human genome decontamination (BWA alignment against GRCh38 followed by read filtering) is no longer required for this bacterial genomics pipeline. Clinical and environmental *Streptococcus agalactiae* (GBS) sequencing runs typically contain negligible human DNA when samples are pure bacterial cultures. Removing this step:

- Eliminates a large reference download (~3 GB GRCh38)
- Reduces per-sample compute time significantly
- Simplifies the raw-read pathway to **QC → Assembly**
- Avoids accidental loss of reads when decontamination filters are overly aggressive

## Frozen files

| File | Original purpose |
|------|------------------|
| `decontam.nf.disabled` | BWA-based removal of human-mapped reads |
| `download_human_ref.nf.disabled` | Download and index GRCh38 human reference |
| `human_ref.yml.disabled` | Conda environment for human reference setup |

## Re-enabling (not recommended)

To restore decontamination, re-include the modules in `pipeline.nf` and re-add the `human_ref_dir` / `human_accession` parameters to `nextflow.config`. See git history for the previous workflow wiring.
