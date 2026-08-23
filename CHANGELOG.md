# Changelog

All notable changes to GBS-Genomics-Pipeline are documented here.

Format follows [Keep a Changelog](https://keepachangelog.com/).

## [1.1.0] - 2026-07-20

### Removed
- Human genome decontamination step removed from active workflow (modules frozen in `modules/frozen/`)
- Raw-read pathway now flows directly: QC → Assembly

### Added
- Standardized output directory structure (`QC/`, `Assembly/`, `AMR/`, `Virulence/`, `MLST/`, `CoreGenome/`, `Phylogeny/`, `Figures/`, `Reports/`, `Logs/`)
- Publication-quality visualization stage (`bin/generate_figures.py`) consuming real pipeline outputs
- Skip reports when figure data are unavailable (no placeholder visualizations)
- Parameter validation at workflow start
- Comprehensive logging in all processes (stage, inputs, outputs, versions, elapsed time)
- Open-source documentation: LICENSE, CITATION.cff, CONTRIBUTING.md, CODE_OF_CONDUCT.md
- Module and parameter reference docs in `docs/`
- Dockerfile additions: Panaroo 1.5.0, ete3 3.1.3, seqkit 2.8.2

### Changed
- `INTEGRATION_VISUALIZATION` is now the final workflow stage
- MLST outputs are collected and passed to visualization
- Improved error handling and graceful skip for single-sample phylogeny
- Professional README rewrite with workflow diagram and troubleshooting

### Fixed
- Removed hard-coded email from config defaults
- Corrected publishDir paths across all modules
- Core genome step handles <2 samples without pipeline failure

## [1.0.0] - 2026-01-01

### Added
- Initial release with QC, decontamination, assembly, annotation, AMR, virulence, MLST, core genome, phylogeny
- Docker, Singularity, and Conda profiles
