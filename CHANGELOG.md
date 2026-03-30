# Changelog

## [v2.1.0] - 2026-03-30

### New features
- CAFE runs in three separate jobs
- CAFE tests which run was the most likely. []. GO running on that.
- CAFE GO run in parrelel.
- Nextlfow schema added (replacing nf-validation)

### Changed
- New readme output file format added. 

## [v2.0.4] - 2026-03-28 (2)

### Fixed
- minor bug fix in config. Updated readme.

## [v2.0.3] - 2026-03-28

### Added
- Software versions now collected via `topic: versions` and written to `pipeline_info/software_versions.yml`
- `EGGNOG_TO_GO` now outputs isoform-level GO annotations (`*.isoform_go.txt`) in addition to gene-level
- `versions` topic emit added to `RENAME_FASTA`, `EGGNOG_TO_GO`, `EGGNOG_TO_OG_GO`, `RESCALE_TREE`, `SUMMARIZE_CHROMO_GO`, `CAFE`, and `CAFE_PLOT` modules

### Fixed
- QUAST now receives genome assembly FASTA instead of protein FASTA
- Removed orphaned heredoc version blocks from `CAFE` script

### Changed
- `EGGNOG_TO_GO` now receives pre-filtered GFF (`AGAT_CONVERTSPGXF2GXF.out.output_gff`) to capture all isoforms
- Replaced deprecated `CUSTOM_DUMPSOFTWAREVERSIONS` module with native `topic: versions` approach
- Legacy unused modules moved to `modules/local/legacy_modules/`

### Removed
- Legacy modules and removed old input config input types

## [v2.0.2] - 2026-03-26

### Added
- Added a summarisation R script for chromoGO to capture the general patterns of GO enrichment across chromosomes.

### Fixed
- Added parrelelisation to CHROMO_GO. 

## [v2.0.1] - 2026-03-22

### New features
- Added function to adapt names of samples from AUGUSTUS type GFFs.
- Fixed missing GO plot errors.
- Schema and config updates to new modules and usage
- New gffmix input test to make sure it works with a variety of gff types.
- Added missing containers for some modules, and with ps to work with nextflow

## [v2.0.0] - 2026-03-20

### New features
- Added EggNOG-mapper integration for GO annotation (`--run_eggnog`)
- Added automatic EggNOG database download (`EGGNOG_DOWNLOAD`)

### Enhancements
- Updated OrthoFinder module to v3.1.3 with corrected output paths
- Added `AGAT_CONVERTSPGXF2GXF` step for robust GFF standardisation
- Improved CAFE R scripts to support both OrthoFinder v2 and v3 output formats
- Simplified workflow — removed unused modules (GET_DATA, GO_ASSIGN, GO_EXPANSION)
- CHROMO_GO now uses ORTHOFINDER_CAFE results instead of waiting for ORTHOFINDER_GO

### Bug fixes
- Fixed meta map handling for NCBIGENOMEDOWNLOAD outputs
- Fixed stop codon (`.`) removal from protein fastas before OrthoFinder/EggNOG
- Fixed EggNOG database download. Using wget instead, which is stable. 
- Fixed OrthoFinder v3 output path changes (`Orthogroups/` vs `WorkingDirectory/Orthogroups/`)
- Fixed cafe_go.pl and cafe_prep.R for OrthoFinder v3 column name changes

### Parameters added
- `--run_eggnog` — enable GO annotation via EggNOG-mapper
- `--eggnog_data_dir` — path to pre-downloaded EggNOG database

## [1.1.0] - 2024-03-15
### Added
- New CAFE plotting module
- Support for gzipped input files
- Tested on nextflow version 25 

### Changed
- Readme to explain v25 is best to use

### Fixed
- Bug where GFF files with spaces in names failed

### Removed
- 

## [1.0.1] - 2023-11-02
### Fixed
- Memory issue in ORTHOFINDER step


For a pipeline like excon, the simplest approach is just to update CHANGELOG.md manually whenever you make a new release.
