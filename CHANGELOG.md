# Changelog

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
