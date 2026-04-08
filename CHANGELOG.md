# Changelog

## [v2.3.0] - 2026-04-07

### Added
- GO enrichment results (`CAFE_GO_RUN`, `CAFE_GO_PREP_LARGE`, `CAFE_GO_RUN_LARGE`, `CHROMO_GO`, `SUMMARIZE_CHROMO_GO`) are now published to subfolders named after the active GO settings (e.g. `cafe_go/weight01_t_cutoff0.05_typenone/`). Re-running with different `--go_algo`, `--go_cutoff`, or `--go_type` values will write to a separate subfolder, preserving results from all parameter combinations.
- `SUMMARIZE_CHROMO_GO` now handles the case where no significant GO terms are found (no `*_res.tab` files produced by `CHROMO_GO`). The process exits cleanly with a message instead of erroring, and its outputs are marked optional so the pipeline continues.
- `CAFE_MODEL_COMPARE` now publishes the winning model's CAFE5 results directory to `results/cafe/best/`, mirroring the layout of `results/cafe/base/`. The complex channel derivation of `ch_best_results` in the workflow is replaced by a direct `CAFE_MODEL_COMPARE.out.best_results` emit.

### Added
- New `OG_ANNOTATION_SUMMARY` module: when `--run_eggnog` is set, produces `OG_annotation_summary.tsv` with one row per orthogroup containing the representative gene ID, description, preferred name, COG category, KEGG KO, and PFAM domains from EggNOG-mapper annotations. Output is written to `results/eggnogmapper/OG_annotation_summary.tsv`.
- New `--eggnog_rep_species` parameter to force the representative species used for orthogroup annotation. When unset, the species with the most annotated genes is chosen automatically.
- `EGGNOGMAPPER` output (`.emapper.annotations` files) is now published to `results/eggnogmapper/`.
- `CAFE_RUN_LARGE` now emits an optional `converged.txt` sentinel. `CAFE_PLOT_LARGE` and `CAFE_GO_PREP_LARGE` only run when convergence was achieved, preventing downstream failures when high-differential families cannot be modelled.
- `CAFE_RUN_LARGE` lambda retry sweep extended: now tries estimated λ, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001, 0.000001, 0.0000001 (previously only descended; now also tries larger values which are more physically appropriate for families with extreme size differentials).
- `EGGNOG_TO_OG_GO` now parallelises GO file reading using `ProcessPoolExecutor` with `task.cpus` workers, reducing runtime on large species sets.

### Changed
- `EGGNOG_TO_OG_GO` label changed from `process_single` to `process_high` (8 CPUs) to support parallel GO file reading.

## [v2.2.0] - 2026-04-06

### Added
- New OrthoFinder algorithm parameters: `--orthofinder_method` (`-M`), `--orthofinder_search` (`-S`), `--orthofinder_msa_prog` (`-A`), and `--orthofinder_tree` (`-T`). These map directly to OrthoFinder command-line flags and are all optional — OrthoFinder defaults are used when unset.
- New `--orthofinder_v2` flag (default `false`) to run OrthoFinder v2.5.5 instead of v3.1.3. v2 uses Hierarchical Orthogroups (`N0.tsv`) which are more appropriate for CAFE5 as they represent gene families traceable to the common ancestor. v3 uses flat orthogroups (`Orthogroups.tsv`) which can have inflated copy-number variance. For large datasets (>30 species), v2 is recommended.
- `ORTHOFINDER_V2_CAFE` results are now published to `results/orthofinder_cafe/` (was only published for v3).
- `CAFE_PREP` now emits `pruned_tree` (the rescaled, species-name-stripped tree) for use by all downstream CAFE runs.
- `CAFE_RUN_LARGE` now retries with progressively smaller lambda values (estimated → 1e-4 → 1e-5 → 1e-6 → 1e-7) when the initial fixed-lambda run fails to converge, as recommended in hahnlab/CAFE5#132.
- CAFE GO enrichment plots now display full GO term text labels and GO IDs.
- New output documentation page `docs/outputs.md` with example figures and detailed descriptions of all output files.

### Changed
- Reverted tree scaling back to the original `RESCALE_TREE` approach (`rescale_tree.py` multiplies branch lengths by `--tree_scale_factor`). The `MAKE_ULTRAMETRIC` module introduced in v2.1.4 is removed.
- `--tree_scale_factor` default changed from `1` back to `1000`.
- `CAFE_PREP` base run and error model estimation now use `pruned_tree` (the rescaled non-ultrametric tree) directly, matching the approach that was validated on large datasets. `SpeciesTree_rooted_ultra.txt` (produced by `chronoMPL()`) is retained for reference only.
- `CAFE_RUN_K` and `CAFE_RUN_BEST` now use `pruned_tree` instead of `SpeciesTree_rooted_ultra.txt`, avoiding convergence failures caused by the ultrametric tree's maximum-possible-lambda constraint.
- `CAFE_RUN_LARGE` failure is now non-fatal — the pipeline continues even if high-differential families cannot be modelled.
- `CAFE_PLOT` (and `CAFE_PLOT_LARGE`) now skip gracefully when CAFE5 did not produce an `*_asr.tre` file, instead of crashing the pipeline.
- `ORTHOFINDER_V2` module now emits `N0.tsv` as `orthologues` (previously emitted `Orthogroups.tsv`). This ensures CAFE_PREP receives hierarchical orthogroups, which have lower copy-number variance and are required for correct CAFE5 analysis.
- `CHROMO_GO` now symlinks `N0.tsv` to `Orthogroups.tsv` when the v2 path is used, so the downstream perl script works regardless of OrthoFinder version.

### Fixed
- Fixed species name mismatch between tree and gene counts in `CAFE_RUN_K`: the `SpeciesTree_rescaled.nwk` retains `.clean` suffixes but `hog_gene_counts.tsv` uses plain names. All CAFE runs now use `pruned_tree` which has suffixes stripped by `sed` in `CAFE_PREP`.
- Fixed `chronos()` convergence failure on large datasets: passing the ×1000 pre-scaled tree to `chronos()` caused a degenerate starting point (~-74 billion log-likelihood). `CAFE_PREP` now receives the original unscaled tree and scaling is applied after any ultrametric correction.

## [v2.1.4] - 2026-04-03

### Fixed
- Replaced `RESCALE_TREE` + `chronos()` with `MAKE_ULTRAMETRIC`, fixing incorrect ultrametric trees for OrthoFinder v3 output. The previous approach multiplied branch lengths by 1000 and then re-estimated them from scratch with `chronos()`, producing near-zero values (0.0001…) that broke CAFE5. The new module uses the OrthoFinder ultrametricization algorithm to adjust branch lengths while preserving relative topology.

### Changed
- `--tree_scale_factor` now sets the root-to-tip distance in the ultrametric output tree (default `1`) rather than a raw branch-length multiplier. Remove any `tree_scale_factor = 1000` from your config — the default of `1` is correct for CAFE5.
- `cafe_prep.R` now asserts the input tree is already ultrametric rather than attempting to correct it with `chronos()`.

## [v2.1.3] - 2026-04-03

### Removed
- `AGAT_CONVERTSPGXF2GXF` module removed — `AGAT_SPKEEPLONGESTISOFORM` already sanitises the GFF, making the conversion step redundant

### Changed
- `AGAT_SPKEEPLONGESTISOFORM` now publishes its sanitised GFF output to `agat/`

## [v2.1.2] - 2026-04-02

### Added
- New `--go_algo` parameter to select the topGO algorithm and statistic combination: `classic_fisher` (default), `weight01_t`, `elim_ks`, or `weight_ks`. Applied to both CAFE GO and chromosome GO analyses.

## [v2.1.1] - 2026-04-01

### Fixed
- Allow Augustus transcript IDs to output to GO hash 
- Force eggnogmapper to process_high only, High memory causes stalling issues

## [v2.1.0] - 2026-03-30

### Added
- CAFE now runs k=1 through k=`cafe_max_k` (default 6) rate categories in parallel (`CAFE_RUN_K`), then selects the best k by AIC (`CAFE_SELECT_K`).
- Best k is re-run with the Poisson birth-death option (`-p`) via `CAFE_RUN_BEST`, and the winner between uniform and Poisson is selected by likelihood (`CAFE_MODEL_COMPARE`).
- CAFE GO enrichment now runs in parallel, one job per species/node and direction (`CAFE_GO_RUN`).
- Nextflow schema added (replacing nf-validation).

### Changed
- Updated output documentation structure in README.

### Fixed
- Fixed flag names for EggNOG parameters.

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
- `EGGNOG_TO_GO` now receives the raw input GFF to capture all isoforms
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
- `AGAT_SPKEEPLONGESTISOFORM` handles GFF sanitisation (making a separate convert step redundant)
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
