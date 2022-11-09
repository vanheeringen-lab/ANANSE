# Changelog

Here, the changes to `ANANSE` will be summarized.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/).

## [Unreleased]

### Added
- `anane view` can now return factory activity with `--activity`
- `anane view` also accepts `-n` for all sub-option (list TFs, list regions and activity)

### Changed
- `ananse network` uses slightly less memory, and should now keep going unless all memory is used.

### Fixed
- `ananse influence` now throws a descriptive error if source and target networks are identical
- `ananse influence` now uses scipy >=1.9, which fixes an error in `mannwhitneyu(method="auto")`
- uploaded correct non-human CAGE-seq models

## [0.4.0] - 2022-06-02

### Added
- CAGE-seq support!
  - for more details, check out https://github.com/vanheeringen-lab/ANANSE-CAGE
- `ananse.binding` can now optionally accept a raw counts file instead of BAM files.
  - `ananse.binding` now has a (case-insensitive) `--columns` argument to filter the counts table by.
- `ananse.view` now optionally accepts `regions` and `tfs`, to filter the output.
- `ananse.view` can also output the header of (up to) `n` TFs & regions.
- `ananse.view` can also output a list of all regions or tfs.
- `ananse network` will check the gene name overlap between the expression files, the motif2factor.txt and the tfs (#120).
  - if overlap is low, and a genomepy genome/annotation is given, the symbols are converted to gene names, similar to `gimme motif2factors`
  - if overlap is none, an informative error is raised
- `ananse network` can now accept one or more column names to use from the expression file (e.g. column names from a counts table)
- `ananse network` now optionally accepts `regions` and `tfs`, to filter the binding.h5 content.
- `ananse network` can (at least internally) output an expression network, binding network or expression-binding network.
- `ananse influence` will also check the gene overlap (between the merged networks and the DEgenes) (#120).
  - if overlap is low, and a genomepy genome/annotation is given, the symbols are converted to gene names, similar to `gimme motif2factors`
  - if overlap is none, an informative error is raised
- `ananse plot` now optionally accepts a file `type` for the output plots (e.g. pdf, png)
- all ananse CLI functions now have a default output file/dir
- **ADDITIONAL UNIT TESTS!**
- **TEST DATA!** See `tests/data/GRCz11_chr9` 
  - including a README.md on the creation and bash command for every function.

### Changed

- Now uses em-dash as gene separator (—), instead of underscore (_). _Should_ solve issues with gene names having underscores, as em-dash is a pretty obscure char.
- `ananse binding` will ignore BLACKLIST_TFS in motif2factors.txt. Currently, this is only the artifact `"NO ORTHOLOGS FOUND"`.
- refactored jaccard stuff
  - reduced jaccard graph size (no more duplicates/self edged)
  - reduced logger messages (max 1 message per motif now)
- `ananse binding` documentation changed to clearly reflect the complexity of non-`hg38`
  - with expandable sections to reduce clutter for the `hg38` gang
- `ananse binding` will now only scan the motifs in the analysis
- `ananse binding` will now only scan the regions overlapping the pfmscorefile and regions (if both are given)
- changed the hardcoded `.txt` file extensions in `ananse plot` to `.tsv` 
- `ananse network` now loads the gene bed into pyranges once (instead of once per chromosome).
- `ananse influence` now uses the top 500.000 TF-gene interactions. By default, filtering for highest positive interactions.
- cleaned up `ananse --help` messages

### Removed

- `ananse binding` now requires a genome fasta (to create `_factor_activity`)
- legacy code (`enhancer_binding` is now `bed`, containing only the relevant code)

### Fixed
- jaccard index with custom pfmfiles
- `ananse view --help`
- issue in `ananse binding` when specified regions have no bam coverage (#153)
- issue in `ananse binding` when probability cant be found for any TF (#147)
- issue in `ananse network` when executing run_network() from python
- `ananse network` merges duplicate genes and transcription factors (#142)
- issue in `ananse influence` when using only 1 network
- `ananse influence` made ~~slightly~~ much faster.
- `ananse influence` hopefully uses less memory now.
- `ananse influence` now skips pvalues for TFs without any targets or non-targets.
- `ananse plot` will no longer warn you incorrectly about your "weight" 
- `ananse plot` error "OSError: Format: "dot" not recognized."


## [0.3.0] - 2021-07-14

### Added

- `ananse view` command to view the `binding.h5` file that is now produced by `ananse binding`.
- Support for region /table file as input to `ananse binding`.
- In-built support for mouse.
- Warning with information if another species than human or mouse is used.
- Warning if annotation files don't match for `ananse influence`.
- Improved logging messages.
- Better checking of input files.
- The -full--output option for both ananse network and ananse influence for additional info in the output files and more visualization options
- An optional jaccard similarity cutoff for `ananse binding` to set a limit on motif overlap for TF binding model selection
- `ananse plot` command to plot the dotplot and a new top TF interaction GRN file.


### Changed

- `ananse binding` produces a HDF5 file (`binding.h5`) which is much smaller on disk.
- Better memory performance of `ananse network`.
- Removed threshold for differential network in `ananse influence`.
- ananse influence now takes the set of  interactions of the top 500.000 of each network instead of the head of each.
- ananse influence now doesnt have plot functionality, moved this to ananse plot instead.


### Fixed

- Gene names don't get capitalized in `ananse influence` (#87).
