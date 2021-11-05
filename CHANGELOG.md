# Changelog

Here, the changes to `ANANSE` will be summarized.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/).

## [Unreleased]

### Added

### Changed

- Now uses em-dash as gene separator (—), instead of underscore (_). _Should_ solve issues with gene names having underscores, as em-dash is a pretty obscure char.
- refactored jaccard stuff
  - reduced jaccard graph size (no more duplicates/self edged)
  - reduced logger messages (max 1 message per motif now)
- `ananse binding` documentation changed to clearly reflect the complexity of non-`hg38`
  - with expandable sections to reduce clutter for the `hg38` gang
- `ananse binding` will now only scan the motifs in the analysis
- cleaned up `ananse binding --help` messages

### Removed

### Fixed
- jaccard index with custom pfmfiles
- `ananse view --help`


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
