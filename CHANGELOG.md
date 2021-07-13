# Changelog

Here, the changes to `ANANSE` will be summarized.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/).

## [Unreleased]

### Added

- `ananse view` command to view the `binding.h5` file that is now produced by `ananse binding`.
- Support for region /table file as input to `ananse binding`.
- In-built support for mouse.
- Warning with information if another species than human or mouse is used.
- Warning if annotation files don't match for `ananse influence`.
- Improved logging messages.
- Better checking of input files.


### Changed

- `ananse binding` produces a HDF5 file (`binding.h5`) which is much smaller on disk.
- Better memory performance of `ananse network`.
- Removed threshold for differential network in `ananse influence`.

### Removed

### Fixed

