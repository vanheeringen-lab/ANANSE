## News

### ANANSE v0.3.0
> 2020-11-25

###### Added

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

###### Changed

- `ananse binding` produces a HDF5 file (`binding.h5`) which is much smaller on disk.
- Better memory performance of `ananse network`.
- Removed threshold for differential network in `ananse influence`.
- ananse influence now takes the set of  interactions of the top 500.000 of each network instead of the head of each.
- ananse influence now doesnt have plot functionality, moved this to ananse plot instead.

###### Fixed

- Gene names don't get capitalized in `ananse influence` (#87).

### ANANSE v0.1.5
> 2020-11-25

* Switch the default enhancer input to H3K27ac data.
* Add enhancer command to establish enhancer file.
* Fix a lot of bugs, thanks [@Branco Heuts](https://www.researchgate.net/profile/Branco_Heuts).
* Simplify the readme file.
* Update the ReadtheDocs for ANANSE.  
* You can download it from GitHub: [ANANSE Release v0.1.5](https://github.com/vanheeringen-lab/ANANSE/releases/tag/v0.1.5).

### ANANSE v0.1.4
> 2020-07-29

* Fix genomepy version.
* Fix weight function.
* Fix TF filter.
* Improve docs.
* You can download it from GitHub: [ANANSE Release v0.1.4](https://github.com/vanheeringen-lab/ANANSE/releases/tag/v0.1.4).

### ANANSE v0.1.3
> 2020-06-02

* Fix genomepy chromosome size bug thanks Siebren Frölich ([@siebrenf](https://github.com/siebrenf)).
* Fix influence calculation path bug.
* Fix promoter and expression network caculation bug.
* Fix a lot of bugs, thanks Siebren Frölich ([@siebrenf](https://github.com/siebrenf)) and Jos Smits ([@JGASmits](https://github.com/JGASmits)) for valuable input, testing and bug hunting.  
* Add ReadtheDocs for ANANSE.  
* You can download it from GitHub: [ANANSE Release v0.1.3](https://github.com/vanheeringen-lab/ANANSE/releases/tag/v0.1.3).

### ANANSE v0.1.2
> 2020-04-07

* Add standard loggers.  
* Add a progress bar.  
* Add mutil threads.  
* You can download it from GitHub: [ANANSE Release v0.1.2](https://github.com/vanheeringen-lab/ANANSE/releases/tag/v0.1.2).

### ANANSE v0.1.1
> 2020-04-02

* A umber of minor bug fixes.  
* Release first bioconda version.  
* You can download it from GitHub: [ANANSE Release 0.1.1](https://github.com/vanheeringen-lab/ANANSE/releases/tag/0.1.1).

### ANANSE v0.1.0
> 2020-01-22

* First stable version of ANANSE has been released with various feature improvements.   
* You can download it from GitHub: [ANANSE Release 0.1.0](https://github.com/vanheeringen-lab/ANANSE/releases/tag/0.1.0).
