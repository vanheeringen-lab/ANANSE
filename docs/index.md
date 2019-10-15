## XXX Network for key transcription factor prediction during cell fate switch analysis

### What is XXX Network?
GimmeMotifs is an analysis framework for transcription factor motif analysis written in Python. It contains command-line scripts to predict de novo motifs, scan for known motifs, identify differential motifs, calculate motif enrichment statistics, plot sequence logos and more. In addition, all this functionality is available from a Python API.

GimmeMotifs is free and open source research software. If you find it useful please cite our paper:

* van Heeringen SJ and Veenstra GJC, GimmeMotifs: a de novo motif prediction pipeline for ChIP-sequencing experiments, Bioinformatics. 2011 Jan 15;27(2):270-1. doi: 10.1093/bioinformatics/btq636.


### Getting started
* The easiest way to install GimmeMotifs is using bioconda on Linux or Mac. From version 0.13.0 only Python 3 (>= 3.4) is supported.
* Have a look at these simple examples to get a taste of what is possible.
* Check out the more detailed tutorials.
* Full command-line reference can be found here.
* There’s also an API documentation.


### Get help
* First, check the FAQ for common issues.
* The preferred way to get support is through the Github issues page.
* Finally, you can reach me by mail or via twitter.


### Full contents
- Installation
 + The easiest way to install
 + Alternative installation
 + Configuration
 + Overview
 + Motif databases
 + Simple examples
 + Install a genome
 + Predict de novo motifs
 + Compare motifs between data sets
 + Create sequence logos
- Tutorials
 + Find de novo motifs
 + Scan for known motifs
 + Find differential motifs
 + Compare two sets with de novo motifs
 + Motif enrichment statistics
* Command-line reference
 + List of tools
 + Input formats
 + Command: gimme motifs
 + Command: gimme maelstrom
 + Command: gimme scan
 + Command: gimme roc
* API documentation
* Examples
 + Working with motifs
 + Motif scanning
 + Finding de novo motifs
 + Motif statistics
 + Maelstrom
* Auto-generated
 + The Motif class
 + Prediction of de novo motifs
 + Motif scanning
* FAQ
 + Sorry, motif prediction tool [X] is not supported
 + I get motifs that have differential scores in gimme maelstrom, however, the number is not different across clusters
 + I have upgraded GimmeMotifs and now it doesn’t find my genome
 + I cannot run gimme index anymore
 + I get ‘RuntimeError: Invalid DISPLAY variable’
 + I get a KeyError when running gimme maelstrom
* Acknowledgments