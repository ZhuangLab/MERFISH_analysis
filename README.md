# MERFISH_analysis
This project contains a series of matlab functions for the analysis of MERFISH as described 
in the publication [*Spatially resolved, highly multiplexed RNA profiling in single cells*](http://www.ncbi.nlm.nih.gov/pubmed/25858977).

Example data and analysis results can be downloaded from the [MERFISH page](http://zhuang.harvard.edu/merfish/). 

## License Information
This software is provided for use to the academic community free of charge. 

Details concerning the licensing of this software for commercial purposes can be found in the license.pdf file. 

## Dependencies
This software requires the following software packages.

1. [matlab-storm](https://github.com/ZhuangLab/matlab-storm): A set of matlab functions for the analysis and visualization of storm-data.
2. [storm-analysis](https://github.com/ZhuangLab/storm-analysis): A collection of STORM analysis functions, primarily in Python.
3. [export-fig](https://github.com/altmany/export_fig): A figure export package for Matlab. 
4. Matlab.
5. Python 2.7 (Only for storm-analysis)

## Installation Instructions

1. Clone the matlab-storm and storm-analysis projects and install.
2. Create or modify a matlab startup script to define all necessary paths for matlab-storm. See the instructions [here](https://github.com/ZhuangLab/matlab-storm/blob/master/README.md).
3. Install the necessary dependencies for export-fig. Instructions can be found [here](https://github.com/altmany/export_fig).

## Examples
Several example scripts are provided to illustrate how to perform a series of basic MERFISH-related tasks.  These scripts can be found in *example_scripts*.

1. analysis_script: This script illustrates the analysis of MERFISH data. It requires [example data](http://zhuang.harvard.edu/merfish/). It produces a series of analyzed data files that contain properties of the decoded barcodes.
2. code_construction_script: This script illustrates the construction of different binary encoding schemes.

## Questions
Contact Jeffrey Moffitt (lmoffitt at mcb.harvard.edu) or Alistair Boettiger (boettiger at fas.harvard.edu).


## Updates
* 2016-07-16, replaced TargetGeneSeqs.fasta in Example MERFISH Data at [http://zhuang.harvard.edu/merfish/](http://zhuang.harvard.edu/merfish/). Previous file listed concatenation of probes, not the sequences of specific isoforms.  