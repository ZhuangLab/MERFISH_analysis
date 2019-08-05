# MERFISH_analysis
This project contains a series of matlab functions for the analysis of MERFISH as described 
in the publication [*Spatially resolved, highly multiplexed RNA profiling in single cells*](http://www.ncbi.nlm.nih.gov/pubmed/25858977) 
as well as the publication [*High-throughput single-cell gene-expression profiling with multiplexed error-robust fluorescence in situ hybridization*](https://www.ncbi.nlm.nih.gov/pubmed/27625426).

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
5. Python 3.6 (Only for storm-analysis)

## Installation Instructions

1. Clone the matlab-storm and storm-analysis projects and install.
2. Create or modify a matlab startup script to define all necessary paths. Instructions for matlab-storm can be found [here](https://github.com/ZhuangLab/matlab-storm/blob/master/README.md). 
Additional paths are required for MERFISH_analysis. An example script is also provided: startup\merfish_startup
3. Install the necessary dependencies for export-fig. Instructions can be found [here](https://github.com/altmany/export_fig).

## Hardware Requirements
32 GB to 64 GB of RAM are recommended for the construction of target regions using the new pipeline (illustrated in the library_design_example script). 

## Examples
Several example scripts are provided to illustrate how to perform a series of basic MERFISH-related tasks.  These scripts can be found in *example_scripts*.

1. analysis_script: This script illustrates the analysis of MERFISH data. It requires [example data](http://zhuang.harvard.edu/merfish/MERFISHData/MERFISH_Examples.zip). It produces a series of analyzed data files that contain properties of the decoded barcodes.
2. code_construction_script: This script illustrates the construction of different binary encoding schemes.
3. library_design_example: This script illustrates the construction of a MERFISH probe library. It requires [example data](http://zhuang.harvard.edu/merfish/MERFISHData/MERFISH_Examples2.zip).

## Updates
1. (December 2016) A more computationally efficient pipeline for the design and construction of MERFISH probes is now included. 
The functions and example scripts associated with the older pipeline are still available but have been moved to folders named 'deprecated'. 
2. (September 2018) The MERFISH analysis pipeline has been updated to allow massively parallel computing on clusters running the SLURM scheduling system. Improved cell segmentation functionality has been added. 
The functions and example scripts associated with older pipelines can be found in folders named 'deprecated'.

## Questions
Contact Xiaowei Zhuang (zhuang at chemistry.harvard.edu) or Jeffrey Moffitt (lmoffitt at mcb.harvard.edu).
