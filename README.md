AIBS fork of code with modifications to support robust design of MERFISH probes. 
Forked by Rusty Nicovich 4/18.  Continued edits to support sequential probes, mouse transcriptome sense, logging, and object-oriented execution.

Key edits and notes:

- Removal of dependence on matlab-storm files.
- Removal of code-initiated edits to MATLAB path.
- Bug fixes to allow execution of code as published on original repo. 
- Collected all apparent parameters into named variables at top of library_design_example.m script.
- Completed implementation for using Ensembl and RefSeq (Entrez) transcriptome files + gene names.
- Sequentially-probed genes (aka smELT) in same or separate panel as barcoded (aka MERFISH) genes supported.
- Bug fix to allow large binary objects to be saved and loaded to/from disk. 
- Support for un-sliced mouse or human transcriptome to pass through analysis.
- Tested with publicly-available mouse transcriptome and non-coding RNA files.
- When not thresholding on any gene abundance (threshold >= 0) then *bulk sequencing is not needed* (proxy file required, but can be made from transcriptome directly). 
- Renamed library_design_example.m to MERFISHProbeDesign.m and moved to .\probe_construction folder.
- Refactored MERFISHProbeDesign to accept probeDesign object (from .\probe_construction\probeDesign.m) in addition to running as script with internal variable values.
- Add output of all probes ([libraryName]_AllOligos.fasta) generated in calculations.
- Add output of message log ([libraryName].log) with error logging.
- Begin implementation of reading on-disk files when specified (preliminary).
- Additional options for filtering generated targetRegions.  Options for filtering by parameter (mimicking default), expanding isoSpecificity to generate desired number of probes, or selecting probes from common sequences across isoforms, when available. 

-------------------------------------------------------------------------------------------

Required inputs for given run:

Transcriptome file (ex for mouse):

ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/cdna/


Non-coding RNA file (ex for mouse):

ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/ncrna/


Codebook file w/ gene name, isoform, and binary codes as to MERFISH +- smELT format

./MERFISH_anlaysis/MERFISH_Examples2/codebookMusmusculusHypothalamus_v01.csv


Readouts.fasta file

./MERFISH_anlaysis/MERFISH_Examples2/readoutsHypothalamus.fasta


Proxy bulk sequencing file (created with ./mFISHcodebookPermutation/makeProxyFPKMfileFromRefSeqfile.py):

./MERFISH_analysis/MERFISH_Examples2/Mus_musculus_proxy.fpkm_tracking


Requires all sub-folders of this repo be on MATLAB path to execute.


-------------------------------------------------------------------------------------------

Most straightforward execution is to collect variables into probeDesign object and then call probeDesign.buildLibrary() method. 
This can be done in script with parameters edited using set-get methods. Convenience methods are included to start from reference log file or from default inputs for given species.
Example:
```
codebookPath = './MERFISH_anlaysis/MERFISH_Examples2/codebookMusmusculusHypothalamus_v01.csv';
pd = probeDesign('hypothalamusLibrary', 'mouse', codebookPath);
set(pd, 'regionGC', [0.43, 0.63], 'regionTm', [66,76], 'isoSpecificity', [0.75, 1], 'specificity', [0.75, 1]);
set(pd, 'FPKMabundanceThreshold', 0, 'numProbesPerGene', 92);
pd.buildLibrary();
```	
Or:
```
% Initialize empty probe design object
pd = probeDesign();
% Match previous log file
logFilePath = 'D:\Data\MERFISH\Homosapiens\SMT-H-1002_isoSpec_70-100\SMT-H-1002_isoSpec_70-100.log';
pd.matchLogFile(logFilePath);
set(pd, 'species', 'Homo sapiens');
set(pd, 'libraryName', 'SMT-H-1002', 'probeSpacing', -20, 'doubleHeadedsmELT', true);
pd.buildLibrary();
```

Details on variables in probeDesign object:

Required variables:
```
                          MERFISHAnalysisPath: String. Path to load directory for rawTranscriptomeFasta, fpkmPath, ncRNAPath, transcriptomeHeaderType.  
												This folder + sub-folders should be on MATLAB path.
                                     basePath: String. Path to save directory for outputs. 
                                  libraryName: String. Name to call library. Prefix for output files.
                                      species: String. Which species this is.  For saving and convenince with default paths of rawTranscriptomeFasta, fpkmPath, ncRNAPath,
												transcriptomeHeaderType and transcriptomeIDType.
                        rawTranscriptomeFasta: String. Path to transcriptome.fasta file.  Format must match transcriptomeHeaderType and transcriptomeIDType.  
												Path must be to valid fasta file. 
                                     fpkmPath: String. Path to isoforms.fpkm_tracking file.  Format must match transcriptomeHeaderType and transcriptomeIDType.  
												Transcript IDs and gene names must match those in rawTranscriptomeFasta.  Path must be to valid file.
                                    ncRNAPath: String. Path to species.GRCh38.ncrna.fa file.  Path must be to valid file. 
                                  readoutPath: String. Path to readouts.fasta file containing readout sequences and names.  Names must match those in Codebook.  
												Path must be to valid fasta file. 
                                 codebookPath: String. Path to input codebook file.  Gene names must match that in fpkmPath.  Gene IDs must match those in rawTranscriptomeFasta. 
												If transcript version numbers do not match rawTranscriptomeFasta, set versionMatch below to false.
												Readout names must match those in readoutPath.  Path must be to valide codebook file. 
												Format of codebook controlled.  See examples and schema for reference.
                      transcriptomeHeaderType: String. 'cufflinks', 'ensembl', or 'refseq'
                          transcriptomeIDType: String. 'NCBI' ['cufflinks' or 'refseq' header type] or 'ENSEMBL' ['ensembl'] 
           penaltyTableExactHomologyToExclude: Double. Used in primer design + rRNA/tRNA penalty table.
                       FPKMabundanceThreshold: Double. Min expression threshold in bulk sequencing data for consideration in isoSpecificity calculations. 
												Set to 0 to include all possible transcripts.
                            useUniformWeights: Boolean. Set to true to include fpkmPath file contents in isoSpecificity calculation. 
												Set to false (and FPKMabundanceThreshold to 0) to ignore contents of fpkmPath file. 
    isoSpecificityTable_lengthOfExactHomology: Double. Length of exact homology used to calculate isoSpecificity penalties.
                                 regionLength: Double.  Length of probing region. Used in trDesigner.DesignTargetRegions.
									 regionGC: [Double Double].  Range for GC content of accepted probe regions. Used in trDesigner.DesignTargetRegions.
                                     regionTm: [Double Double].  Range for Tm value of accepted probe regions. Used in trDesigner.DesignTargetRegions.
                               isoSpecificity: [Double Double].  Range of uniqueness of regions to probe in selected transcript(s). Measured between isoforms of same gene.  
												Can include abundance information from fpkmPath file if useUniformWeights set to true.  Used in trDesigner.DesignTargetRegions.
                                  specificity: [Double Double].  Range of uniqueness for selected transcript(s) against all other genes.  Used in trDesigner.DesignTargetRegions.
                  monovalentSaltConcentration: Double. Concentration (in mol/L) of monovalent salt for Tm calc.  Used in trDesigner.DesignTargetRegions and PrimerDesigner.
                           probeConcentration: Double. Concentration (in mol/L) of probe for Tm calc.  Used in trDesigner.DesignTargetRegions and PrimerDesigner.
                                 probeSpacing: Double. Gap between adjacent probes on target transcripts, in nT.  Used in trDesigner.DesignTargetRegions.
                             numProbesPerGene: Double.  Number of probes to retain per target gene. Total across all transcripts given for that gene. 
                           nPrimersToGenerate: Double.  Number of primers to generate before culling to acceptable selected pool. Used in PrimerDesigner.
                                 primerLength: Double.  Length of primers, in nT. Used in PrimerDesigner.
                                 cutPrimersTm: [Double Double].  Range of Tm for generated primer regions. Used in PrimerDesigner.cutPrimers.
                                 cutPrimersGC: [Double Double].  Range of GC content for accepted primer regions. Used in PrimerDesigner.cutPrimers.
                    cutPrimersMaxHomologySelf: Double. Max self-homology length accepted in primers. Used in PrimerDesigner.RemoveSelfCompPrimers.
                   cutPrimersMaxHomologyCross: Double. Max cross-homology length accepted in primers. Used in PrimerDesigner.RemoteHomologousPrimers.
                                 versionMatch: Boolean.  If true, codebook and transcriptome must have exact match between transcript ID and version number.  
												If false, only transcript ID need match.  Ex : 'ENST00000374472.4' must match if true, but only 'ENST00000374472' if false.
												Used in oligomers build step.
                            doubleHeadedsmELT: Boolean.  If true, append same readout to both ends of transcript if number of probes for this gene < nProbesPerGene AND
												Codebook has 1 positive bit for gene.  Used in oligomers build step. 
						keepAllPossibleProbes: Boolean.  If true, output allOligos.fasta file containing all possible probe sequences.
									debugMode: Boolean.  If true, assign useful variables to base workspace from MERFISHProbeDesign function.
					 readoutPermuteBySequence: Boolean.  If true, use modulo of sequence to assign readouts to probe.  If false, use randperm().
							   spaceOutProbes: Boolean.  If true, and probes > numProbesPerGene, select probes which are most spaced out.  If false, select randomly. 
							  specifyReadouts: Boolean.  If true, use readout names specified in header of codebook file.  If false, use order from readouts.fasta file.
						geneIsoformListSource: String. Source of gene + isoform pairs to filter in targetRegions filters. 'default', 'allGenes', or 'codebook'.
                               tRFilterMethod: String. Method for filtering targetRegions. 'default', 'parameter','relaxIsospecificity', or 'commonRegions'.
                                tRFilterField: String. If tRFilterMethod = 'parameter', which parameter field to filter over. 'regionLength','GC', 'Tm', 
												'specificity', 'isoSpecificity', 'none'
                           tRFilterParameters: [Double Double]. If tRFilterMethod = 'parameter', range for filtering targetRegions. 
								rRNAtRNAPath: String.  File path to previously-generated rRNAtRNAPath file.  Will be loaded if inputs in object match those of
													saved.             	 
                            transcriptomePath: String.  File path to previously-generated transcriptome object file.  Will be loaded if inputs in object match those of saved.
                         specificityTablePath: String.  File path to previously-generated specificity table object file.  Will be loaded if inputs in object match those of saved.
                      isoSpecificityTablePath: String.  File path to previously-generated isoSpecificity table object file.  Will be loaded if inputs in object match those of saved.
                               trDesignerPath: String.  File path to previously-generated trDesigner object folder.  Will be loaded if inputs in object match those of saved.
							   
```
Defaults in GitHub repo code:
```
                          MERFISHAnalysisPath: ['C:\Users\Jeff.Morgan0\Dropbox\ZhuangLab\MERFISH_Public\MERFISH_analysis\']; Specified in .\startup\merfish_startup.m line 56.
                                     basePath: [MERFISHAnalysisPath '\MERFISH_Examples2\']; 			Specified in .\example_scripts\library_design.m line 15.
                                  libraryName: ['L1E1']; 												Specified in .\example_scripts\library_design.m line 294.
                        rawTranscriptomeFasta: [basePath 'transcripts.fasta']; 							Specified in .\example_scripts\library_design.m line 20.
                                     fpkmPath: [basePath 'isoforms.fpkm_tracking']; 					Specified in .\example_scripts\library_design.m line 21.
                                    ncRNAPath: [basePath 'Homo_sapiens.GRCh38.ncrna.fa']; 				Specified in .\example_scripts\library_design.m line 22. 
                                  readoutPath: [basePath 'readouts.fasta']; 							Specified in .\example_scripts\library_design.m line 23. 
                                 codebookPath: [basePath 'codebook.csv']; 								Specified in .\example_scripts\library_design.m line 24. 
                      transcriptomeHeaderType: 'cufflinks'; 											Specified in .\probe_construction\Transcriptome.m line 71. 
                          transcriptomeIDType: 'NCBI'; 													Specified in .\probe_construction\Transcriptome.m line 72.
           penaltyTableExactHomologyToExclude: 15; 														Specified in .\example_scripts\library_design.m line 195.
                       FPKMabundanceThreshold: 1e-2; 													Specified in .\example_scripts\library_design.m line 211.
                            useUniformWeights: False; 							Specified in .\probe_construction\OTTable.m line 87, .\example_scripts\library_design.m line 136.
    isoSpecificityTable_lengthOfExactHomology: 17; 														Specified in .\example_scripts\library_design.m line 134.
                                 regionLength: 30; 														Specified in .\example_scripts\library_design.m line 235.
									 regionGC: [0.43, 0.63]; 											Specified in .\example_scripts\library_design.m line 236.
                                     regionTm: [66, 76]; 												Specified in .\example_scripts\library_design.m line 237.
                               isoSpecificity: [0.75, 1]; 												Specified in .\example_scripts\library_design.m line 238.
                                  specificity: [0.75, 1]; 												Specified in .\example_scripts\library_design.m line 239.
                  monovalentSaltConcentration: 0.3; 													Specified in .\probe_construction\TRDesigner.m line 800.
                           probeConcentration: 5e-9; 													Specified in .\probe_construction\TRDesigner.m line 801.
                                 probeSpacing: 0; 														Specified in .\probe_construction\TRDesigner.m line 804 as 'threePrimeSpace'.
                             numProbesPerGene: 92; 														Specified in .\example_scripts\library_design.m line 293.
                           nPrimersToGenerate: 1000; 													Specified in .\example_scripts\library_design.m line 437.
                                 primerLength: 20 														Specified in .\example_scripts\library_design.m line 438.
                                 cutPrimersTm: [70, 72] 												Specified in .\example_scripts\library_design.m line 444.
                                 cutPrimersGC: [0.5, 0.65] 												Specified in .\example_scripts\library_design.m line 445.
                    cutPrimersMaxHomologySelf: 6 														Specified in .\example_scripts\library_design.m line 448.
                   cutPrimersMaxHomologyCross: 8 														Specified in .\example_scripts\library_design.m line 449.
```
-------------------------------------------------------------------------------------------

Original README continues below:

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
5. Python 2.7 (Only for storm-analysis)

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
