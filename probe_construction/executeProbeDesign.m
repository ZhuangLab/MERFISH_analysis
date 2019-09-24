% Generate probeDesign object

% Initialize probe design object for human panel
% pd = probeDesign(libraryName, species, codebookPath);

% Initialize empty probe design object
% pd = probeDesign();
% Match previous log file
% logFilePath = 'D:\Data\MERFISH\Homosapiens\SMT-H-1002_isoSpec_70-100\SMT-H-1002_isoSpec_70-100.log';
% pd.matchLogFile(logFilePath);
% set(pd, 'species', 'Homo sapiens');
% set(pd, 'libraryName', 'SMT-H-1002_isoSpec_70-100_overlap-20+doubleHead', 'probeSpacing', -20, 'doubleHeadedsmELT', true);
% pd.buildLibrary();

% % Default from demo
% pd = probeDesign('L1E4', 'human', 'C:\Users\ScanningLabAnalysis\Documents\MATLAB\MERFISH_analysis\MERFISH_Examples2\codebook.csv');
% set(pd, 'regionGC', [0.43, 0.63], 'regionTm', [66,76], 'isoSpecificity', [0.75, 1], 'specificity', [0.75, 1], 'probeSpacing', 0);
% set(pd, 'FPKMabundanceThreshold', 1e-2, 'numProbesPerGene', 92);
% set(pd, 'rawTranscriptomeFasta', 'C:\Users\ScanningLabAnalysis\Documents\MATLAB\MERFISH_analysis\MERFISH_Examples2\transcripts.fasta', ...
%         'fpkmPath', 'C:\Users\ScanningLabAnalysis\Documents\MATLAB\MERFISH_analysis\MERFISH_Examples2\isoforms.fpkm_tracking', ...
%         'ncRNAPath', 'C:\Users\ScanningLabAnalysis\Documents\MATLAB\MERFISH_analysis\MERFISH_Examples2\Homo_sapiens.GRCh38.ncrna.fa', ...
%         'readoutPath', 'C:\Users\ScanningLabAnalysis\Documents\MATLAB\MERFISH_analysis\MERFISH_Examples2\readouts.fasta');
%     
%  set(pd, 'transcriptomeHeaderType', 'cufflinks', 'useUniformWeights', false, 'versionMatch', true);
%     
 
 
pd = probeDesign('MO4E4v2', 'mouse', 'C:\Users\ScanningLabAnalysis\Documents\MATLAB\MERFISH_analysis\MERFISH_Examples2\codebook.csv');
set(pd, 'regionGC', [0.43, 0.63], 'regionTm', [66,76], 'isoSpecificity', [0.75, 1], 'specificity', [0.75, 1], 'numProbesPerGene', 92)
set(pd, 'probeSpacing', -20);
set(pd, 'FPKMabundanceThreshold', 0)
 
set(pd, 'codebookPath', 'D:\Data\MERFISH\Musmusculus\MO4E4\M1_codebook_AIBSFormat.csv')


pd.buildLibrary();
