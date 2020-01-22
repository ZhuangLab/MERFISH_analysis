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
 
 
% pd = probeDesign('Human_40gene_smELT_test201911217', 'human', 'C:\Users\ScanningLabAnalysis\Documents\MATLAB\MERFISH_analysis\MERFISH_Examples2\SMT-H-1002c_codebook.csv');
% set(pd, 'regionGC', [0.43, 0.63], 'regionTm', [66,76], 'isoSpecificity', [0, 1], 'specificity', [0.75, 1], 'numProbesPerGene', 48);
% set(pd, 'probeSpacing', -20);
% set(pd, 'ncRNAPath', 'D:\Data\MERFISH\Homosapiens\Homo_sapiens.GRCh38.ncrna.fa')
% set(pd, 'FPKMabundanceThreshold', 0);
% set(pd, 'rawTranscriptomeFasta', 'D:\Data\MERFISH\Homosapiens\Homo_sapiens_GRCh38_latest_rna.fna');
% set(pd, 'fpkmPath', 'D:\Data\MERFISH\Homosapiens\Homo_sapiens_proxyRandomFPKM.fpkm_tracking');
% set(pd, 'versionMatch', true);
% 
% pd.buildLibrary()


%try
 %   pd = probeDesign();
  %  logFilePath = 'D:\Data\MERFISH\Homosapiens\Human_MTG_Sequential_20191217_DHS\Human_MTG_Sequential_20191217_DHS.log';
   % pd.matchLogFile(logFilePath);
    %set(pd, 'libraryName', 'MTG_20191220seq_fromLog', 'species', 'Homo sapiens');
    %pd.buildLibrary();
    %display('Completed sequential');
%catch mError
 %   display('Error on MTG_20191220_fromLog');
%end

%try
 %   pd = probeDesign('Human_MTG_Barcoded_20191220', 'human', 'D:\Data\MERFISH\Homosapiens\Human_MTG_Panel1_barcoded2.csv');
  %  pd.buildLibrary();
   % display('Completed barcoded');
%catch mError
 %   display('Error on MTG_20191220_fromLog2');
%end
% 
%pd = probeDesign('Mouse_VISp_Barcodedalt_20200118', 'mouse', 'D:\Data\MERFISH\Musmusculus\Mus_musculus_VISp152JLC_barcoded_altered20200115.csv');
%set(pd, 'regionGC', [0.43, 0.63], 'regionTm', [66,76], 'isoSpecificity', [0, 1], 'specificity', [0.75, 1], 'numProbesPerGene', 92);
%set(pd, 'probeSpacing', -20);  
%set(pd, 'ncRNAPath', 'D:\Data\MERFISH\Musmusculus\Mus_musculus.GRCm38.ncrna.fa')
%set(pd, 'FPKMabundanceThreshold', 0);
%set(pd, 'rawTranscriptomeFasta','D:\Data\MERFISH\Musmusculus\Mus_musculus.GRCm38.cdna.all.fa');
%set(pd, 'fpkmPath', 'D:\Data\MERFISH\Musmusculus\Mus_musculus_proxyRandomFPKM.fpkm_tracking');



pd = probeDesign();
logFilePath='D:\Data\MERFISH\Musmusculus\Mouse_VISp_Barcodedalt_20200118\Mouse_VISp_Barcodedalt_20200118.log'
pd.matchLogFile(logFilePath);
set(pd, 'libraryName', 'Mouse_VISp_Barcodealt_from010120log');
pd.buildLibrary()


 















