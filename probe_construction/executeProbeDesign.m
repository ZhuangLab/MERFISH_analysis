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
  %  logFilePath = 'D:\Data\MERFISH\Homosapiens\Human_MTG_corrected3_barcoded_20200111b_actually sequential\Human_MTG_corrected3_barcoded_20200111b.log';
   % pd.matchLogFile(logFilePath);
    %set(pd, 'libraryName', 'Human_MTG_corrected3_sequential_20200121', 'species', 'Homo sapiens');
    %pd.buildLibrary();
    %display('Completed sequential');

%catch mError
 %   display('Error on Human_MTG_corrected3_sequential_20200121');
%end

%try 
 %   pd = probeDesign();
  %  logFilePath = 'D:\Data\MERFISH\Homosapiens\Human_MTG_barcoded_20200114\Human_MTG_barcoded_20200114.log';
   % pd.matchLogFile(logFilePath);
    %set(pd, 'libraryName', 'Human_MTG_barcoded_202000121', 'species', 'human');
    %pd.buildLibrary();
    %display('Completed barcoded');
%catch mError
 %   display('Error on Human_MTG_barcoded_202000121');
%end
% 
% pd = probeDesign();
% pd.matchLogFile('D:\Data\MERFISH\Musmusculus\Mouse_VISp_Barcodealt_from010120log\Mouse_VISp_Barcodealt_from010120log.log');
% set(pd, 'libraryName', 'MouseVISp_memErrorTesting_v2', 'specifyReadouts', true, 'keepAllPossibleProbes', true, 'debugMode', true);
% pd.buildLibrary()

bf = 'C:\Users\ScanningLabAnalysis\Documents\MATLAB\MERFISHProbeDesign\targetRegionsCalcs';
cdbkName = 'fakeMTGcodebook.csv';
pd = probeDesign();
pd.matchLogFile('D:\Data\MERFISH\Homosapiens\Human_MTG_barcoded_20200114\Human_MTG_barcoded_20200114.log');
set(pd, 'rawTranscriptomeFasta', 'D:\Data\MERFISH\Homosapiens\LIMSReferenceFiles\human.transcript.genesymbol.fa');
set(pd, 'ncRNAPath', 'D:\Data\MERFISH\Homosapiens\LIMSReferenceFiles\human.transcript.genesymbol.nc.fa');
set(pd, 'codebookPath', fullfile(bf, 'fakeMTGcodebook.csv'));
set(pd, 'fpkmPath', fullfile(bf, 'MTG.cleaned.isoforms.fpkm_tracking'));
set(pd, 'libraryName', 'HumanMTG_TargetRegionsTester_intTest', ...
        'specifyReadouts', true, ...
        'keepAllPossibleProbes', true, ...
        'debugMode', true, ...
        'useUniformWeights', false, ...
        'FPKMabundanceThreshold', 0.01, ...
        'transcriptomeHeaderType', 'cufflinks', ...
        'geneIsoformListSource', 'allGenes',...
        'tRFilterMethod', 'commonRegions', ...
        'spaceOutProbes', true);
pd.buildLibrary()

%%
pd = probeDesign();
pd.matchLogFile('D:\Data\MERFISH\Homosapiens\HumanMTG_TargetRegionsTester_intTest\HumanMTG_TargetRegionsTester_intTest.log');




 















