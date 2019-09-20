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

pd = probeDesign('MO4E4', 'mouse', 'D:\Data\MERFISH\Musmusculus\MO4E4\M1_codebook_AIBSFormat.csv');
set(pd, 'regionGC', [0.43, 0.63], 'regionTm', [66,76], 'isoSpecificity', [0.75, 1], 'specificity', [0.75, 1], 'probeSpacing', -20);
pd.buildLibrary();
