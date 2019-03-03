%%-------------------------------------------------------------------------
%% Script for running MERFISH analysis
%%-------------------------------------------------------------------------

%% initialize data paths
repositoryPath              = '/n/regal/zhuang_lab/lsepulveda';                                 
experimentPath              = '181029_BC037_MERFISH';
samplePath                  = 'sample_02';
dataPath                    = [repositoryPath '/' experimentPath];
settingsPath                = [dataPath '/Settings/' samplePath];

%% initialize merfish analysis paths. This variables are harcoded, variable names should not be changed
rDPath                      = [dataPath '/' samplePath ];                                       % raw data path
nDPath                      = [dataPath '/normalized_data/' samplePath '-new-pip-test/' ];      % normalize data path, 
aPath                       = ['~/Software/code-matlab/merfish/fpkm/FPKMDataPublished.matb'];   % U2OS Bulk RNAseq data, from Moffit et al, PNAS, 2016 
positionsPath               = [settingsPath '/positions_full-1000_1024_60x_MERFISH.txt'];       % positions file used to acquire images, each column is X,Y
dataOrganizationPath        = [settingsPath '/data_organization.csv'];                          % csv file describing the organization of image files (filenames, number of frames, etc)
codebookPath                = [settingsPath '/L26E1_codebook.csv'];                             % mapping between gene names and barcodes
scratchPath                 = [repositoryPath '/scratch/lsepulvedaduran'];                      % location of scratch directory

%% define name of merfish python environment
merfish_ENV                 = 'merfish_analysis';

%% Initialize file to save the output from merfish run
% Create diary
if ~exist([nDPath 'log'], 'dir')
     mkdir([nDPath 'log']);
     disp(['Created ' nDPath 'log']);
end
diary([nDPath 'log/mDecoder.log']);
diary on;

% Archive analysis
PageBreak();
display(['Creating MERFISHDecoder for ' rDPath]);
display(['Normalized data will be saved in ' nDPath]);
display(['Started at ' datestr(now)]);
scriptTimer = tic;

%% Initialize parameters for the decoder
parameters                      = [];
parameters.pixelSize            = 107.4;
parameters.imageSize            = [2048 2048];
parameters.lowPassKernelSize    = 1;
parameters.crop                 = 40;                          
parameters.minBrightness        = 1;                            % No threshold applied
parameters.minArea              = 1;                            % Minimum area for saving barcodes
parameters.areaThresh           = 4;                            % Area threshold for optimization
parameters.minBrightness        = 1;                            % All decoded barcodes will be saved...no initial cuts
parameters.stageOrientation     = [-1 1];                       % This value may not be correct....
parameters.overwrite            = false;
parameters.codebookPath         = codebookPath;                 % path to the codebook file
parameters.dataOrganizationPath = dataOrganizationPath;         % Path to a data organization file
parameters.hal_version          = 'hal2';  

%% Run the merfish scheduler
% define which analysis to perform:
% m = decoder;   w = warp/preprocess;    o = optimize;             d = decode;   
% s = segment;   c = combine boundaries; p = parse;                f = perFormance; 
% n = calc Nums; r = sum Raw data;       i = combIne raw sum data; l =low resolution mosaics
% b = barcode metadata; u = doUblet score
aControl = 'mwodscpfnlbu';
MERFISHScheduler;

%% Archive analysis
PageBreak();
disp(['...completed in ' num2str(toc(scriptTimer)) ' s']);
disp(['Completed at ' datestr(now)]);


