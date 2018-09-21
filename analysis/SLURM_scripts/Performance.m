% Calculate MERFISH performance on Odyssey
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 21, 2017
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2018.
% -------------------------------------------------------------------------
% Purpose: 1) Run the basic MERFISH performance metrics function
% -------------------------------------------------------------------------

if ~exist('nDPath')
	error('A normalized data path must be provided.');
end
if ~exist('aPath')
	error('A normalized data path must be provided.');
end

%% Create diary
if ~exist([nDPath 'log'], 'dir')
     mkdir([nDPath 'log']);
	 display(['Created ' nDPath 'log']); 
end
diary([nDPath 'log' filesep 'performance.log']);
diary on;

%% Get slurm properties
PageBreak();
nodeID = getenv('SLURM_NODELIST');

display(['Running on ' nodeID]);
display(['Started at ' datestr(now)]);
scriptTimer = tic;

%% Load MERFISH Decoder
display(['Loading MERFISH Decoder from ' nDPath]);
mDecoder = MERFISHDecoder.Load(nDPath);

%% Confirm validity of all barcode files
bFiles = BuildFileStructure([nDPath filesep 'barcodes' filesep 'barcode_fov' filesep], ...
    'regExp', ['fov_(?<fov>[0-9]+)_blist'], ...
    'fileExt', 'bin', ...
    'fieldNames', {'fov'}, ...
    'fieldConv', {@str2num});

display(['Found ' num2str(length(bFiles)) ' barcode files']);

isCorrupt = false(1, length(bFiles));
for b=1:length(bFiles)
    try
        header = ReadBinaryFileHeader(bFiles(b).filePath);
    catch
        isCorrupt(b) = header.isCorrupt;
    end
    if header.isCorrupt
        display(['Barcode file for fov ' num2str(bFiles(b).fov) ' is corrupt']);
    end
end
if all(~isCorrupt)
    display(['All barcode files appear to be complete and uncorrupted!']);
end

% Determine if some barcode files are missing
missingFovIds = setdiff(mDecoder.fovIDs, [bFiles.fov]);
if ~isempty(missingFovIds)
    warning('Discovered missing fov ids!');
    display(missingFovIds);
end

%% Setup parallel pool
% create a local cluster object
pc = parcluster('local');

% explicitly set the JobStorageLocation to the temp directory that was
% created in your sbatch script
pc.JobStorageLocation = strcat('/scratch/jeffmoffitt/', getenv('SLURM_JOB_ID'))

% start the parallel pool
p = parpool(pc, str2num(getenv('SLURM_NTASKS')))

%% Create parameters object
performanceParameters = [];
performanceParameters.parallel = p;
performanceParameters.brightnessThreshold = mDecoder.parameters.quantification.minimumBarcodeBrightness;
performanceParameters.areaThreshold = mDecoder.parameters.quantification.minimumBarcodeArea;
performanceParameters.stageOrientation = mDecoder.parameters.decoding.stageOrientation; 
performanceParameters.abundDataPath = aPath;
performanceParameters.verbose = true;
performanceParameters.outputPath = [mDecoder.normalizedDataPath mDecoder.reportPath 'performance' filesep];

%% Archive parameters in log
PageBreak();
disp(['Using the following parameters']);
foundFields = setdiff(fields(performanceParameters), 'parallel');
for f=1:length(foundFields)
    if ischar(performanceParameters.(foundFields{f}))
        disp(['   ' foundFields{f} ': ' num2str(performanceParameters.(foundFields{f}))]);
    else
        disp(['   ' foundFields{f} ': ' num2str(performanceParameters.(foundFields{f}))]);
    end
end

%% Calculate peformance
MERFISHPerformanceMetrics(nDPath, 'parameters', performanceParameters);

%% Archive analysis
PageBreak();
display(['...completed in ' num2str(toc(scriptTimer)) ' s']);
display(['Completed at ' datestr(now)]);

% Turn off diary
diary off;

%% Cleanup parallel pool
delete(p);

