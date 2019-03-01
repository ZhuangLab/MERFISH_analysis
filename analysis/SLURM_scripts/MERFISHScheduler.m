% A master script for coordinating the many aspects of MERFISH analysis
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 21, 2017
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2018.
% -------------------------------------------------------------------------
% Purpose: 1) To act as a metascheduler, interfacing matlab with SLURM
% 2) To automate the many steps associated with analyzing a MERFISH dataset
% 3) To better handle the rare, unexplained, cluster failure
% -------------------------------------------------------------------------

%% ------------------------------------------------------------------------
% Prepare workspace 
% -------------------------------------------------------------------------
%% Check for information on what analysis to perform
if ~exist('aControl')
    % m = decoder; w = warp/preprocess; o = optimize; d = decode;
    % s = segment; c = combine boundaries; p = parse; 
    % f= perFormance; n = calculate Numbers;
    % r=sum Raw data; i = combIne raw sum data
    % l=low resolution mosaics
    % b=barcode metadata
    % u=doUblet score
    aControl = 'mwodscpfnrilbu';
end

%% Check for required paths
if ~exist('nDPath')
    error('A normalized data path must be provided');
end

if ismember('m', aControl) && ~exist('rDPath')
    error('A raw data path must be provided to create a MERFISHDecoder instance');
elseif ~ismember('m', aControl)
    rDPath = '';
end

if ismember('f', aControl) && ~exist('aPath')
    error('A path to abundance data must be provided to calculate MERFISH performance metrics');
end

%% Check for skip decoding construction parameter
if ~exist('skipDecoderConstruction')
    skipDecoderConstruction = false;
end

%% Check for queue information
if ~exist('smallJobQueue')
    %smallJobQueue = 'shared,general,zhuang,serial_requeue';
    smallJobQueue = 'shared,zhuang';
end
if ~exist('bigJobQueue')
    bigJobQueue = 'zhuang,general,shared';
end

%% Check for an email
if ~exist('email')
    warning('No email address was provided. So no status updates will be sent');
    email = '';
end

%% Check for provided parameters
if ~exist('parameters')
	parameters = [];
end

%% Check for matlab and python load commands
if ~exist('loadMatlabCommand', 'var')
	loadMatlabCommand = 'module load matlab/R2017a-fasrc02';
end
if ~exist('loadPythonCommand', 'var')
	loadPythonCommand = 'module load python/2.7.14-fasrc01';
end
if ~exist('activateEnvironmentCommand', 'var')
	activateEnvironmentCommand = ['source activate ' merfish_ENV];
end

%% Create normalized data path if it does not exist
if ~exist(nDPath, 'dir')
    mkdir(nDPath);
end

%% Create log path
logPath = [nDPath filesep 'log' filesep];
if ~exist(logPath, 'dir')
    mkdir(logPath);
end

%% Create scheduler path if it does not exist
schedulerPath = [nDPath filesep 'scheduler' filesep];
if ~exist(schedulerPath, 'dir')
    mkdir(schedulerPath);
end

%% Create and start analysis log
diary([nDPath 'log' filesep 'analysis_scheduler.log']);
diary on;

%% Archive start to log file
PageBreak();
disp(['Running MERFISH Scheduler: ']);
disp(['Raw data path: ' rDPath]);
disp(['Analyzed data path: ' nDPath]);
disp(['Requested analysis: ' aControl]);
disp(['Small jobs queue: ' smallJobQueue]);
disp(['Big jobs queue: ' bigJobQueue]);
disp(['Provided parameters: ']);
disp(parameters);
disp(['Started at ' datestr(now)]);
schedulerTimer = tic;

%% ------------------------------------------------------------------------
% Create messenger
% -------------------------------------------------------------------------
% Determine the name of messenger from the nDPath
pathParts = strsplit(nDPath, filesep);
pathParts = pathParts(~cellfun(@isempty, pathParts));

% Construct a morpheus instance
morpheus = Morpheus(email, 'name', pathParts{end}, 'maxNumErrors', 10);

%% ------------------------------------------------------------------------
% Create vrykolakas to keep the shell alive
% -------------------------------------------------------------------------
if ~exist('vryObj') % Only make one if another does not exist
    vryObj = Vrykolakas('timerPeriod', 20*60, ... % 20 minutes seems sufficient to overcome ssh timeouts
        'aliveString', ' ');
end

%% ------------------------------------------------------------------------
% Build decoder 
% -------------------------------------------------------------------------
if ismember('m', aControl) && ~skipDecoderConstruction
    % Confirm that the decoder does not yet exist
    if exist([nDPath filesep 'mDecoder'], 'dir')
        error('A decoder already exists.');
    end
   
    % Mark the start of decoder construction
    PageBreak();
    disp(['Creating MERFISHDecoder']);
    disp(['Started at ' datestr(now)]);
    localTimer = tic;
    
    % Build the decoder
    mDecoder = MERFISHDecoder(rDPath, nDPath, ...
        'parameters', parameters);
        
    % Save the decoder
    mDecoder.Save(); % Assume the default location

    % Mark completion
    disp(['...completed MERFISHDecoder construction in ' num2str(toc(localTimer)) ' s']);
end

%% Load Decoder: to determine the fov ids
mDecoder = MERFISHDecoder.Load(nDPath);

%% ------------------------------------------------------------------------
% Create jobs and job arrays 
% -------------------------------------------------------------------------
%% Warp and preprocess data: Create job array
if ismember('w', aControl)
    % Create local log path
    localLogPath = [nDPath filesep 'log' filesep 'process' filesep];
    if ~exist(localLogPath, 'dir')
        mkdir(localLogPath);
    end

    % Create job array
    for f=1:length(mDecoder.fovIDs)
        % Create job for each fov
        preProcessJobs(f) = SLURMJob({loadMatlabCommand, ...
            loadPythonCommand, ...
			activateEnvironmentCommand, ...
            ['matlab-default -nosplash -nodesktop -r ' ...
            '"basePath = ''' nDPath '''; arrayID = ' num2str(mDecoder.fovIDs(f)) '; ' ...
            'ProcessFOV; exit;"']}, ...
            'name', ['Preprocess task ' num2str(f) ' for fov ' num2str(mDecoder.fovIDs(f))], ...
            'scriptPath', [schedulerPath 'preprocess' filesep], ...
            'scriptName', ['w_f' num2str(mDecoder.fovIDs(f)) '.slurm'], ...
            'outputLog', [localLogPath 'fov_process_' num2str(mDecoder.fovIDs(f)) '.out'], ...
            'errorLog', [localLogPath 'fov_process_' num2str(mDecoder.fovIDs(f)) '.err'], ...
            'openMode', 'append', ...
            'constraint', 'holyib', ...
            'timeLimit', 12*60, ...
            'memoryLimit', 12000, ...
            'cpusPerTask', 1, ... 
            'partition', smallJobQueue, ...
            'exclude', 'holyzhuang01,holy2c14407', ... % Zhuang lab 'big job' nodes
            'timerPeriod', 5*60, ...
            'numberResubmit', 5, ...
            'completeFcn', @(~)mDecoder.CheckStatus(mDecoder.fovIDs(f), 'process'), ...
            'verbose', true); 
    end
    % Add to job array
    preprocessJobArray = SLURMJobArray(preProcessJobs, ...
        'name', 'preprocess', ...
        'verbose', true);
    
    % Link to morpheus
    morpheus.AddSuccessListener(preprocessJobArray, 'JobArrayComplete');
    morpheus.AddErrorListener(preprocessJobArray, 'JobArrayFailed');
    
end
    
%% Optimize thresholds: Create job
if ismember('o', aControl)    
    % Define command
    command = {loadMatlabCommand, ...                   % Load matlab
        'mkdir -p /scratch/$USER/$SLURM_JOB_ID', ...    % Make a local scratch directory for parallel job
        ['matlab-default -nosplash -nodesktop -r "overwrite = true; nDPath = ''' nDPath ...
        '''; Optimize;"'], ...                          % Define call to optimize
        'rm -rf /scratch/$USER/$SLURM_JOB_ID'};         % Cleanup scratch directory

    % Create job for the optimization process
    optimizeJob = SLURMJob(command, ...
        'name', 'Optimize', ...
        'scriptPath', [schedulerPath 'optimize' filesep], ...
        'scriptName', ['optimize.slurm'], ...
        'outputLog', [logPath 'optimize.out'], ...
        'errorLog', [logPath 'optimize.err'], ...
        'openMode', 'append', ...
        'constraint', 'holyib', ...
        'timeLimit', 24*60, ...
        'memoryLimit', 255000, ...
        'ntasks', 10, ...
        'partition', bigJobQueue, ...
        'timerPeriod', 20*60, ...
        'completeFcn', @(~)all(mDecoder.CheckStatus([], 'optimize')), ...
        'verbose', true); 
    
    % Add to job array (to allow triggered submit)
    optimizeJobArray = SLURMJobArray(optimizeJob, ...
        'name', 'optimize', ...
        'verbose', true);
    
    % Link to morpheus
    morpheus.AddSuccessListener(optimizeJobArray, 'JobArrayComplete');
    morpheus.AddErrorListener(optimizeJobArray, 'JobArrayFailed');
    
end

%% Create decoding jobs
if ismember('d', aControl)
    % Create local log path
    localLogPath = [nDPath filesep 'log' filesep 'decoding' filesep];
    if ~exist(localLogPath, 'dir')
        mkdir(localLogPath);
    end

    % Create job array
    for f=1:length(mDecoder.fovIDs)
        % Create job for each fov
        decodingJobs(f) = SLURMJob({loadMatlabCommand, ...
            ['matlab-default -nosplash -nodesktop -r ' ...
            '"basePath = ''' nDPath '''; arrayID = ' num2str(mDecoder.fovIDs(f)) '; ' ...
            'DecodeFOV; exit;"']}, ...
            'name', ['Decoding task ' num2str(f) ' for fov ' num2str(mDecoder.fovIDs(f))], ...
            'scriptPath', [schedulerPath 'decoding' filesep], ...
            'scriptName', ['d_f' num2str(mDecoder.fovIDs(f)) '.slurm'], ...
            'outputLog', [localLogPath 'fov_decode_' num2str(mDecoder.fovIDs(f)) '.out'], ...
            'errorLog', [localLogPath 'fov_decode_' num2str(mDecoder.fovIDs(f)) '.err'], ...
            'openMode', 'append', ...
            'constraint', 'holyib', ...
            'timeLimit', 6*60, ...
            'memoryLimit', 32000, ...
            'cpusPerTask', 1, ... 
            'partition', smallJobQueue, ...
            'exclude', 'holyzhuang01,holy2c14407', ... % Zhuang lab 'big job' nodes
            'timerPeriod', 5*60, ...
            'numberResubmit', 5, ...
            'completeFcn', @(~)mDecoder.CheckStatus(mDecoder.fovIDs(f), 'decode'), ...
            'verbose', true);  
    end
    
    % Combine into job array
    decodingJobArray = SLURMJobArray(decodingJobs, ...
        'name', 'decoding', ...
        'verbose', true);
    
    % Link to morpheus
    morpheus.AddSuccessListener(decodingJobArray, 'JobArrayComplete');
    morpheus.AddErrorListener(decodingJobArray, 'JobArrayFailed');

end

%% Create segmentation jobs
if ismember('s', aControl)
    % Create local log path
    localLogPath = [nDPath filesep 'log' filesep 'segment' filesep];
    if ~exist(localLogPath, 'dir')
        mkdir(localLogPath);
    end

    % Create job array
    for f=1:length(mDecoder.fovIDs)
        % Create job for each fov
        segmentJobs(f) = SLURMJob({loadMatlabCommand, ...
            ['matlab-default -nosplash -nodesktop -r ' ...
            '"basePath = ''' nDPath '''; arrayID = ' num2str(mDecoder.fovIDs(f)) '; ' ...
            'SegmentFOV; exit;"']}, ...
            'name', ['Segmentation task ' num2str(f) ' for fov ' num2str(mDecoder.fovIDs(f))], ...
            'scriptPath', [schedulerPath 'segmentation' filesep], ...
            'scriptName', ['s_f' num2str(mDecoder.fovIDs(f)) '.slurm'], ...
            'outputLog', [localLogPath 'fov_segment_' num2str(mDecoder.fovIDs(f)) '.out'], ...
            'errorLog', [localLogPath 'fov_segment_' num2str(mDecoder.fovIDs(f)) '.err'], ...
            'openMode', 'append', ...
            'constraint', 'holyib', ...
            'timeLimit', 30, ...
            'memoryLimit', 8000, ...
            'cpusPerTask', 1, ... 
            'partition', smallJobQueue, ...
            'exclude', 'holyzhuang01,holy2c14407', ... % Zhuang lab 'big job' nodes
            'timerPeriod', 15*60, ...
            'numberResubmit', 5, ...
            'completeFcn', @(~)mDecoder.CheckStatus(mDecoder.fovIDs(f), 'segment'), ...
            'verbose', true); 
    end
    
    % Combine into job array
    segmentJobArray = SLURMJobArray(segmentJobs, ...
        'name', 'segment', ...
        'verbose', true);
    
    % Link to morpheus
    morpheus.AddSuccessListener(segmentJobArray, 'JobArrayComplete');
    morpheus.AddErrorListener(segmentJobArray, 'JobArrayFailed');
    
end

%% Create combine features job
if ismember('c', aControl)    
    % Define command
    command = {loadMatlabCommand, ... % Load matlab
        ['matlab-default -nosplash -nodesktop -r "overwrite = true; nDPath = ''' nDPath ...
        '''; CombineFoundFeatures();"']};         % Issue combine boundaries command

    % Create job for the combination process
    combineJob = SLURMJob(command, ...
        'name', 'Combine', ...
        'scriptPath', [schedulerPath 'combine' filesep], ...
        'scriptName', ['combine.slurm'], ...
        'outputLog', [logPath 'combine.out'], ...
        'errorLog', [logPath 'combine.err'], ...
        'openMode', 'append', ...
        'constraint', 'holyib', ...
        'timeLimit', 12*60, ...
        'memoryLimit', 64000, ...
        'ntasks', 1, ... 
        'partition', bigJobQueue, ...
        'timerPeriod', 5*60, ...
        'completeFcn', @(~)all(mDecoder.CheckStatus([], 'combine')), ...
        'verbose', true); 
    
    % Add to job array (to allow triggered submit)
    combineJobArray = SLURMJobArray(combineJob, ...
        'name', 'combine', ...
        'verbose', true);
    
    % Link to morpheus
    morpheus.AddSuccessListener(combineJobArray, 'JobArrayComplete');
    morpheus.AddErrorListener(combineJobArray, 'JobArrayFailed');
    
end

%% Create low resolution mosaic job
if ismember('l', aControl)    
    % Define command
    command = {loadMatlabCommand, ... % Load matlab
        ['matlab-default -nosplash -nodesktop -r "nDPath = ''' nDPath ...
        '''; LowResMosaic();"']};         % Issue low resolution mosaic command

    % Create job for the combination process
    mosaicJob = SLURMJob(command, ...
        'name', 'Mosaic', ...
        'scriptPath', [schedulerPath 'mosaic' filesep], ...
        'scriptName', ['mosaic.slurm'], ...
        'outputLog', [logPath 'mosaic.out'], ...
        'errorLog', [logPath 'mosaic.err'], ...
        'openMode', 'append', ...
        'constraint', 'holyib', ...
        'timeLimit', 12*60, ...
        'memoryLimit', 127000, ...
        'ntasks', 1, ... 
        'partition', bigJobQueue, ...
        'timerPeriod', 5*60, ...
        'verbose', true, ...
        'numberResubmit', 3, ...
        'completeFcn', @(~)all(mDecoder.CheckStatus(1, 'mosaic')), ...
        'preCheck', true);
    
    % Add to job array (to allow triggered submit)
    mosaicJobArray = SLURMJobArray(mosaicJob, ...
        'name', 'mosaic', ...
        'verbose', true);
    
    % Link to morpheus
    morpheus.AddSuccessListener(mosaicJobArray, 'JobArrayComplete');
    morpheus.AddErrorListener(mosaicJobArray, 'JobArrayFailed');
    
end

%% Create parse jobs
if ismember('p', aControl)
    % Create local log path
    localLogPath = [nDPath filesep 'log' filesep 'parse' filesep];
    if ~exist(localLogPath, 'dir')
        mkdir(localLogPath);
    end
    
    % Create job array
    for f=1:length(mDecoder.fovIDs)
        % Create job for each fov
        parseJobs(f) = SLURMJob({loadMatlabCommand, ...
            ['matlab-default -nosplash -nodesktop -r ' ...
            '"basePath = ''' nDPath '''; arrayID = ' num2str(mDecoder.fovIDs(f)) '; ' ...
            'ParseFOV; exit;"']}, ...
            'name', ['Parse task ' num2str(f) ' for fov ' num2str(mDecoder.fovIDs(f))], ...
            'scriptPath', [schedulerPath 'parse' filesep], ...
            'scriptName', ['p_f' num2str(mDecoder.fovIDs(f)) '.slurm'], ...
            'outputLog', [localLogPath 'fov_parse_' num2str(mDecoder.fovIDs(f)) '.out'], ...
            'errorLog', [localLogPath 'fov_parse_' num2str(mDecoder.fovIDs(f)) '.err'], ...
            'openMode', 'append', ...
            'constraint', 'holyib', ...
            'timeLimit', 3*60, ...
            'memoryLimit', 24000, ...
            'cpusPerTask', 1, ... 
            'partition', smallJobQueue, ...
            'exclude', 'holyzhuang01,holy2c14407', ... % Zhuang lab 'big job' nodes
            'timerPeriod', 5*60, ...
            'numberResubmit', 5, ...
            'completeFcn', @(~)mDecoder.CheckStatus(mDecoder.fovIDs(f), 'parse'), ...
            'verbose', true); 
    end
    % Combine into job array
    parseJobArray = SLURMJobArray(parseJobs, ...
        'name', 'parse', ...
        'verbose', true);
        
    % Link to morpheus
    morpheus.AddSuccessListener(parseJobArray, 'JobArrayComplete');
    morpheus.AddErrorListener(parseJobArray, 'JobArrayFailed');
    
end

%% Calculate performance
if ismember('f', aControl)    
    % Define command
    command = {loadMatlabCommand, ...                   % Load matlab
        'mkdir -p /scratch/$USER/$SLURM_JOB_ID', ...    % Make a local scratch directory for parallel job
        ['matlab-default -nosplash -nodesktop -r "overwrite = true; nDPath = ''' nDPath '''; ' ....
        'aPath = ''' aPath '''; ' ...
        'Performance;"'], ...                           % Define call to optimize
        'rm -rf /scratch/$USER/$SLURM_JOB_ID'};         % Cleanup scratch directory

    % Create job for the optimization process
    performanceJob = SLURMJob(command, ...
        'name', 'Performance', ...
        'scriptPath', [schedulerPath 'performance' filesep], ...
        'scriptName', ['performance.slurm'], ...
        'outputLog', [logPath 'performance.out'], ...
        'errorLog', [logPath 'performance.err'], ...
        'openMode', 'truncate', ...
        'constraint', 'holyib', ...
        'timeLimit', 4*60, ...
        'memoryLimit', 255000, ...
        'ntasks', 12, ... 
        'partition', bigJobQueue, ...
        'timerPeriod', 5*60, ...
        'verbose', true, ...
        'preCheck', false); 
    
    % Add to job array (to allow triggered submit)
    performanceJobArray = SLURMJobArray(performanceJob, ...
        'name', 'performance', ...
        'verbose', true);
    
    % Link to morpheus
    morpheus.AddSuccessListener(performanceJobArray, 'JobArrayComplete');
    morpheus.AddErrorListener(performanceJobArray, 'JobArrayFailed');
    
end

%% Calculate numbers
if ismember('n', aControl)    
    % Define command
    command = {loadMatlabCommand, ...                   % Load matlab
        'mkdir -p /scratch/$USER/$SLURM_JOB_ID', ...    % Make a local scratch directory for parallel job
        ['matlab-default -nosplash -nodesktop -r "overwrite = true; nDPath = ''' nDPath ...
        '''; CalculateNumbers();"'], ...                          % Define call to optimize
        'rm -rf /scratch/$USER/$SLURM_JOB_ID'};         % Cleanup scratch directory

    % Create job for the optimization process
    numbersJob = SLURMJob(command, ...
        'name', 'Numbers', ...
        'scriptPath', [schedulerPath 'numbers' filesep], ...
        'scriptName', ['numbers.slurm'], ...
        'outputLog', [logPath 'numbers.out'], ...
        'errorLog', [logPath 'numbers.err'], ...
        'openMode', 'truncate', ...
        'constraint', 'holyib', ...
        'timeLimit', 10*60, ...
        'memoryLimit', 128000, ...
        'ntasks', 12, ... 
        'partition', bigJobQueue, ...
        'timerPeriod', 5*60, ...
        'verbose', true, ...
		'preCheck', false);
    
    % Add to job array (to allow triggered submit)
    numbersJobArray = SLURMJobArray(numbersJob, ...
        'name', 'numbers', ...
        'verbose', true);
    
    % Link to morpheus
    morpheus.AddSuccessListener(numbersJobArray, 'JobArrayComplete');
    morpheus.AddErrorListener(numbersJobArray, 'JobArrayFailed');
    
end

%% Export barcode metadata
if ismember('b', aControl)    
    % Define command
    command = {loadMatlabCommand, ...                   % Load matlab
        'mkdir -p /scratch/$USER/$SLURM_JOB_ID', ...    % Make a local scratch directory for parallel job
        ['matlab-default -nosplash -nodesktop -r "overwrite = true; nDPath = ''' nDPath ...
        '''; ExportBarcodeMetadata();"'], ...                          % Define call to optimize
        'rm -rf /scratch/$USER/$SLURM_JOB_ID'};         % Cleanup scratch directory

    % Create job for the optimization process
    metadataJob = SLURMJob(command, ...
        'name', 'Export Barcodes', ...
        'scriptPath', [schedulerPath 'barcode_metadata' filesep], ...
        'scriptName', ['barcode_metadata.slurm'], ...
        'outputLog', [logPath 'barcode_metadata.out'], ...
        'errorLog', [logPath 'barcode_metadata.err'], ...
        'openMode', 'truncate', ...
        'constraint', 'holyib', ...
        'timeLimit', 10*60, ...
        'memoryLimit', 250000, ...
        'ntasks', 12, ... 
        'partition', bigJobQueue, ...
        'timerPeriod', 5*60, ...
        'verbose', true, ...
		'preCheck', false);
    
    % Add to job array (to allow triggered submit)
    metadataJobArray = SLURMJobArray(metadataJob, ...
        'name', 'export barcodes', ...
        'verbose', true);
    
    % Link to morpheus
    morpheus.AddSuccessListener(metadataJobArray, 'JobArrayComplete');
    morpheus.AddErrorListener(metadataJobArray, 'JobArrayFailed');
    
end

%% Calculate doublet score values
if ismember('u', aControl)    
    % Define command
    command = {loadMatlabCommand, ...                   % Load matlab
        ['matlab-default -nosplash -nodesktop -r "overwrite = true; nDPath = ''' nDPath ...
        '''; CalculateDoubletScore();"']};         % Cleanup scratch directory

    % Create job for the optimization process
    doubletScoreJob = SLURMJob(command, ...
        'name', 'Calculate Doublet Score', ...
        'scriptPath', [schedulerPath 'doublet_score' filesep], ...
        'scriptName', ['doublet_score.slurm'], ...
        'outputLog', [logPath 'doublet_score.out'], ...
        'errorLog', [logPath 'doublet_score.err'], ...
        'openMode', 'truncate', ...
        'constraint', 'holyib', ...
        'timeLimit', 10*60, ...
        'memoryLimit', 32000, ...
        'ntasks', 1, ... 
        'partition', bigJobQueue, ...
        'timerPeriod', 5*60, ...
        'verbose', true, ...
		'preCheck', false);
    
    % Add to job array (to allow triggered submit)
    doubletScoreJobArray = SLURMJobArray(doubletScoreJob, ...
        'name', 'doublet score', ...
        'verbose', true);
    
    % Link to morpheus
    morpheus.AddSuccessListener(doubletScoreJobArray, 'JobArrayComplete');
    morpheus.AddErrorListener(doubletScoreJobArray, 'JobArrayFailed');
    
end


%% Sum raw data: Create job array
if ismember('r', aControl)
    % Create local log path
    localLogPath = [nDPath filesep 'log' filesep 'sum' filesep];
    if ~exist(localLogPath, 'dir')
        mkdir(localLogPath);
    end
    
    % Create job array
    for f=1:length(mDecoder.fovIDs)
        % Create job for each fov
        sumJobs(f) = SLURMJob({loadMatlabCommand, ...
            ['matlab-default -nosplash -nodesktop -r ' ...
            '"basePath = ''' nDPath '''; arrayID = ' num2str(mDecoder.fovIDs(f)) '; ' ...
            'SumFOV; exit;"']}, ...
            'name', ['Sum task ' num2str(f) ' for fov ' num2str(mDecoder.fovIDs(f))], ...
            'scriptPath', [schedulerPath 'sum' filesep], ...
            'scriptName', ['r_f' num2str(mDecoder.fovIDs(f)) '.slurm'], ...
            'outputLog', [localLogPath 'fov_sum_' num2str(mDecoder.fovIDs(f)) '.out'], ...
            'errorLog', [localLogPath 'fov_sum_' num2str(mDecoder.fovIDs(f)) '.err'], ...
            'openMode', 'append', ...
            'constraint', 'holyib', ...
            'timeLimit', 3*60, ...
            'memoryLimit', 24000, ...
            'cpusPerTask', 1, ... 
            'partition', smallJobQueue, ...
            'exclude', 'holyzhuang01,holy2c14407', ... % Zhuang lab 'big job' nodes
            'timerPeriod', 5*60, ...
            'numberResubmit', 5, ...
            'completeFcn', @(~)mDecoder.CheckStatus(mDecoder.fovIDs(f), 'sum'), ...
            'verbose', true); 
    end
    % Add to job array
    sumJobArray = SLURMJobArray(sumJobs, ...
        'name', 'sum', ...
        'verbose', true);
    
    % Link to morpheus
    morpheus.AddSuccessListener(sumJobArray, 'JobArrayComplete');
    morpheus.AddErrorListener(sumJobArray, 'JobArrayFailed');
    
end

%% Combine raw sum: Create job array
if ismember('i', aControl)    
    % Define command
    command = {loadMatlabCommand, ...   % Load matlab
        ['matlab-default -nosplash -nodesktop -r "overwrite = true; nDPath = ''' nDPath ...
        '''; CombineSum();"']};         % Issue combine boundaries command

    % Create job for the optimization process
    combineSumJob = SLURMJob(command, ...
        'name', 'Combine sum', ...
        'scriptPath', [schedulerPath 'combine_sum' filesep], ...
        'scriptName', ['combine_sum.slurm'], ...
        'outputLog', [logPath 'combine_sum.out'], ...
        'errorLog', [logPath 'combine_sum.err'], ...
        'openMode', 'append', ...
        'constraint', 'holyib', ...
        'timeLimit', 3*60, ...
        'memoryLimit', 64000, ...
        'ntasks', 1, ... 
        'partition', bigJobQueue, ...
        'timerPeriod', 5*60, ...
        'completeFcn', @(~)all(mDecoder.CheckStatus([], 'combine_sum')), ...
        'verbose', true); 
    
    % Add to job array (to allow triggered submit)
    combineSumJobArray = SLURMJobArray(combineSumJob, ...
        'name', 'combine sum', ...
        'verbose', true);
    
    % Link to morpheus
    morpheus.AddSuccessListener(combineSumJobArray, 'JobArrayComplete');
    morpheus.AddErrorListener(combineSumJobArray, 'JobArrayFailed');
    
end

%% ------------------------------------------------------------------------
% Submit jobs and job arrays 
% -------------------------------------------------------------------------
%% Submit performance job array
if ismember('f', aControl)
    % Configure submission dependencies
    if ismember('d', aControl)
        performanceJobArray.AddSubmitTrigger(decodingJobArray, 'JobArrayComplete');
    else
        performanceJobArray.Submit();
    end
end

%% Submit doublet score job
if ismember('u', aControl)
    % Configure submission dependencies
    if ismember('b', aControl)
        doubletScoreJobArray.AddSubmitTrigger(metadataJobArray, 'JobArrayComplete');
    else
        doubletScoreJobArray.Submit();
    end
end

%% Submit barcode metadata export job
if ismember('b', aControl)
    % Configure submission dependencies
    if ismember('p', aControl)
        metadataJobArray.AddSubmitTrigger(parseJobArray, 'JobArrayComplete');
    else
        metadataJobArray.Submit();
    end
end

%% Submit numbers job array
if ismember('n', aControl)
    % Configure submission dependencies
    if ismember('p', aControl)
        numbersJobArray.AddSubmitTrigger(parseJobArray, 'JobArrayComplete');
    else
        numbersJobArray.Submit();
    end
end

%% Submit combine raw sum job array
if ismember('i', aControl)
    % Configure submission dependencies
    if ismember('r', aControl) && ismember('l', aControl)
        combineSumJobArray.AddSubmitTrigger(mosaicJobArray, 'JobArrayComplete');
        combineSumJobArray.AddSubmitTrigger(sumJobArray, 'JobArrayComplete');
    elseif ismember('r', aControl) && ~ismember('l', aControl)
        combineSumJobArray.AddSubmitTrigger(sumJobArray, 'JobArrayComplete');
    elseif ~ismember('r', aControl) && ismember('l', aControl)
        combineSumJobArray.AddSubmitTrigger(mosaicJobArray, 'JobArrayComplete');
    else
        combineSumJobArray.Submit();
    end
end

%% Submit raw sum job array
if ismember('r', aControl)
    % Configure submission dependencies
    if ismember('c', aControl)
        sumJobArray.AddSubmitTrigger(combineJobArray, 'JobArrayComplete');
    else
        sumJobArray.Submit();
    end
end

%% Submit parse job array
if ismember('p', aControl)
    % Configure submission dependencies
    if ismember('c', aControl) && ismember('d', aControl)
        parseJobArray.AddSubmitTrigger(combineJobArray, 'JobArrayComplete');
        parseJobArray.AddSubmitTrigger(decodingJobArray, 'JobArrayComplete');
    elseif ismember('c', aControl) && ~ismember('d', aControl)
        parseJobArray.AddSubmitTrigger(combineJobArray, 'JobArrayComplete');
    elseif ~ismember('c', aControl) && ismember('d', aControl)
        parseJobArray.AddSubmitTrigger(decodingJobArray, 'JobArrayComplete');
    else
        parseJobArray.Submit();
    end
end

%% Submit combine features job array
if ismember('c', aControl)
    % Configure submission dependencies
    if ismember('s', aControl) && ismember('l', aControl)
        combineJobArray.AddSubmitTrigger(mosaicJobArray, 'JobArrayComplete');
        combineJobArray.AddSubmitTrigger(segmentJobArray, 'JobArrayComplete');
    elseif ~ismember('s', aControl) && ismember('l', aControl)
        combineJobArray.AddSubmitTrigger(mosaicJobArray, 'JobArrayComplete');
    elseif ismember('s', aControl) && ~ismember('l', aControl)
        combineJobArray.AddSubmitTrigger(segmentJobArray, 'JobArrayComplete');
    else
        combineJobArray.Submit();
    end
end

%% Submit segment job array
if ismember('s', aControl)
    % Configure submission dependencies
    if ismember('w', aControl)
        segmentJobArray.AddSubmitTrigger(preprocessJobArray, 'JobArrayComplete');
    else
        segmentJobArray.Submit();
    end
end

%% Submit decoding job array
if ismember('d', aControl)
    % Configure submission dependencies
    if ismember('o', aControl)
        decodingJobArray.AddSubmitTrigger(optimizeJobArray, 'JobArrayComplete');
    else
        decodingJobArray.Submit();
    end
end

%% Submit optimize job array
if ismember('o', aControl)
    % Configure submission dependencies
    if ismember('w', aControl)
        optimizeJobArray.AddSubmitTrigger(preprocessJobArray, 'JobArrayComplete');
    else
        optimizeJobArray.Submit();
    end
end

%% Submit mosaic job array
if ismember('l', aControl)
    % Configure submission dependencies
    if ismember('w', aControl)
        mosaicJobArray.AddSubmitTrigger(preprocessJobArray, 'JobArrayComplete');
    else
        mosaicJobArray.Submit();
    end
end

%% Submit preprocess job array
if ismember('w', aControl)
    preprocessJobArray.Submit();
end
