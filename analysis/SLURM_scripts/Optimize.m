% Optimize data sets on Odyssey
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 21, 2017
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2018.
%--------------------------------------------------------------------------
% Purpose: 1) Optimize a data set on Odyssey
% -------------------------------------------------------------------------

if ~exist('nDPath')
	error('A normalized data path must be provided.');
end
%% Create diary
if ~exist([nDPath 'log'], 'dir')
     mkdir([nDPath 'log']);
	 display(['Created ' nDPath 'log']); 
end
diary([nDPath 'log' filesep 'optimize.log']);
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

%% Aggregate warp data and generate warp report
mDecoder.GenerateWarpReport();

%% Combine pixel histograms and prepare initial scale factors
mDecoder.InitializeScaleFactors();

%% Setup parallel pool
% create a local cluster object
pc = parcluster('local')

% explicitly set the JobStorageLocation to the temp directory that was
% created in your sbatch script
pc.JobStorageLocation = strcat('/scratch/', getenv('USER'), '/', getenv('SLURM_JOB_ID'))

% start the parallel pool
p = parpool(pc, str2num(getenv('SLURM_NTASKS')))

%% Add parallel pool object to MERFISHDecoder
mDecoder.SetParallel(p);

%% Check to see if overwrite has been defined
if ~exist('overwrite') 
	overwrite = false;
end

mDecoder.OptimizeScaleFactors(25, 'overwrite', overwrite, 'useBlanks', false); % Normally 50 CHANGE ME BACK

%% Archive analysis
mDecoder.Save();

%% Archive analysis
PageBreak();
display(['...completed in ' num2str(toc(scriptTimer)) ' s']);
display(['Completed at ' datestr(now)]);

% Turn off diary
diary off;

%% Cleanup parallel pool
delete(p);

