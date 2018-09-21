% Export barcode metadata
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% October 19, 2017
% -------------------------------------------------------------------------
% Purpose: 1) Export a list of barcodes with limited metadata
% -------------------------------------------------------------------------
% This work is licenses under CC BY NC SA 
% -------------------------------------------------------------------------
if ~exist('nDPath')
	error('A normalized data path must be provided.');
end
%% Create diary
if ~exist([nDPath 'log'], 'dir')
     mkdir([nDPath 'log']);
	 display(['Created ' nDPath 'log']); 
end
diary([nDPath 'log' filesep 'calculate.log']);
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

%% Setup parallel pool
% Create a local cluster object
pc = parcluster('local')

% Explicitly set the JobStorageLocation to the temp directory that was
% created in your sbatch script
pc.JobStorageLocation = strcat('/scratch/jeffmoffitt/', getenv('SLURM_JOB_ID'))

% Start the parallel pool and set it to never timeout
p = parpool(pc, str2num(getenv('SLURM_NTASKS')), 'IdleTimeout', Inf)

%% Add parallel pool object to MERFISHDecoder
mDecoder.SetParallel(p);

%% Calculate feature counts
mDecoder.BarcodesToCSV();

%% Calculate the doublet score metadata features
CalculateDoubletScore();

%% Archive analysis
PageBreak();
display(['...completed in ' num2str(toc(scriptTimer)) ' s']);
display(['Completed at ' datestr(now)]);

% Turn off diary
diary off;

%% Cleanup parallel pool
delete(p);

