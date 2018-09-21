% Segment all FOV
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 21, 2017
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2018.
% -------------------------------------------------------------------------
% Purpose: 1) To segment a MERFISH data set
% -------------------------------------------------------------------------

if ~exist('nDPath')
	error('A normalized data path must be provided.');
end

%% Record properties of the run
nodeID = getenv('SLURM_NODELIST');
jobID = getenv('SLURM_JOBID');
PageBreak();
display(['Running ' mfilename ' at ' datestr(now)]);
display(['Running on ' nodeID]);
display(['With job id ' jobID]);

%% Load MERFISH Decoder
mDecoder = MERFISHDecoder.Load(nDPath);

%% Setup parallel pool
% create a local cluster object
pc = parcluster('local');

% explicitly set the JobStorageLocation to the temp directory that was
% created in your sbatch script
pc.JobStorageLocation = strcat('/scratch/jeffmoffitt/', getenv('SLURM_JOB_ID'));

% start the parallel pool
p = parpool(pc, str2num(getenv('SLURM_NTASKS')))

%% Add parallel pool object to MERFISHDecoder
mDecoder.SetParallel(p);

%% Segment all FOV
mDecoder.SegmentFOV([]);

%% Summarize progress
PageBreak();
display(['Completed decoding of fov ' num2str(fovID) ' at ' datestr(now)]);

