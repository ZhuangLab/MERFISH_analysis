% Sum the raw signal within features in single FOV
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 21, 2017
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2018.
% -------------------------------------------------------------------------
% Purpose: 1) To sum the raw signals within features in a MERFISH dataset
% -------------------------------------------------------------------------

%% Load MERFISH Decoder object
if ~exist('arrayID')
    fovID = str2num(getenv('SLURM_ARRAY_TASK_ID')); % SLURM array id is the FOV
else
    fovID = arrayID; % Allow for alternative methods to set this fov id
end
nodeID = getenv('SLURM_NODELIST');
jobID = getenv('SLURM_JOBID');
PageBreak();
display(['Running SumFOV on  ' num2str(fovID) ' at ' datestr(now)]);
display(['Running on ' nodeID]);
display(['With job id ' jobID]);

mDecoder = MERFISHDecoder.Load(basePath);
mDecoder.overwrite = false; % Allow graceful restart of processes

if isempty(fovID)
	error;
end

mDecoder.SumRawSignalFOV(fovID);

PageBreak();
display(['Completed summation of raw data in fov ' num2str(fovID) ' at ' datestr(now)]);

