% Segment a FOV
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 21, 2017
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2018.
% -------------------------------------------------------------------------
% Purpose: 1) To segment a MERFISH data set
% -------------------------------------------------------------------------


%% Record properties of the run
if ~exist('arrayID')
    fovID = str2num(getenv('SLURM_ARRAY_TASK_ID')); % SLURM array id is the FOV
else
    fovID = arrayID; % Allow for alternative methods to set this fov id
end
nodeID = getenv('SLURM_NODELIST');
jobID = getenv('SLURM_JOBID');
PageBreak();
display(['Running ' mfilename ' on  ' num2str(fovID) ' at ' datestr(now)]);
display(['Running on ' nodeID]);
display(['With job id ' jobID]);

%% Confirm FOV info
if isempty(fovID)
	error;
end

%% Load MERFISH Decoder
mDecoder = MERFISHDecoder.Load(basePath);
mDecoder.overwrite = false; % Allow graceful restart of processes

%% Decode FOV
mDecoder.SegmentFOV(fovID);

%% Summarize progress
PageBreak();
display(['Completed segmentation of fov ' num2str(fovID) ' at ' datestr(now)]);

