% Combine Segmentation Boundaries
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 21, 2017
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2018.
% -------------------------------------------------------------------------
% Purpose: 1) Combined found features on Odyssey
% -------------------------------------------------------------------------

if ~exist('nDPath')
	error('A normalized data path must be provided.');
end

%% Get slurm properties
PageBreak();
nodeID = getenv('SLURM_NODELIST');

display(['Running on ' nodeID]);
display(['Started at ' datestr(now)]);
scriptTimer = tic;

%% Load MERFISH Decoder
display(['Loading MERFISH Decoder from ' nDPath]);
mDecoder = MERFISHDecoder.Load(nDPath);

%% Check to see if overwrite has been defined
if ~exist('overwrite') 
	overwrite = false;
end

%% Control for possible overwrite
mDecoder.overwrite = overwrite;

%% Issue command to combine features
mDecoder.CombineFeatures();

%% Generate report
mDecoder.GenerateFoundFeaturesReport()

%% Export the found features to a csv file
mDecoder.FoundFeaturesToCSV('downSampleFactor', 10, 'zIndex', 4);

%% Archive analysis
PageBreak();
display(['...completed in ' num2str(toc(scriptTimer)) ' s']);
display(['Completed at ' datestr(now)]);
