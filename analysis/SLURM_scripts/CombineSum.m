% Combine Raw Sum
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 21, 2017
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2018.
% -------------------------------------------------------------------------
% Purpose: 1) Combined summation boundaries on Odyssey
% -------------------------------------------------------------------------


if ~exist('nDPath')
	error('A normalized data path must be provided.');
end

%% Create diary
if ~exist([nDPath 'log'], 'dir')
     mkdir([nDPath 'log']);
	 display(['Created ' nDPath 'log']); 
end
diary([nDPath 'log' filesep 'combine_sum.log']);
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

%% Check to see if overwrite has been defined
if ~exist('overwrite') 
	overwrite = false;
end

%% Control for possible overwrite
mDecoder.overwrite = overwrite;

%% Run the raw sum combination
mDecoder.CombineRawSum();

%% Create the summation report
mDecoder.GenerateSummationReport();

%% Archive analysis
PageBreak();
display(['...completed in ' num2str(toc(scriptTimer)) ' s']);
display(['Completed at ' datestr(now)]);

% Turn off diary
diary off;