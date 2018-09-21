function [dataOrg, metaData, parameters] = LoadDataOrganization(dataOrgPath, varargin)
% dataOrg = LoadDataOrganization(dataOrgPath, varargin)
% This function opens a data organization file
%
%--------------------------------------------------------------------------
% Necessary Inputs: 
%   dataOrgPath -- A valid path to a data organization file
%--------------------------------------------------------------------------
% Outputs: 
%   --None
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
% 
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 21, 2017
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2018.
%-----------------------------------------------------------------------------------------------------------------------------------
% This function loads a data organization csv file

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);

% Parameters for displaying progress
defaults(end+1,:) = {'verbose', 'boolean', false};      % Display progress?

% Parameters for handling internal delimiters in fields
defaults(end+1,:) = {'internalDelimiter', 'char', ';'}; % The internal delimiter for individual fields

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1 || ~exist(dataOrgPath)
    error('matlabFunctions:invalidArguments', ...
        'A valid path to a data organization file must be provided');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Identify file type
% -------------------------------------------------------------------------
[~, ~, fileType] = fileparts(dataOrgPath);

% -------------------------------------------------------------------------
% Load data organization file
% -------------------------------------------------------------------------
switch fileType
    case '.csv'
        % Load csv file using the table2struct approach
        dataOrg = table2struct(readtable(dataOrgPath, 'Delimiter', ','));
        
    otherwise
        error('matlabFunctions:invalidFileType', ...
            'Only csv files are supported.');
end

% -------------------------------------------------------------------------
% Handle internal delimiters
% -------------------------------------------------------------------------
possibleFields = {'frame', 'zPos'}; % Fields that could have internal delimiters
for f=1:length(possibleFields)
    for d=1:length(dataOrg)
        if isfield(dataOrg(d), possibleFields{f}) % Confirm that the field exists
            if ischar(dataOrg(d).(possibleFields{f})) % Check to see if it is a char--a flag that it could not be converted to numbers
                dataOrg(d).(possibleFields{f}) = cell2mat(cellfun(@str2num, ...
                    strsplit(dataOrg(d).(possibleFields{f}), ';'), 'UniformOutput', false));
            end
        end
    end
end

% -------------------------------------------------------------------------
% Build metadata structure
% -------------------------------------------------------------------------
% The number of data channels
metaData.numDataChannels = length(dataOrg);

% The number and location of z-stacks
if isfield(dataOrg, 'zPos')
    uniqueZPos = unique([dataOrg.zPos]);
else
    uniqueZPos = 0;
	% Add zPos to each dataOrg entry
	for d=1:length(dataOrg)
		dataOrg(d).zPos = 0;
	end
end
metaData.zPos = sort(uniqueZPos, 'ascend'); % Order the z from small to large
metaData.numZPos = length(uniqueZPos);

% -------------------------------------------------------------------------
% Display data organization file if requested
% -------------------------------------------------------------------------
if parameters.verbose
    PageBreak();
    display(['Loaded data organization file: ' dataOrgPath]);
    display(['Found ' num2str(metaData.numDataChannels) ' data channels']);
    display(['Found ' num2str(metaData.numZPos) ' z-stacks']);
end




