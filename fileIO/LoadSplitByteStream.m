function [data, sbsInfo, parameters] = LoadSplitByteStream(filePath, varargin)
% ------------------------------------------------------------------------
% data = LoadSplitByteStream(filePath)
% This function loads a bytestream version of a matlab object as a '.matb' 
%   file and reconstitutes it into a variable.  In particular, this function 
%   is designed to handle split bytestreams as produced by SaveSplitByteStream().
%
%--------------------------------------------------------------------------
% Necessary Inputs
% filePath/A string to a bytestream file
%--------------------------------------------------------------------------
% Outputs
% variable/The reconstituted byteStream variable. 
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 21, 2017
%--------------------------------------------------------------------------
% Based on
% http://undocumentedmatlab.com/blog/serializing-deserializing-matlab-data
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Parse necessary inputs
%--------------------------------------------------------------------------
if any(nargin < 1)
    error('matlabFunctions:invalidArguments', 'Both a file path and an object are required');
end

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true};

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

%--------------------------------------------------------------------------
% Print updates
%--------------------------------------------------------------------------
if parameters.verbose
    tic;
    display(['Loading split bytestream info from ' filePath]);
end

%--------------------------------------------------------------------------
% Load split byte stream info
%--------------------------------------------------------------------------
sbsInfo = LoadSplitByteStreamHeader(filePath);

%--------------------------------------------------------------------------
% Find base path
%--------------------------------------------------------------------------
[basePath] = fileparts(filePath);

%--------------------------------------------------------------------------
% Load based on version
%--------------------------------------------------------------------------
switch sbsInfo.version
    case {'1.0', '1.1'}
        data = [];
        for i=1:sbsInfo.numBlocks
            data = cat(1, data, LoadByteStream([basePath filesep sbsInfo.fileNames{i}], 'verbose', parameters.verbose));
        end
        
    otherwise
        error('matlabFunctions:invalidFile', 'The specified split bytestream version is not valid');
end

%--------------------------------------------------------------------------
% Reshape output
%--------------------------------------------------------------------------
data = reshape(data, sbsInfo.size);
