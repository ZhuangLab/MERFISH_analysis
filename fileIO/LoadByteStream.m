function [data, parameters] = LoadByteStream(filePath, varargin)
% ------------------------------------------------------------------------
% data = LoadByteStream(filePath)
% This function loads a bytestream version of a matlab object as a '.matb' 
%   file and reconstitutes it into a variable.  These objects can be saved
%   with SaveAsByteStream().
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
    display(['Loading ' filePath]);
end


    %--------------------------------------------------------------------------
    % Open file
    %--------------------------------------------------------------------------
    fid = fopen(filePath, 'r');

    if fid < 0
        error('matlabFunctions:invalidArguments', 'Invalid file path');
    end

    %--------------------------------------------------------------------------
    % Read byteStream
    %--------------------------------------------------------------------------

    byteStream = fread(fid, inf, '*uint8');

%--------------------------------------------------------------------------
% Convert to matlab object
%--------------------------------------------------------------------------
data = getArrayFromByteStream(byteStream);

%--------------------------------------------------------------------------
% Close file id
%--------------------------------------------------------------------------
fclose(fid);

%--------------------------------------------------------------------------
% Print updates
%--------------------------------------------------------------------------
if parameters.verbose
    display(['.... finished in ' num2str(toc)]);
end