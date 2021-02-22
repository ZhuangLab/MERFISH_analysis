function [savePath, parameters] = SaveAsByteStream(filePath, variable, varargin)
% ------------------------------------------------------------------------
% savePath = SaveAsByteStream(filePath, variable)
% This function saves a bytestream version of a matlab object as a '.matb' 
%   file.  It can be loaded and 'reconstituted' by the function
%   LoadByteStream().
%--------------------------------------------------------------------------
% Necessary Inputs
% filePath/A string to a the path at which the variable will be saved.
% variable/Any matlab object, e.g. class, struct, array, etc. 
%--------------------------------------------------------------------------
% Outputs
% savePath/The path to the saved file.  Empty if saving was unsuccesful. 
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
if any(nargin < [1 2])
    error('matlabFunctions:invalidArguments', 'Both a file path and an object are required');
end

if ~isstr(filePath)
    error('matlabFunctions:invalidArguments', 'The first argument must be a file path');
end

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true};
defaults(end+1,:) = {'overwrite', 'boolean', true};

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

%--------------------------------------------------------------------------
% Check extension
%--------------------------------------------------------------------------
[fileDir, fileName, ext] = fileparts(filePath);

if ~strcmp(ext, '.matb')
    warning('matlabFunctions:invalidArguments', 'It is strongly recommended to save bytestreams as .matb files');
end

%--------------------------------------------------------------------------
% Overwrite protection
%--------------------------------------------------------------------------
if ~parameters.overwrite
    while exist(filePath,'file') ~= 0
        filePath = IncrementSaveName(filePath);
        if parameters.verbose
            warning('detected existing file.  Increased counter on save file.'); 
        end
    end
end

%--------------------------------------------------------------------------
% Open file
%--------------------------------------------------------------------------
fid = fopen(filePath, 'w');

if fid < 0
    savePath = [];
    error('matlabFunctions:invalidArguments', 'Invalid file path');
end

%--------------------------------------------------------------------------
% Print updates
%--------------------------------------------------------------------------
if parameters.verbose
    tic;
    display(['Saving ' filePath]);
end

%--------------------------------------------------------------------------
% Create byte stream
%--------------------------------------------------------------------------
try
    byteStream = getByteStreamFromArray(variable);
catch
    fclose(fid);
    delete(filePath);
    assignin('base', 'bigVariable', variable);
    error('matlabFunctions:byteStreamError', 'Unable to create bytestream. It is likely that the file is too large');
end

%--------------------------------------------------------------------------
% Save byte stream
%--------------------------------------------------------------------------
fwrite(fid, byteStream, 'uint8');

%--------------------------------------------------------------------------
% Close file id
%--------------------------------------------------------------------------
fclose(fid);
savePath = filePath;

%--------------------------------------------------------------------------
% Print updates
%--------------------------------------------------------------------------
if parameters.verbose
    display(['.... finished in ' num2str(toc)]);
end

