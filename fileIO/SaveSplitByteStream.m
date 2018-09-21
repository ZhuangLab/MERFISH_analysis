function [saveBase, parameters] = SaveSplitByteStream(filePath, variable, blockSize, varargin)
% ------------------------------------------------------------------------
% savePath = SaveSplitByteStream(filePath, variable, blockSize)
% This function saves a bytestream version of a matlab object as a '.matb' 
%   file. Unlinke SaveAsByteStream, this function splits arrays into blockSize, 
%   so as to overcome the issues associated with matlab's undocumented 
%   byte stream conversion. The full array can be loaded with
%   LoadSplitByteStream().
%--------------------------------------------------------------------------
% Necessary Inputs
% filePath/A string to a the path at which the variable will be saved.
% variable/A matlab array of anything 
% blockSize/The number of array alements to save in each block
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
if any(nargin < [1 2 3])
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
defaults(end+1,:) = {'splitPostFix', 'char', '_'};

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

if isempty(ext)
    ext = '.matb';
end

%--------------------------------------------------------------------------
% Build block indices
%--------------------------------------------------------------------------
indices = [0:blockSize:(numel(variable)-1) numel(variable)];

if parameters.verbose
    disp(['Splitting file into ' num2str(length(indices)) ' chunks']);
end

%--------------------------------------------------------------------------
% Construct reload information
%--------------------------------------------------------------------------
sbsInfo.size = size(variable);                  % How to reshape the flatten array
sbsInfo.numBlocks = length(indices)-1;            % How many blocks were saved
sbsInfo.splitPostFix = parameters.splitPostFix; % A postfix name for the split files
sbsInfo.version = '1.1';                        % Version information to handle future upgrades
sbsInfo.isByteStream = 'Y';                     % A flag that is likely redundant but may help handle corrupted/misidentified files
sbsInfo.fileNames = {};                         % A list of the file names to assist in rapidly loading and reconstruction
sbsInfo.uuID = char(java.util.UUID.randomUUID); % A unique ID that allows rapid comparison of copies of a byte stream

%--------------------------------------------------------------------------
% Flatten variable
%--------------------------------------------------------------------------
variable = reshape(variable, [numel(variable) 1]);

%--------------------------------------------------------------------------
% Save individual blocks
%--------------------------------------------------------------------------
saveBase = [fileName parameters.splitPostFix];
for i=1:(length(indices)-1)
    % Generate name
    localName = [saveBase num2str(i) ext];
    sbsInfo.fileNames{end+1} = localName;
    
    % Save byte stream
    SaveAsByteStream([fileDir filesep localName], variable((indices(i)+1):indices(i+1)), 'verbose', parameters.verbose);
end

%--------------------------------------------------------------------------
% Save reload information
%--------------------------------------------------------------------------
SaveAsByteStream([fileDir filesep fileName ext], sbsInfo);

