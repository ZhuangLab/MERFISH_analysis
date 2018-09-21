function parameters = CatBinaryFiles(filePath1, filePath2, varargin)
% ------------------------------------------------------------------------
% CatBinaryFiles(filePath1, filePath2, varargin)
% This function combines the contents of two binary files. It reads the
% binary file from file path 2 and appends it to file path 1. 
%
% These binary files must have matching headers to be combined. 
%--------------------------------------------------------------------------
% Necessary Inputs: 
%   filePath1 -- A valid path to binary file to which another binary file
%       will be appended. If this file does not exist, the original file
%       will simply be copied. 
%   filePath2 -- A valid path to binary file that will be appended to the
%       original file
%--------------------------------------------------------------------------
% Outputs: 
%   --None
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
% --verbose (Boolean): Display progress/output?
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 21, 2017
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2018.
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);

% Parameters for displaying progress
defaults(end+1,:) = {'verbose', 'boolean', false};  % Display progress?

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 2 || ~exist(filePath2)
    error('matlabFunctions:invalidArguments', ...
        'Two paths must be provided and the path to the file to add must be valid');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Check to see if the file to write to exists
% -------------------------------------------------------------------------
if ~exist(filePath1)
    copyfile(filePath2, filePath1, 'f');
    if parameters.verbose
        display(['File 1 does not exist. So File 2 was copied to File 1']);
    end
    return;
end

% -------------------------------------------------------------------------
% Confirm that the files have the same headers
% -------------------------------------------------------------------------
[fileHeader1, fileLayout1] = ReadBinaryFileHeader(filePath1);
[fileHeader2, fileLayout2] = ReadBinaryFileHeader(filePath2);

% Abort if either are corrupt
if fileHeader1.isCorrupt || fileHeader2.isCorrupt
    error('matlabFunctions:corruptFiles', 'One of the binary files is corrupt');
end

% Check that both are the same version
if ~(fileHeader1.version == fileHeader2.version)
    error('matlabFunctions:incompatibleFiles', 'The two binary files do not have the same version');
end

% Check layout (for a specific version)
switch fileHeader1.version
    case 1
        % Check header properties
        if ~(fileHeader1.headerLength == fileHeader2.headerLength) || ... % Length
                ~(all(size(fileLayout1) == size(fileLayout2))) || ... % Number of entries
                ~all(cellfun(@(x,y)strcmp(x,y),fileLayout1(:,1), fileLayout2(:,1))) || ... % Variable type
                ~all(cellfun(@(x,y)all(x==y),fileLayout1(:,2), fileLayout2(:,2))) || ... % Size
                ~all(cellfun(@(x,y)strcmp(x,y),fileLayout1(:,3), fileLayout2(:,3)))

            error('matlabFunctions:incompatibleFiles', 'The two binary files do not have the same layout');
        end
    otherwise
        error('matlabFunctions::binaryFileVersion', 'This version is not yet supported.');
end


% -------------------------------------------------------------------------
% Open the second file, read, and concatenate to the first
% -------------------------------------------------------------------------
% Mark the first file as open
fid1 = fopen(filePath1, 'r+');
fseek(fid1, 1, 'bof');
fwrite(fid1, 1, 'uint8');
fwrite(fid1, 0, 'uint32');
fclose(fid1);

% Open the second file to read in the data
fid2 = fopen(filePath2, 'r');

% Advance the file pointer for the second file past the header
fseek(fid2, fileHeader2.headerLength + 10, 'bof'); % 10 = 1 + 1 + 4 + 4 
                                                   % (version, is corrupt, numEntries, headerLength)  
data = fread(fid2, inf, 'uint8'); % Read the entire file 
fclose(fid2); % Close the first file

% Append the data to the first file
fid1 = fopen(filePath1, 'A'); % Reopen in buffered append mode
fwrite(fid1, data, 'uint8');
fclose(fid1); % Close the first file

% Open the first file to mark as close and update the number of entries
fid1 = fopen(filePath1, 'r+');
fseek(fid1, 1, 'bof');
fwrite(fid1, 0, 'uint8');
fwrite(fid1, fileHeader1.numEntries + fileHeader2.numEntries, 'uint32');
fclose(fid1);
