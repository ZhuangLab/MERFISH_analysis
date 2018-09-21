function [flatHeader, fileLayout] = ReadBinaryFileHeader(filePath)
% ------------------------------------------------------------------------
% [flatHeader, fileLayout] = ReadBinaryFileHeader(filePath) 
% This function reads the header of the generic binary file format defined
% in WriteBinaryFile.
%
% This header contains two basic structures: 1) A fixed flat header which
% contains information such as the number of entries, the binary file
% format version, etc (returned as the structure flatHeader) and 2) a
% extensible, variable header which describes the name, type, and shape of
% all data elements stored in the binary file (returned as fileLayout).
%
% File layout is returned in a format that can be interpreted by
% memmapfile. See this documentation for the organization of this cell
% array.
%

%--------------------------------------------------------------------------
% Necessary Inputs: 
%   filePath -- A valid path to the a binary file
%--------------------------------------------------------------------------
% Outputs: 
%
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
% 
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 21, 2017
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2018.

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1 || ~exist(filePath)
    error('matlabFunctions:invalidArguments', 'A valid path must be provided.');
end

% -------------------------------------------------------------------------
% Open file
% -------------------------------------------------------------------------
fid = fopen(filePath, 'r');

if fid < 0
    error('matlabFunctions:invalidArguments', 'Could not open file');
end

% -------------------------------------------------------------------------
% Read version
% -------------------------------------------------------------------------
flatHeader.version = fread(fid, 1, 'uint8'); % Read version

% -------------------------------------------------------------------------
% Switch reader based on version
% -------------------------------------------------------------------------
switch flatHeader.version
    case 1
        % -------------------------------------------------------------------------
        % Reader fixed header
        % -------------------------------------------------------------------------
        flatHeader.isCorrupt = fread(fid, 1, 'uint8');
        flatHeader.numEntries = fread(fid, 1, 'uint32');
        flatHeader.headerLength = fread(fid, 1, 'uint32');
        
        % -------------------------------------------------------------------------
        % Read variable header, i.e. file format 
        % -------------------------------------------------------------------------
        layoutString = fread(fid, flatHeader.headerLength, 'uchar=>char');
        layoutParts = strsplit(layoutString', ',');
        fileLayout = reshape(layoutParts, [3 length(layoutParts)/3])';
        fileLayout(:,2) = cellfun(@str2num, fileLayout(:,2) , 'UniformOutput', false);
        fileLayout(:, [1 3]) = fileLayout(:, [3 1]);
        % -------------------------------------------------------------------------
        % Close file
        % -------------------------------------------------------------------------
        fclose(fid);
                
    otherwise
        error('matlabFunctions:invalidFile', 'The version is not supported.');
end
