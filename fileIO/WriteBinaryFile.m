function parameters = WriteBinaryFile(filePath, structArray, varargin)
% ------------------------------------------------------------------------
% WriteBinaryFile(filePath, struct, varargin) 
% This function writes a structure array to a custom binary file format
% defined by the fields in the structure array. 
%
% Structure fields may be arrays of any size and format but every instance
% of this field in each structure in the array MUST be the same size and
% type. 
%--------------------------------------------------------------------------
% Necessary Inputs: 
%   filePath -- A valid path to the file that will be saved
%   structArray -- The array to be saved. 
%--------------------------------------------------------------------------
% Outputs: 
%   --None
%--------------------------------------------------------------------------
% File organization
% Version number -- A uint8 that specifies the reader/writer version number
% Corrupt -- A uint8 that is set to 1 when the file is opened and 0 to when
%    closed. If the file is improperly closed, a 1 will indicate that it is
%    corrupt
% Number of entries -- A uint32 that specifies the number of entries
% Header length -- A uint32 that specifies the length of the following
% header
% Header -- A character array that can be parsed to
%    determine the layout of the file. Each entry is written as follows
%    field name ,  field dimensions ,  field class ,
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
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);

% Parameters for displaying progress
defaults(end+1,:) = {'verbose', 'boolean', false};  % Display progress?
defaults(end+1,:) = {'overwrite', 'boolean', true}; % Overwrite file? If not, append data
defaults(end+1,:) = {'append', 'boolean', false};    % Append to file, but only if overwrite is not true

% -------------------------------------------------------------------------
% Version number
% -------------------------------------------------------------------------
version = 1;

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 2 || ~isstruct(structArray)
    error('matlabFunctions:invalidArguments', 'A path and a structure array must be provided.');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Parse structure
% -------------------------------------------------------------------------
% Determine fields and number of entries
fieldsToWrite = fields(structArray);
numEntries = length(structArray);

% Determine class types
dataTypes = {};
dataSize = {};
for f=1:length(fieldsToWrite)
    dataTypes{f} = class(structArray(1).(fieldsToWrite{f}));
    dataSize{f} = num2str(size(structArray(1).(fieldsToWrite{f})));
end

% Check validity of dataTypes
if ~all(ismember(dataTypes, {'int8', 'int16', 'int32', 'int64', ...
        'uint8', 'uint16', 'uint32', 'uint64', 'single', 'double'}))
    error('matlabFunctions:invalidArguments', 'One field has an unsupported data type');
end

% Prepare header string
headerString = {};
deliminter = repmat({','}, [1 length(fieldsToWrite)]);
headerString(1:6:6*length(fieldsToWrite)) = fieldsToWrite;
headerString(2:6:6*length(fieldsToWrite)) = deliminter;
headerString(3:6:6*length(fieldsToWrite)) = dataSize;
headerString(4:6:6*length(fieldsToWrite)) = deliminter;
headerString(5:6:6*length(fieldsToWrite)) = dataTypes;
headerString(6:6:6*length(fieldsToWrite)) = deliminter;
headerString = cat(2, headerString{:});
headerString = headerString(1:(end-1)); % Remove final delimiter
headerLength = length(headerString);

% -------------------------------------------------------------------------
% Check to see if file exists
% -------------------------------------------------------------------------
isAppend = false;
if exist(filePath)
    if parameters.overwrite && ~parameters.append % If overwriting, delete file
        delete(filePath);
        if parameters.verbose
            display(['Deleting existing file: ' filePath]);
        end
    elseif ~parameters.append % If not overwriting, and not appending generate error
        error('matlabFunctions:existingFile', 'Found existing file.');
    else % The file exists and the user has explicitly requested to append
        % Read previous file header
        [flatHeader, fileLayout] = ReadBinaryFileHeader(filePath);
        
        % Confirm that the file is not corrupt
        if flatHeader.isCorrupt
            error('matlabFunctions:invalidFile', 'File appears to be corrupt');
        end
        
        % Confirm that the file layout is identical
        if flatHeader.headerLength ~= headerLength | ...
                ~all(strcmp(fieldsToWrite, fileLayout(:,3))) | ...
                ~all(cellfun(@(x,y)all(str2num(x)==y), dataSize', fileLayout(:,2))) | ...
                ~all(strcmp(dataTypes', fileLayout(:,1)))
            error('matlabFunctions:incompatibleHeaders', 'Cannot append binary data with a different organization');
        end
        
        % Set the isAppend flag
        isAppend = true;
    end
    % Run checks
end

% -------------------------------------------------------------------------
% Handle the appending/writing scenarios
% -------------------------------------------------------------------------
if ~isAppend
    % -------------------------------------------------------------------------
    % Open file
    % -------------------------------------------------------------------------
    fid = fopen(filePath, 'W');

    if fid < 0
        error('matlabFunctions:invalidArguments', 'Could not open file');
    end

    % -------------------------------------------------------------------------
    % Write header
    % -------------------------------------------------------------------------
    % Fixed header
    fwrite(fid, version, 'uint8');      % version
    fwrite(fid, 1, 'uint8');            % Mark file as open until closed properly
    fwrite(fid, 0, 'uint32');           % number of entries
    fwrite(fid, headerLength, 'uint32');% length of following header 

    % Variable header (file layout)
    fwrite(fid, headerString, 'uchar');
    
else % Appending
    % -------------------------------------------------------------------------
    % Open file to mark the corruption flag
    % -------------------------------------------------------------------------
    fid = fopen(filePath, 'r+');
    if fid < 0
        error('matlabFunctions:invalidArguments', 'Could not open file');
    end

    fseek(fid, 1, 'bof');
    fwrite(fid, 0, 'uint8'); % Mark file as open or corrupted.
    
    fclose(fid);
    
    % -------------------------------------------------------------------------
    % Reopen the file in append mode with automatic buffer flushing
    % -------------------------------------------------------------------------
    fid = fopen(filePath, 'A');
    if fid < 0
        error('matlabFunctions:invalidArguments', 'Could not open file');
    end
end

% -------------------------------------------------------------------------
% Write data one entry at a time
% -------------------------------------------------------------------------
for i=1:numEntries
    for f=1:length(fieldsToWrite)
        fwrite(fid, structArray(i).(fieldsToWrite{f}), dataTypes{f});
    end
end

% -------------------------------------------------------------------------
% Update the number of entries if appending
% -------------------------------------------------------------------------
if isAppend
    numEntries = numEntries + flatHeader.numEntries;
end

% -------------------------------------------------------------------------
% Finalize file by marking the number of entries
% -------------------------------------------------------------------------
% Close and reopen to force buffer/flush write
fclose(fid);
fid = fopen(filePath, 'r+'); % Must open in read/write option to modify the header
if fid < 0
    error('matlabFunctions:invalidArguments', 'Could not reopen the file to finalize the header.');
end
fseek(fid, 1, 'bof');              % Skip version
fwrite(fid, 0, 'uint8');           % Mark file as closed
fwrite(fid, numEntries, 'uint32'); % Update number of entries

% -------------------------------------------------------------------------
% Close file
% -------------------------------------------------------------------------
fclose(fid);
    
    
