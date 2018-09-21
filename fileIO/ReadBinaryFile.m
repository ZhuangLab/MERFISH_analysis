function [structArray, flatHeader, parameters] = ReadBinaryFile(filePath, varargin)
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
%   filePath -- A valid path to the a binary file
%--------------------------------------------------------------------------
% Outputs: 
%   --None
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
defaults(end+1,:) = {'verbose', 'boolean', false};      % Display progress?
defaults(end+1,:) = {'format', ...
    {'array', 'structure', 'structureArray'}, 'array'}; % The format of the returned data
defaults(end+1,:) = {'first', 'positive', []};          % The first index to load
defaults(end+1,:) = {'last', 'positive', []};          % The last index to load
defaults(end+1,:) = {'fieldsToLoad', 'cell', {}};       % The fields to load
% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1 || ~exist(filePath)
    error('matlabFunctions:invalidArguments', 'A valid path must be provided.');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Read header
% -------------------------------------------------------------------------
[flatHeader, fileLayout] = ReadBinaryFileHeader(filePath);
% Confirm that the file is not corrupt
if flatHeader.isCorrupt
    error('matlabFunctions:invalidFile', 'File appears to be corrupt');
end

% -------------------------------------------------------------------------
% Coerce and check load indices
% -------------------------------------------------------------------------
if isempty(parameters.first)
    parameters.first = 1;
end
if isempty(parameters.last)
    parameters.last = flatHeader.numEntries;
end
if parameters.first < 1
    error('matlabFunctions:invalidArguments', 'The first load index cannot be smllaer than 1');
end
if parameters.last > flatHeader.numEntries
   error('matlabFunctions:invalidArguments', 'The last load index cannot be larger than the number of entries in the file');
end  

% -------------------------------------------------------------------------
% Check requested fields
% -------------------------------------------------------------------------
if ~isempty(parameters.fieldsToLoad)
    foundFields = fileLayout(:,3);
    if ~isempty(setdiff(parameters.fieldsToLoad, foundFields))
        error('matlabFunctions::invalidArguments', 'Some of the requested fields are not present in the binary file');
    end
end

% -------------------------------------------------------------------------
% Determine the data block size
% -------------------------------------------------------------------------
% Create useful lookup table
dataSize = cat(1,{'uint8', 1}, ...
    {'uint16', 2}, ...
    {'uint32', 4}, ...
    {'uint64', 8}, ...
    {'int8', 1}, ...
    {'int16', 2}, ...
    {'int32', 4}, ...
    {'single', 4}, ...
    {'double', 8});

% Determine data types
[validDataType, typeInd] = ismember(fileLayout(:,1), dataSize(:,1));

% Check for valid data types
if ~all(validDataType)
    error('matlabFunctions:invalidDataType', 'Found an invalid data type');
end
%Determine the number of bytes per entry data type
numBytes = cat(1,dataSize{typeInd, 2});

% Determine the base range for each data set
blockSize = numBytes.* cellfun(@prod, fileLayout(:,2));
blockOffset = cumsum([0 blockSize']);
totalBlockSize = blockOffset(end); % In bytes

% -------------------------------------------------------------------------
% Display properties of file if requested
% -------------------------------------------------------------------------
if parameters.verbose
    PageBreak();
    display(['Loading: ' filePath]);
    display(['    Writer version: ' num2str(flatHeader.version)]);
    display(['    Number of entries: ' num2str(flatHeader.numEntries)]);
    display(['    Loading: ' num2str(parameters.first) ' to ' num2str(parameters.last)]);
    display(['    Found fields: ']);
    for i=1:size(fileLayout,1)
        display(['        ' fileLayout{i,3} ': ' fileLayout{i,1} ...
            ' of size ' num2str(fileLayout{i,2})]);
    end
    loadTimer = tic;
end 
% -------------------------------------------------------------------------
% Switch reader based on version
% -------------------------------------------------------------------------
switch flatHeader.version
    case 1
        switch parameters.format
            case 'array'; % Simple, but slow and memory inefficient
                if isempty(parameters.fieldsToLoad)
                    % -------------------------------------------------------------------------
                    % Create memory map 
                    % -------------------------------------------------------------------------
                    memoryMap = memmapfile(filePath, ...
                        'Format', fileLayout, ...
                        'Writable', false, ...
                        'offset', 10 + flatHeader.headerLength + totalBlockSize*(parameters.first-1), ...
                        'Repeat', (parameters.last - parameters.first + 1));

                    % Return contents of memory map
                    structArray = memoryMap.Data;
                else
                    warning('Under construction... very slow as currently implemented....');
                    % Open file for reading
                    fid = fopen(filePath, 'r');
                                        
                    % Load data fields for structure
                    varArgIn = cell(1, 2*length(parameters.fieldsToLoad));
                    varArgIn(1:2:end) = parameters.fieldsToLoad; % Define field names
                    % Load data for each field
                    for f=1:length(parameters.fieldsToLoad)
                        localField = parameters.fieldsToLoad{f};
                        fieldID = find(strcmp(fileLayout(:,3), localField));
                        localData = zeros([1 prod(fileLayout{fieldID,2})*(parameters.last - parameters.first + 1)], ...
                            fileLayout{fieldID,1});
                        for d=parameters.first:parameters.last
                            fileOffset = 10 + flatHeader.headerLength + totalBlockSize*(d-1) + blockOffset(fieldID);
                            fseek(fid, fileOffset, 'bof');
                            localData(d) = fread(fid, fileLayout{fieldID,2}, fileLayout{fieldID,1}, totalBlockSize);
                        end
                        varArgIn{2*f} = localData;
                    end

                    % Close file
                    fclose(fid);
                    
                    % Create structure
                    structArray = struct(varArgIn{:}); 
                end

            case 'structureArray'
                error('Under construction!');

                % -------------------------------------------------------------------------
                % Create memory map 
                % -------------------------------------------------------------------------
                memoryMap = memmapfile(filePath, ...
                    'Format', fileLayout, ...
                    'Writable', false, ...
                    'offset', 10 + flatHeader.headerLength, ...
                    'Repeat', flatHeader.numEntries);

                % Return contents of memory map
                structArray = memoryMap.Data;
                
                
                % Prepare data for creating structure array
                argIn = cell(1, 2*size(fileLayout,1));
                argIn(1:2:end) = fileLayout(:,3); % Set field names
                
                for b=2:2:2*size(fileLayout,1)
                    argIn{b} = cat(1,structArray.(argIn{b-1}));
                end
                
                structArray = StructureArray(argIn{:});
                
                %%%% THE FOLLOWING CODE IS A WORK IN PROGRESS: AND PROMISES TO BE VERY FAST
%                 
%                 % -------------------------------------------------------------------------
%                 % Determine basic mapping for loading file quickly
%                 % -------------------------------------------------------------------------
%                 [validDataType, typeInd] = ismember(fileLayout(:,1), dataSize(:,1));
% 
%                 % Check for valid data types
%                 if ~all(validDataType)
%                     error('matlabFunctions:invalidDataType', 'Found an invalid data type');
%                 end
% 
%                 % Determine the number of bytes per dataType
%                 numBytes = cat(1,dataSize{typeInd, 2});
% 
%                 % Determine the base range for each data set
%                 blockSize = numBytes.* cellfun(@prod, fileLayout(:,2));
%                 blockOffset = cumsum([0 blockSize']);
%                 totalBlockSize = blockOffset(end);
%                 
%                 % -------------------------------------------------------------------------
%                 % Load data as single byte units, parse, and typecast
%                 % -------------------------------------------------------------------------
%                 % Create memory map
%                 memoryMap = memmapfile(filePath, ...
%                     'Format', 'uint8', ...
%                     'Writable', false, ...
%                     'offset', 10 + flatHeader.headerLength, ...
%                     'Repeat', inf);
%        
%                 % Extract data
%                 byteData = memoryMap.Data;
%                 
%                 % Prepare data for creating structure array
%                 argIn = cell(1, 2*size(fileLayout,1));
%                 argIn(1:2:end) = fileLayout(:,3); % Set field names
%                 
%                 % Faster indexing
%                 boolInd = false(1, length(byteData));
%                 
%                 % Extract, typecast, and reshape data
%                 for b=1:size(fileLayout(:,3))
%                     
%                     blockInd = repmat([1:blockSize(b)]', [1 flatHeader.numEntries]) + ...
%                         repmat(totalBlockSize*([1:flatHeader.numEntries] - 1), [blockSize(b) 1]);
%                     blockInd = blockInd(:);
%                     boolInd(blockInd) = true;
%                     
%                     argIn{2*b} = squeeze( reshape(typecast( byteData(boolInd), fileLayout{b,1}), ...
%                         [fileLayout{b,2} flatHeader.numEntries]) );
%                         
%                 end
%                 
%                 % Create structure array
%                 structArray = StructureArray(argIn{:});
                
            case 'structure'
                error('This format is not yet supported');
        end
                    
    otherwise
        error('matlabFunctions:invalidFile', 'The version is not supported.');
end
% -------------------------------------------------------------------------
% Load data as single byte units, parse, and typecast
% -------------------------------------------------------------------------
if parameters.verbose
    display(['...completed in ' num2str(toc(loadTimer)) ' s']);
end


