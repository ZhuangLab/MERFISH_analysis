function [parsedFileStruct, parameters] = BuildFileStructure(folderPath, varargin)
% ------------------------------------------------------------------------
% [parsedFileStruct] = BuildFileStructure(folderPath, varargin)
% This function parses the names of all files within a dataPath (and included
% directories) to produce a structure array of properties specified by
% those names. 
%--------------------------------------------------------------------------
% Necessary Inputs
% folderPath/ A path to the desired folder
%--------------------------------------------------------------------------
% Outputs
% parsedFileStruct/ A structure array containing a variety of default
%   and user specified fields
%   -- name: The local name of the file
%   -- filePath: The full path to the file
%   -- Additional fields can be specified
%--------------------------------------------------------------------------
% Variable Inputs
% 'fileExt'/string ('*'): The extension of files to be returned. By default all
%   files are returned.
% 'delimiters'/cell of strings ({'_'}): The characters used to split a
%   complex file name into parts. The default is '_'.
% 'fieldNames'/cell of strings ({'field1'}): The names assigned to portion of a
%   split string. Any entry without a field name will not be assigned.  
% 'fieldConv'/cell of conversion functions ({'char'}): The function used to
%   convert the parsed entry to a data type. If not specified for a field,
%   it will remain a string.
% 'appendExtraFields'/boolean ('False'): If there are extra fields beyond
%   the specified fields, they will be combined to generate the final field. 
%--------------------------------------------------------------------------
% Example:
% Consider a name STORM_01_03.dax, where the delimiter is '_' the first
% entry represents the name
%-------------------------------------------------------------------------- 
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 5, 2014
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'fileExt', 'string', '*'}; % File extension to return
defaults(end+1,:) = {'delimiters', 'cell', {'_'}}; % Delimiters to use to split
defaults(end+1,:) = {'fieldNames', 'cell', {'field1'}}; % FieldNames
defaults(end+1,:) = {'fieldConv', 'cell', {@char}}; % Conversion functions
defaults(end+1,:) = {'appendExtraFields', 'boolean', false}; 
defaults(end+1,:) = {'excludeFlags', 'cell', {}}; % A cell array of strings that are excluded
defaults(end+1,:) = {'requireFlag', 'string', ''}; % A string that is required
defaults(end+1,:) = {'requireExactMatch','boolean',false};
defaults(end+1,:) = {'containsDelimiters','positive', []}; % An integer specifying a field that might have internal delimiters
defaults(end+1,:) = {'regExp', 'string', []}; % A regular expression with tokens that match fieldNames

% -------------------------------------------------------------------------
% Parse and normalize necessary input
% -------------------------------------------------------------------------
if nargin < 1 || ~(exist(folderPath) == 7) % 7=folder
    error('matlabFunctions:invalidArguments', 'A valid folder is required.');
end
if folderPath(end) ~= filesep
    folderPath(end+1) = filesep;
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Check for field names that will overwrite hard coded field names
% -------------------------------------------------------------------------
if any(ismember({'name', 'filePath'}, parameters.fieldNames))
    error('matlabFunctions:invalidArguments', 'name and filePath are protected names.');
end

% -------------------------------------------------------------------------
% Define old functionality of internal delimiters
% -------------------------------------------------------------------------
if isempty(parameters.containsDelimiters)
    parameters.containsDelimiters =length(parameters.fieldConv);
end

% -------------------------------------------------------------------------
% Check contains delimiters
% -------------------------------------------------------------------------
if length(parameters.containsDelimiters) > 1
    error('matlabFunctions:invalidDelimiters', 'Only one field can have internal delimiters.');
end

% -------------------------------------------------------------------------
% Additional parsing of parameters
% -------------------------------------------------------------------------
% Add extension delimiter
parameters.delimiters{end+1} = '.';
% Coerce unfilled conversion functions
for i=(length(parameters.fieldConv)+1):length(parameters.fieldNames)
    parameters.fieldConv{i} = @char;
end

% -------------------------------------------------------------------------
% Find files
% -------------------------------------------------------------------------
if ~isempty(parameters.requireFlag)
    fileData = dir([folderPath '*' parameters.requireFlag '*.' parameters.fileExt]);
else
    fileData = dir([folderPath '*.' parameters.fileExt]);
end

% -------------------------------------------------------------------------
% Exclude files with any of the excluded strings
% -------------------------------------------------------------------------
if ~isempty(parameters.excludeFlags);
    excludeIdx =false(1,length(fileData));
    for e=1:length(parameters.excludeFlags)
         excludeIdx  = ~cellfun(@isempty, strfind({fileData.name},parameters.excludeFlags{e})) | excludeIdx;
    end
    fileData(excludeIdx) = [];
end

% -------------------------------------------------------------------------
% Parse Names
% -------------------------------------------------------------------------
parsedFileStruct = [];
count = 1;
for i=1:length(fileData)
    % Switch on whether or not a regular expression was provided
    if ~isempty(parameters.regExp)
        % Apply regular expression and capture token values
        localStruct = regexp(fileData(i).name, parameters.regExp, 'names');
        
        % Handle case in which regexp did not match
        if isempty(localStruct)
            continue;
        end
        
        % Apply field conversions
        for f=1:length(parameters.fieldNames)
            if ~isfield(localStruct, parameters.fieldNames{f}); % Add empty field if not found
                localStruct.(parameters.fieldNames{f}) = [];
            else
                localStruct.(parameters.fieldNames{f}) = ...
                    parameters.fieldConv{f}(localStruct.(parameters.fieldNames{f}));
            end
        end
        localStruct.name = fileData(i).name;
        localStruct.filePath = [folderPath fileData(i).name];
        localStruct.regExp = parameters.regExp;

        parsedFileStruct= [parsedFileStruct localStruct];        
    else
        % Split text
        [splitText, foundDelimiters] = strsplit(fileData(i).name, parameters.delimiters);

        % Remove extension
        splitText = splitText(1:(end-1));
        foundDelimiters = foundDelimiters(1:(end-1));

        % Combine final field if appropriate
        if length(splitText) > length(parameters.fieldNames) && parameters.appendExtraFields
            % Identify field with internal delimiters and recombine split string
            combinedString = {};
            lengthDiff = length(splitText) - length(parameters.fieldNames) + 1;
            startInd = parameters.containsDelimiters;
            finishInd = length(splitText) - (length(parameters.fieldNames)- startInd);
            combinedString(1:2:2*lengthDiff) = splitText(startInd:finishInd);
            combinedString(2:2:(2*lengthDiff-1)) = foundDelimiters(startInd:(finishInd-1));

            % Replace split text with the revised values
            oldSplitText = splitText;
            splitText = {};
            splitText(1:(parameters.containsDelimiters-1)) = oldSplitText(1:(parameters.containsDelimiters-1));
            splitText{parameters.containsDelimiters} = [combinedString{:}];
            splitText((parameters.containsDelimiters+1):(length(parameters.fieldNames))) = ...
                oldSplitText((parameters.containsDelimiters+lengthDiff):end);
        end

        if parameters.requireExactMatch
            splitCondition = length(splitText) == length(parameters.fieldNames);
        else
            splitCondition = length(splitText) <= length(parameters.fieldNames);
        end

        % Parse split text
        if splitCondition
            parsedFileStruct(count).name = fileData(i).name;
            parsedFileStruct(count).filePath = [folderPath fileData(i).name];
            for j=1:length(splitText)
                parsedFileStruct(count).(parameters.fieldNames{j}) = ...
                    parameters.fieldConv{j}(splitText{j});
            end
            parsedFileStruct(count).delimiters = parameters.delimiters;
            count = count + 1;
        end
    end
end

% -------------------------------------------------------------------------
% Issue warning if some files did not fit the pattern
% -------------------------------------------------------------------------
if length(parsedFileStruct) < length(fileData)
    warning('matlabFunctions:unparsedFiles', 'Some files did not fit the specified pattern.');
end
