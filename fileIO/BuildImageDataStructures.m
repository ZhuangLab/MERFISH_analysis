function [imageData, parameters] = BuildImageDataStructures(folderPath, varargin)
% ------------------------------------------------------------------------
% [imageData] = BuildImageDataStructures(folderPath, varargin)
% This function creates imageData structures for all files in the
% folderPath that satisfy the provided criteria.  It is a wrapper for the 
% function BuildFileStructure.
%--------------------------------------------------------------------------
% Necessary Inputs
% folderPath/ A path to the desired folder
%--------------------------------------------------------------------------
% Outputs
% imageData/ A structure array containing a structure for each file found
% in the folder.  See CreateImageDataStructure for field information.
%--------------------------------------------------------------------------
% Variable Inputs
% See BuildFileStructure
%-------------------------------------------------------------------------- 
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 26, 2014
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
defaults(end+1,:) = {'appendExtraFields', 'boolean', false}; % Conversion functions
defaults(end+1,:) = {'requireFlag', 'string', ''}; % Conversion functions

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
% Find Files
% -------------------------------------------------------------------------
foundFiles = BuildFileStructure(folderPath, 'parameters', parameters);

% -------------------------------------------------------------------------
% Transfer undefined fields
% -------------------------------------------------------------------------
defaultImageDataStruct = CreateImageDataStructure(1);
extraFields = setdiff(fieldnames(defaultImageDataStruct), fieldnames(foundFiles));

imageData = foundFiles;
for i=1:length(imageData)
    for j=1:length(extraFields)
        imageData(i).(extraFields{j}) = defaultImageDataStruct.(extraFields{j});
    end
end

% -------------------------------------------------------------------------
% Reorder field names
% -------------------------------------------------------------------------
imageData = orderfields(imageData, defaultImageDataStruct);
