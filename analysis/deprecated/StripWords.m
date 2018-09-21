function words = StripWords(words, varargin)
% ------------------------------------------------------------------------
% words = StripWords(words, varargin)
% StripWords(words) returns a truncated from of the words structure
% containing either a default set of fields (see below) or a specified set
% of fields. 
%--------------------------------------------------------------------------
% Necessary Inputs
% words/A word data structure: See CreateWordsStructure. 
%--------------------------------------------------------------------------
% Outputs
% strippedWords/ A word data structure with all but a defined set of fields
% removed. 
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
% fieldsToKeep/cell array of strings: A list of the fields to keep in words
%   The default fields are: intCodeword, geneName, isExactMatch, isCorrectedMatch, imageX,
%   imageY, cellID, wordCentroidX, wordCentroidY
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% October 20, 2014
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    dataPath = '\\cajal\TSTORMdata\Datasets.xls';
end

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'fieldsToKeep', 'cell', ...
    {'intCodeword', 'geneName', 'isExactMatch', 'isCorrectedMatch', ...
    'imageX', 'imageY', 'wordCentroidX', 'wordCentroidY', 'cellID'}}; 

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if (nargin < 1 || ...
        ~isstruct(words) || ...
        ~isempty(setdiff(fields(CreateWordsStructure(0,0)), fields(words)) ))
    error('matlabFunctions:invalidArguments', 'Invalid words structure.');
end

% -------------------------------------------------------------------------
% Check fieldsToKeep
% -------------------------------------------------------------------------
if ~isempty(setdiff(parameters.fieldsToKeep, fields(words)))
    warning('matlabFunctions:invalidArguments', 'Some desired fields are not present');
end

% -------------------------------------------------------------------------
% Truncate words
% -------------------------------------------------------------------------
parameters.fieldsToKeep = parameters.fieldsToKeep(ismember(parameters.fieldsToKeep, ...
    fields(words)));

fieldsToRemove = setdiff(fields(words), parameters.fieldsToKeep);

words = rmfield(words, fieldsToRemove);