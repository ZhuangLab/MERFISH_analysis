function [words, parameters] = DecodeWords(words, exactMap, correctableMap, varargin)
% ------------------------------------------------------------------------
% words = DecodeWords(words, exactMap, correctableMap, varargin)
% This function decodes words based on the maps provided in the
%   containers.Map objects, exactMap and correctableMap. 
%--------------------------------------------------------------------------
% Necessary Inputs
% words/A structure array of found words. These structures must contain the
%   following fields:
%   --codeword: A logical array containing the desired word. 
% exactMap/A containers.Map object containing keys and values corresponding
%   to 'correct' codewords.  If only correctable matches are desired, then
%   [] can be passed for this argument. 
% correctableMap/A containers.Map object containing keys and values
%   corresponding to codewords that are not exact matches to entries in the
%   codebook but which can be associated with a codeword via some method,
%   e.g. an error correcting code.  If only exact matches are desired, then
%   [] can be passed for this argument. 
%--------------------------------------------------------------------------
% Outputs
% words/The same structure array with the addition of the following fields
%   --geneName: The string entry in the codebookMap corresponding to 
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
% --keyType ('int', 'binStr'): A string specifying the type of the key for
%   the providedcontainers.Map
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 8, 2014
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'keyType', {'int', 'binStr'}, 'binStr'}; % The type of key for the containers.Map

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 3
    error('matlabFunctions:invalidArguments', 'A word structure and two codebook maps must be provided.');
end

if ~isstruct(words) || ~ismember('codeword', fields(words))
    error('matlabFunctions:invalidArguments', 'The first argument must be a word structure or array');
end

if ~isempty(exactMap) & ~strcmp(class(exactMap), 'containers.Map')
    error('matlabFunctions:invalidArguments', 'The second argument must be a containers.Map object');
end

if ~isempty(correctableMap) & ~strcmp(class(correctableMap), 'containers.Map')
    error('matlabFunctions:invalidArguments', 'The third argument must either be empty or be a containers.Map object');
end

% -------------------------------------------------------------------------
% Determine Key Conversion Function
% -------------------------------------------------------------------------
removeWs = @(x) x(~isspace(x)); % Useful short hand for removing whitespace
switch parameters.keyType
    case 'int'
        keyConv= @(x) uint16(bin2dec(removeWs(num2str(x))));
    case 'binStr'
        keyConv = @(x) removeWs(num2str(x));
end

% -------------------------------------------------------------------------
% Decode words
% -------------------------------------------------------------------------
for i=1:length(words)
    % ---------------------------------------------------------------------
    % Decode exact matches
    % ---------------------------------------------------------------------
    if ~isempty(exactMap)
        try
            words(i).geneName = exactMap(keyConv(words(i).codeword));
            words(i).isExactMatch = true;
        end
    end
    % ---------------------------------------------------------------------
    % Decode exact matches
    % ---------------------------------------------------------------------
    if ~isempty(correctableMap)
        try
            words(i).geneName = correctableMap(keyConv(words(i).codeword));
            words(i).isCorrectedMatch = true;
        end
    end 
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
