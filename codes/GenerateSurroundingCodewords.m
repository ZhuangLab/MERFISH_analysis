function codewords = GenerateSurroundingCodewords(codeword, hammDist, varargin)
% ------------------------------------------------------------------------
% codewords = GenerateSurroundingCodewords(codeword, hammDist, varargin)
% This function generates a cell array of all codewords exactly the hamming
% distance (hammDist) of the specified codeword. 
%
% The codeword can be provided as a logical array, a string, or an integer.
%--------------------------------------------------------------------------
% Necessary Inputs
% --codebook/The codeword used to find surrounding codewords. 
% --hammDist/The hamming distance between the 'central' codeword and the 
%   returned codewords. 
% 
%--------------------------------------------------------------------------
% Outputs
% --codewords/A cell array of 1xN logicals corresponding to the surrounding
%   codewords
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 7, 2014
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if any(nargin < [1 2])
    error('matlabFunctions:invalidArguments', ...
        'Both a codeword and a hamming distance are required.');
end

if length(hammDist) ~= 1
    error('matlabFunctions:invalidArguments', 'Incorrect hamming distance');
end

if ~islogical(codeword)
    error('matlablFunctions:invalidArguments', 'Codeword must be a logical');
end

% -------------------------------------------------------------------------
% Allow logical output (faster) 
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'logical', 'boolean', false};

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);


% -------------------------------------------------------------------------
% Generate surrounding words
% -------------------------------------------------------------------------
C = nchoosek(1:length(codeword), hammDist);
codewords = repmat(codeword, [length(C) 1]);

for i=1:length(C)
    codewords(i, C(i,:)) = ~codewords(i, C(i,:));
end

% -------------------------------------------------------------------------
% Convert to cell
% -------------------------------------------------------------------------
if ~parameters.logical
    codewords = num2cell(codewords, 2);
end

