function codewords = SECDEDCorrectableWords(codeword, varargin)
% ------------------------------------------------------------------------
% codewords = SECDEDCorrectableWords(codeword, varargin)
% This function generates a cell array of all codewords that would be
% corrected to the given codeword using a SECDED code, i.e. hamming
% distance = 1. 
%--------------------------------------------------------------------------
% Necessary Inputs
% --codebook/The codeword used to find surrounding codewords. 
% 
%--------------------------------------------------------------------------
% Outputs
% --codewords/A cell array of 1xN logicals corresponding to the surrounding
%   codewords
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 8, 2014
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabFunctions:invalidArguments', ...
        'A codeword must be provided.');
end

if ~islogical(codeword)
    error('matlablFunctions:invalidArguments', 'Codeword must be a logical');
end

% -------------------------------------------------------------------------
% Generate surrounding words
% -------------------------------------------------------------------------
codewords = GenerateSurroundingCodewords(codeword, 1);



