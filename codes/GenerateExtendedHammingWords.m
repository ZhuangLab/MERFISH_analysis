function [words, gen, numParityBits] = GenerateExtendedHammingWords(numDataBits, varargin)
% ------------------------------------------------------------------------
% [words, generator, numParityBits] = GenerateExtendedHammingWords(numDataBits)
% This function returns the minimum extended hamming code words for the
% number of data bits specified. It also returns the generator matrix.
%--------------------------------------------------------------------------
% Necessary Inputs
%--------------------------------------------------------------------------
% Outputs
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% October 7, 2014
%
% Updated Dec 5, 2014
% removed loop over words to accelerate code.  (200x)
%
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1 || numDataBits < 1
    error('matlabFunctions:invalidArguments', 'A valid number of data bits is required.');
end

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'numOn','positive', 0};        % Returns a code with a fixed number of on bits
defaults(end+1,:) = {'parallel', 'parallel', []};   % A parallel.pool object can be provided to speed some calculations

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Handle number of parallel workers
% -------------------------------------------------------------------------
if isempty(parameters.parallel)
    numPar = 0;
else
    numPar = parameters.parallel.NumWorkers;
end

% -------------------------------------------------------------------------
% Determine the number of required parity bits for the hamming code
% -------------------------------------------------------------------------
numParityBits = ceil(fzero(@(x)2^(x) - x -1 - numDataBits, max(log2(numDataBits),1)));

% -------------------------------------------------------------------------
% Create the hamming code parity matrix
% -------------------------------------------------------------------------
par = hammgen(max(ceil(log2(numDataBits + numParityBits)), 3));

% -------------------------------------------------------------------------
% Create the extended hamming parity check matrix
% -------------------------------------------------------------------------
par(:,end+1) = zeros(size(par,1),1);
par(end+1,:) = ones(1, size(par,2));
par = mod(rref(par),2);

% -------------------------------------------------------------------------
% Create the generator
% -------------------------------------------------------------------------
gen = gen2par(par);

% -------------------------------------------------------------------------
% Create the shortened generator
% -------------------------------------------------------------------------
excessBits = size(gen,1) - numDataBits;
gen = gen(1:(end-excessBits), 1:(end-excessBits));

% -------------------------------------------------------------------------
% Create the words
% -------------------------------------------------------------------------
if parameters.numOn == 0 % Return all. Memory inefficient for large codes
    wordDec = (0:(2^numDataBits-1))';
    words = de2bi(wordDec,size(gen,1));
    words = mod(words*gen,2);
else
    % Determine number of data words
    numWords = 2^numDataBits;
    numOn = parameters.numOn; % Make local copy for smpd to copy
    
    % Decide the number of parallel workers to use
    numParToUse = min(numWords, numPar);
    
    % Initialize words for parallel processing
    words = Composite(numParToUse);
    for p=1:length(words)
        words{p} = zeros(0, size(gen,2));
    end
    % Loop over words and keep only those that have the correct number of
    % words using parallel processing if requested
    spmd (numParToUse)
        for i=labindex:numlabs:numWords
            % Covert to binary
            localWord = de2bi(i-1, size(gen,1));
            % Calculate extended hamming code word
            localWord = mod(localWord*gen,2);
            % Look for threshold
            if sum(localWord) == parameters.numOn
                words = cat(1,words, localWord);
            end
        end
    end
    
    % Construct words from all parallel workers
    words = cat(1, words{:});

end
    
