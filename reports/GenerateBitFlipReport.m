function [bitFlipReport, parameters] = GenerateBitFlipReport(words, exactMap, varargin)
% ------------------------------------------------------------------------
% [figHandles, parameters] = GenerateBitFlipReport(words, varargin)
% This function creates an error report for bit flip probabilities
%--------------------------------------------------------------------------
% Necessary Inputs
% words/A structure array with an element for each word. See
%   CreateWordsStructure for information about fields. 
%--------------------------------------------------------------------------
% Outputs
% figHandles/Handles for the generated figures. 
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Jeffrey Moffitt 
% lmoffitt@mcb.harvard.edu
% October 20, 2014
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.

% -------------------------------------------------------------------------
% Default Reports to Generate
% -------------------------------------------------------------------------
defaultReports = cell(0,2);
defaultReports(end+1,:) = {'bitFlipProbabilitiesAverage', 'on'}; %
defaultReports(end+1,:) = {'bitFlipProbabilitiesAllGenes', 'on'}; %

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'reportsToGenerate','cell', defaultReports};
defaults(end+1,:) = {'printedUpdates', 'boolean', true};
defaults(end+1,:) = {'saveAndClose', 'boolean', false};
defaults(end+1,:) = {'useSubFolderForCellReport', 'boolean', true};
defaults(end+1,:) = {'overwrite', 'boolean', true};
defaults(end+1,:) = {'figFormats', 'cell', {'png', 'fig'}};
defaults(end+1,:) = {'probToUse', {'exact', 'firstOrder'}, 'exact'};
defaults(end+1,:) = {'errCorrFunc', 'function', []};
defaults(end+1,:) = {'numHybs', 'nonnegative', 16};

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Defined required words fields
% -------------------------------------------------------------------------
requiredWordsFields = {'intCodeword'};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if (nargin < 1 || ...
        ~isstruct(words) || ...
        ~isempty(setdiff(requiredWordsFields, fields(words))))
    error('matlabFunctions:invalidArguments', 'Words is missing some required fields.');
end

% -------------------------------------------------------------------------
% Check error check function
% -------------------------------------------------------------------------
if isempty(parameters.errCorrFunc)
    warning('matlabFunctions:invalidArguments', ...
        'An error correction function must be provided to estimate bit flip probabilities');
    bitFlipReport = [];
    return
end

% -------------------------------------------------------------------------
% Normalize codeword type
% -------------------------------------------------------------------------
geneNames = values(exactMap);
codewords = keys(exactMap);

if ischar(codewords{1})
    codewords = cellfun(@(x)(x=='1'), codewords, 'UniformOutput', false);
else
    codewords = cellfun(@(x)de2bi(x,parameters.numHybs), codewords, 'UniformOutput', false);
end

% -------------------------------------------------------------------------
% Compute Word Histogram
% -------------------------------------------------------------------------
[n, x] = hist([words.intCodeword], 1:2^parameters.numHybs);

% -------------------------------------------------------------------------
% Generate bit flip probabilities
% -------------------------------------------------------------------------
firstOrderbitFlipProbabilities = nan(length(geneNames), parameters.numHybs, 2);
numCounts = zeros(1, length(codewords));
for i=1:length(codewords)
    % Identify words
    correctWord = bi2de(fliplr(codewords{i}));
    surroundingWords = [cellfun(@(x)bi2de(fliplr(x)), parameters.errCorrFunc(codewords{i}))'];
    
    % Compute first order estimates of probability: Assumes that no more
    % than one error is possible
    totalCounts = sum(n([correctWord surroundingWords]));
    localFirstOrderProb = n(surroundingWords)/totalCounts;  
    numCounts(i) = totalCounts;
    
    % Compute exact probability: assumes that any number of errors are
    % present but that no cross contamination between words occurs
    alpha = n(surroundingWords)/n(correctWord); 
    localExactProb = alpha./(1+alpha);
    
    % Identify different transition types
    oneToZeroInds = find(surroundingWords < correctWord);
    zeroToOneInds = find(surroundingWords > correctWord);
    
    % Sort on order (left bits > right bits)
    [~, oneToZeroSind] = sort(surroundingWords(oneToZeroInds), 'ascend');
    [~, zeroToOneSind] = sort(surroundingWords(zeroToOneInds), 'ascend');

    % Record probabilities: first order
    firstOrderbitFlipProbabilities(i, oneToZeroInds(oneToZeroSind), 1) = ...
        localFirstOrderProb(oneToZeroInds(oneToZeroSind));
    firstOrderbitFlipProbabilities(i, zeroToOneInds(zeroToOneSind), 2) = ...
        localFirstOrderProb(zeroToOneInds(zeroToOneSind));
    
    % Record probabilities: exact
    exactBitFlipProbabilities(i, oneToZeroInds(oneToZeroSind), 1) = ...
        localExactProb(oneToZeroInds(oneToZeroSind));
    exactBitFlipProbabilities(i, zeroToOneInds(zeroToOneSind), 2) = ...
        localExactProb(zeroToOneInds(zeroToOneSind));    
end

% -------------------------------------------------------------------------
% Archive results
% -------------------------------------------------------------------------
bitFlipReport.geneNames = geneNames;
bitFlipReport.counts = numCounts;
bitFlipReport.firstOrderProbabilities = firstOrderbitFlipProbabilities;
bitFlipReport.exactProbabilities = exactBitFlipProbabilities;

% -------------------------------------------------------------------------
% Determine method for calculating probability
% -------------------------------------------------------------------------
switch parameters.probToUse
    case 'exact'
        bitFlipReport.probabilities = bitFlipReport.exactProbabilities;
    case 'firstOrder'
        bitFlipReport.probabilities = bitFlipReport.firstOrderProbabilities;
end

% Average probabilities scaled by gene
bitFlipReport.hybProb = squeeze(nanmean(bitFlipReport.probabilities,1));
bitFlipReport.hybProbErr = squeeze(nanstd(bitFlipReport.probabilities,1));
bitFlipReport.numCounts = numCounts;

% Compute weights
weigthVec = numCounts'/sum(numCounts);
weights = repmat(weigthVec, [1 parameters.numHybs 2]);

% Average probabilities scaled by counts per gene
bitFlipReport.scaledHybProb = squeeze(nansum(bitFlipReport.probabilities.*weights,1));
bitFlipReport.scaledHybProbErr = squeeze(sqrt(nanvar(bitFlipReport.probabilities, weigthVec, 1)))/ ...
    sqrt(length(weigthVec));

% -------------------------------------------------------------------------
% Display Probabilites
% -------------------------------------------------------------------------
figCount = 1;
reportID = find(strcmp('bitFlipProbabilitiesAverage', parameters.reportsToGenerate(:,1)));
if ~isempty(reportID)
    bitFlipReport.figHandles(figCount) = figure('Name',['bitFlipProbs'], ...
        'visible', parameters.reportsToGenerate{reportID, 2});
    titles = {'1->0', '0->1'};
    for i=1:2
        subplot(1,2,i);
        bar(1:parameters.numHybs, squeeze(bitFlipReport.scaledHybProb(:,i))); hold on;
        errorbar(1:parameters.numHybs, squeeze(bitFlipReport.scaledHybProb(:,i)), ...
            squeeze(bitFlipReport.scaledHybProbErr(:,i)), 'k.');
        xlabel('Hyb');
        ylabel('Probability');
        title(titles{i});
        xlim([0 parameters.numHybs + 1]);
    end
    
    PresentationPlot();

    if parameters.saveAndClose
        SaveFigure(bitFlipReport.figHandles(figCount),'overwrite',parameters.overwrite,...
            'formats',parameters.figFormats);
        close(bitFlipReport.figHandles(figCount));
    end

    figCount = figCount + 1;
end

% -------------------------------------------------------------------------
% Display Probabilites for All Genes
% -------------------------------------------------------------------------
reportID = find(strcmp('bitFlipProbabilitiesAllGenes', parameters.reportsToGenerate(:,1)));
if ~isempty(reportID)
    bitFlipReport.figHandles(figCount) = figure('Name',['bitFlipProbsAllGenes'], ...
        'visible', parameters.reportsToGenerate{reportID, 2});
    titles = {'1->0', '0->1'};
    for i=1:2
        subplot(1,2,i);
        imagesc(squeeze(bitFlipReport.probabilities(:,:,i)));
        xlabel('Hyb');
        ylabel('Genes');
        title(titles{i});
        colorbar;
        inds = cellfun(@str2num, cellstr(get(gca, 'YTickLabel')));
        inds = inds(inds>1 & inds <= length(bitFlipReport.geneNames));
        set(gca, 'YTickLabel', bitFlipReport.geneNames(inds));
    end
    
    PresentationPlot('FontSize', 10);

    if parameters.saveAndClose
        SaveFigure(bitFlipReport.figHandles(figCount),'overwrite',parameters.overwrite,...
            'formats',parameters.figFormats);
        close(bitFlipReport.figHandles(figCount));
    end

    figCount = figCount + 1;
end
