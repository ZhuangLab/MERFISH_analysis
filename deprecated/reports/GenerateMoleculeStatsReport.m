function [moleculeStats, parameters] = GenerateMoleculeStatsReport(words, varargin)
% ------------------------------------------------------------------------
% moleculeStats = GenerateMoleculeStatsReport(words, varargin)
% This function generates a series of basic statistics for all words and
% also cell by cell.  
%--------------------------------------------------------------------------
% Necessary Inputs
% words/A structure array with an element for each found word. 
%   See CreateWordsStructure for information on field names
%--------------------------------------------------------------------------
% Outputs
% moleculeStats/A structure containing information from the report
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Jeffrey Moffitt 
% lmoffitt@mcb.harvard.edu
% October 1, 2014
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.

% -------------------------------------------------------------------------
% Default Reports to Generate
% -------------------------------------------------------------------------
defaultReports = cell(0,2);
defaultReports(end+1,:) = {'molStats', 'on'}; % {report name, 'off'/'on' do not/do display figure}
defaultReports(end+1,:) = {'molHist', 'on'};
defaultReports(end+1,:) = {'molDistStats', 'on'};

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

defaults(end+1,:) = {'brightHistBins', 'positive', 100};
defaults(end+1,:) = {'distHistBins', 'positive', 100};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if (nargin < 1 || ...
        ~isstruct(words) || ...
        ~isempty(setdiff(fields(CreateWordsStructure(0,0)), fields(words))))
    error('matlabFunctions:invalidArguments', 'Invalid words structure.');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Generate basic statistics: Molecule Properties
% -------------------------------------------------------------------------
mListFields = {'a', 'bg', 'h'};
numHyb = words(1).numHyb;
sem = @(x)std(x)/sqrt(length(x));
for j=1:length(mListFields)
    data = [words.(mListFields{j})];
    for i=1:numHyb
        localData = data(i:numHyb:end);
        localData = localData(~isnan(localData));
        [n, x] = hist(localData, parameters.brightHistBins);
        moleculeStats.([mListFields{j} 'Hist'])(i,1,:) = x;
        moleculeStats.([mListFields{j} 'Hist'])(i,2,:) = n;
        moleculeStats.([mListFields{j} 'Mean'])(i) = mean(localData);
        moleculeStats.([mListFields{j} 'STD'])(i) = std(localData);
        moleculeStats.([mListFields{j} 'SEM'])(i) = sem(localData);
        moleculeStats.([mListFields{j} 'Median'])(i) = median(localData);
        moleculeStats.([mListFields{j} 'IQR'])(i) = iqr(localData);
        moleculeStats.([mListFields{j} 'N'])(i) = length(localData);
    end
end

% -------------------------------------------------------------------------
% Generate basic statistics: Molecule distances from word center
% -------------------------------------------------------------------------
xPos = [words([words.numOnBits]>1).xc];
xPos(xPos==0) = nan;
xDist = reshape(xPos, [numHyb sum([words.numOnBits]>1)]) - ...
    repmat([words([words.numOnBits]>1).wordCentroidX], [numHyb 1]);

yPos = [words([words.numOnBits]>1).yc];
yPos(yPos==0) = nan;
yDist = reshape(yPos, [numHyb sum([words.numOnBits]>1)]) - ...
    repmat([words([words.numOnBits]>1).wordCentroidY], [numHyb 1]);

moleculeStats.xDist = xDist;
moleculeStats.yDist = yDist;

posNames = {'x', 'y'};
for i=1:numHyb
    for j=1:length(posNames)
        data = moleculeStats.([posNames{j} 'Dist'])(i,:);
        moleculeStats.([posNames{j} 'Mean'])(i) = mean(data(~isnan(data)));
        moleculeStats.([posNames{j} 'STD'])(i) = std(data(~isnan(data)));
        moleculeStats.([posNames{j} 'MeanAbs'])(i) = mean(abs(data(~isnan(data))));
        [n, x] = hist(data(~isnan(data)), parameters.distHistBins);

        moleculeStats.([posNames{j} 'Hist'])(i,1,:) = x;
        moleculeStats.([posNames{j} 'Hist'])(i,2,:) = n;
    end
end

% -------------------------------------------------------------------------
% Plot Molecule Stats
% -------------------------------------------------------------------------
reportID = find(strcmp('molStats', parameters.reportsToGenerate(:,1)));
if ~isempty(reportID)
    moleculeStats.molStatsFigHandle = figure(...
        'Name', 'molStats', ...
        'visible', parameters.reportsToGenerate{reportID,2});
    for i=1:length(mListFields)
        subplot(length(mListFields), 2, 2*(i-1)+1);
        errorbar(1:numHyb, moleculeStats.([mListFields{i} 'Mean']), ...
            moleculeStats.([mListFields{i} 'SEM']), 'b.');
        xlabel('Hybe Number');
        xlim([0 numHyb+1]);
        ylabel(mListFields{i});
        
        subplot(length(mListFields), 2, 2*(i));
        bar(1:numHyb, moleculeStats.([mListFields{i} 'N']));
        xlabel('Hybe Number');
        ylabel('Counts');
        xlim([0 numHyb+1]);
    end
        
    PresentationPlot();
    
    if parameters.saveAndClose
        SaveFigure(moleculeStats.molStatsFigHandle, 'overwrite', parameters.overwrite, ...
                'formats', parameters.figFormats);
        close(moleculeStats.molStatsFigHandle);
        moleculeStats.molStatsFigHandle = [];
    end
end

% -------------------------------------------------------------------------
% Plot Molecule Stats: Position
% -------------------------------------------------------------------------
reportID = find(strcmp('molDistStats', parameters.reportsToGenerate(:,1)));
if ~isempty(reportID)
    moleculeStats.molDistStatsFigHandle = figure(...
        'Name', 'molDistStats', ...
        'visible', parameters.reportsToGenerate{reportID,2});
    for i=1:length(posNames)
        subplot(length(posNames), 2, 2*(i-1)+1);
        errorbar(1:numHyb, moleculeStats.([posNames{i} 'Mean']), ...
            moleculeStats.([posNames{i} 'STD']), 'b.');
        xlabel('Hybe Number');
        xlim([0 numHyb+1]);
        ylabel(posNames{i});
        title('Average/STD Distances');
        
        subplot(length(posNames), 2, 2*(i));
        localData = moleculeStats.([posNames{i} 'Dist']);
        localData = localData(~isnan(localData));
        [n, x] = hist(localData, parameters.distHistBins);
        stairs(x,n);
        xlabel('Position');
        ylabel('Counts');
        title([posNames]);
    end
        
    PresentationPlot();
    
    if parameters.saveAndClose
        SaveFigure(moleculeStats.molDistStatsFigHandle, 'overwrite', parameters.overwrite, ...
                'formats', parameters.figFormats);
        close(moleculeStats.molDistStatsFigHandle);
    end
end

% -------------------------------------------------------------------------
% Plot Molecule Histograms
% -------------------------------------------------------------------------
reportID = find(strcmp('molHist', parameters.reportsToGenerate(:,1)));
if ~isempty(reportID)
    for i=1:length(mListFields)
        moleculeStats.molHistFigHandle(i) = figure(...
            'Name', ['molHist_' mListFields{i}], ...
            'visible', parameters.reportsToGenerate{reportID,2}, ...
            'Position', [8 51 1076 635]);
   
        maxValueX = max(moleculeStats.([mListFields{i} 'Hist'])(:,1,end));
        maxValueY = max(max(moleculeStats.([mListFields{i} 'Hist'])(:,2,:)));
        for j=1:numHyb
            subplot(4, ceil(numHyb/4), j);
            stairs(squeeze(moleculeStats.([mListFields{i} 'Hist'])(j,1,:)), ...
                squeeze(moleculeStats.([mListFields{i} 'Hist'])(j,2,:)));
            xlabel(mListFields{i});
            xlim([0 maxValueX]);
            ylabel('Counts');
            ylim([0 1.5*maxValueY]);
            title(['Hyb ' num2str(j)]);
        end
        PresentationPlot('FontSize', 10);
    
        if parameters.saveAndClose
            SaveFigure(moleculeStats.molHistFigHandle(i), 'overwrite', parameters.overwrite, ...
                'formats', parameters.figFormats);
            close(moleculeStats.molHistFigHandle(i));
        end
    end
end

% -------------------------------------------------------------------------
% Plot Distances
% -------------------------------------------------------------------------
