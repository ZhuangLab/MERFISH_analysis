function [FPKMReport, parameters] = GenerateFPKMReport(words, FPKMData, varargin)
% ------------------------------------------------------------------------
% [FPKMData, parameters] = GenerateFPKMReport(words, FPKMData, varargin)
% This function generates an FPKM report from the provided words.
%--------------------------------------------------------------------------
% Necessary Inputs
% words/A structure array with an element for each found word. 
%   See CreateWordsStructure for information on field names
% FPKMdata/A structure array with the following fields
%   --geneName: A string specifying the name of the RNA isform
%   --FPKM: The FPKM for that gene
%--------------------------------------------------------------------------
% Outputs
% FPKMData/A structure containing information from the report
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Jeffrey Moffitt 
% lmoffitt@mcb.harvard.edu
% September 26, 2014
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.

% -------------------------------------------------------------------------
% Default Reports to Generate
% -------------------------------------------------------------------------
defaultReports = cell(0,2);
defaultReports(end+1,:) = {'totalFPKMReport', 'on'}; % {report name, 'off'/'on' do not/do display figure}
%defaultReports(end+1,:) = {'cellByCellFPKMReport', 'on'}; Not included as
%a default report

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'reportsToGenerate','cell', defaultReports};
defaults(end+1,:) = {'printedUpdates', 'boolean', true};
defaults(end+1,:) = {'FPKMReportExactMatchOnly', 'boolean', false};
defaults(end+1,:) = {'FPKMReportEmbedNames', 'boolean', true};
defaults(end+1,:) = {'showNames', 'boolean', false};
defaults(end+1,:) = {'embedNames', 'boolean', false};
defaults(end+1,:) = {'saveAndClose', 'boolean', false};
defaults(end+1,:) = {'useSubFolderForCellReport', 'boolean', true};
defaults(end+1,:) = {'overwrite', 'boolean', true};
defaults(end+1,:) = {'figFormats', 'cell', {'png', 'fig'}};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if (nargin < 2 || ...
        ~isstruct(words) || ...
        ~ismember('geneName', fields(words)) ) || ...
        ~isstruct(FPKMData) || ...
        ~all(ismember({'geneName', 'FPKM'}, fields(FPKMData)))
    error('matlabFunctions:invalidArguments', 'Invalid words or FPKMdata structures.');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Identify valid words
% -------------------------------------------------------------------------
if parameters.FPKMReportExactMatchOnly
    validWords = [words.isExactMatch];
    figNameModifier = 'ExactOnly';
else
    validWords = [words.isExactMatch] | [words.isCorrectedMatch];
    figNameModifier = 'ExactAndCorrected';
end

if parameters.printedUpdates
    display('---------------------------------------------------------------');
    display(['Found ' num2str(sum(validWords)) ' valid words in ' ...
        num2str(length(words)) ' total words']);
end

FPKMReport.wordInds = find(validWords); % Save for report

% -------------------------------------------------------------------------
% Determine target gene names
% -------------------------------------------------------------------------
targetNames = {FPKMData.geneName};

% -------------------------------------------------------------------------
% Generate total report
% -------------------------------------------------------------------------
reportID = find(strcmp('totalFPKMReport', parameters.reportsToGenerate(:,1)));
if ~isempty(reportID)
    
    FPKMReport.totalReportFigHandle = figure(...
        'Name', ['FPKMReportAllCells_' figNameModifier], ...
        'visible', parameters.reportsToGenerate{reportID,2});
    
    counts = GenerateGeneCounts({words(validWords).geneName}, ...
        targetNames); % Record number of names not in list

    FPKMReport.geneNames = targetNames;
    FPKMReport.geneNames{end+1} = 'unknown name';
    FPKMReport.counts = counts;
    FPKMReport.countsWOUnknown = counts(1:(end-1));
    FPKMReport.FPKM = [FPKMData.FPKM];
    
    if parameters.showNames
        FPKMReport.pearsonCorr = PlotCorr2([FPKMData.FPKM], counts(1:(end-1)), ...
            'figHandle', FPKMReport.totalReportFigHandle, ...
            'pointNames', targetNames);
        set(FPKMReport.totalReportFigHandle, 'Name', ...
            [get(FPKMReport.totalReportFigHandle, 'Name') 'withNames']);
    else
        FPKMReport.pearsonCorr = PlotCorr2([FPKMData.FPKM], counts(1:(end-1)), ...
            'figHandle', FPKMReport.totalReportFigHandle);
    end
    if parameters.embedNames
        ImbedNamesInFigure(FPKMReport.totalReportFigHandle, targetNames);
    end
    xlabel('FPKM');
    ylabel('Counts');
    PresentationPlot();
    
    if parameters.saveAndClose
        SaveFigure(FPKMReport.totalReportFigHandle, 'overwrite', true, ...
            'formats', {'fig', 'png'});
        close(FPKMReport.totalReportFigHandle);
    end
    
    if parameters.printedUpdates
        display(['All Cells p log10: ' num2str(FPKMReport.pearsonCorr.log10rho)]);
    end
end

% -------------------------------------------------------------------------
% Generate Individual Cell Report
% -------------------------------------------------------------------------
reportID = find(strcmp('cellByCellFPKMReport', parameters.reportsToGenerate(:,1)));
if ~isempty(reportID)
    cellIDs = unique([words.cellID]);
    for c=1:length(cellIDs)
        FPKMReport.cellReportFigHandles(c) = figure(...
            'Name', ['FPKMReportCell_' num2str(cellIDs(c))], ...
            'visible', parameters.reportsToGenerate{reportID,2});
        
        localInds = validWords & [words.cellID] == cellIDs(c);
        counts = GenerateGeneCounts({words(localInds).geneName}, ...
            {FPKMData.geneName});
        FPKMReport.cellReportCounts(c,:) = counts;
        
        
        if parameters.showNames
            FPKMReport.cellByCellPearsonCorr(c) = PlotCorr2([FPKMData.FPKM], counts(1:(end-1)), ...
                'figHandle', FPKMReport.cellReportFigHandles(c), ...
                'pointNames', targetNames);
        else
            FPKMReport.cellByCellPearsonCorr(c) = PlotCorr2([FPKMData.FPKM], counts(1:(end-1)), ...
                'figHandle', FPKMReport.cellReportFigHandles(c));
        end
        if parameters.embedNames
            ImbedNamesInFigure(FPKMReport.cellReportFigHandles(c), targetNames);
        end
        xlabel('FPKM');
        ylabel('Counts');
        PresentationPlot();

        if parameters.saveAndClose
            if parameters.useSubFolderForCellReport
                SaveFigure(FPKMReport.cellReportFigHandles(c), 'overwrite', parameters.overwrite, ...
                    'formats', parameters.figFormats, ...
                    'subFolder', ['Cell_' num2str(cellIDs(c))]);
            else
                SaveFigure(FPKMReport.cellReportFigHandles(c), 'overwrite', parameters.overwrite, ...
                    'formats', parameters.figFormats);
            end
            close(FPKMReport.cellReportFigHandles(c));
        end
        
        if parameters.printedUpdates
            display(['Cell ' num2str(cellIDs(c)) ' p log10: ' num2str(FPKMReport.cellByCellPearsonCorr(c).log10rho)]);
        end
    end
end

end

% -------------------------------------------------------------------------
% Internal Function Definitions: Simple name-based histogram
% -------------------------------------------------------------------------
function counts = GenerateGeneCounts(geneNames, targetNames)

    counts = zeros(1, length(targetNames)+1);
    for j=1:length(targetNames)
        inds = strcmp(geneNames, targetNames{j});
        counts(j) = sum(inds);
        geneNames = geneNames(~inds);
    end
    counts(end) = length(geneNames);
end