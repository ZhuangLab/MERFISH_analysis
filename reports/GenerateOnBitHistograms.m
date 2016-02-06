function [figHandles, parameters] = GenerateOnBitHistograms(words, varargin)
% ------------------------------------------------------------------------
% [figHandles, parameters] = GenerateOnBitHistograms(words, varargin)
% This function creates on bit histograms for the provided words for all
% words or by cell. 
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
% October 3, 2014
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.

% -------------------------------------------------------------------------
% Default Reports to Generate
% -------------------------------------------------------------------------
defaultReports = cell(0,2);
defaultReports(end+1,:) = {'numOnBitsHistByCell', 'on'}; %
defaultReports(end+1,:) = {'numOnBitsHistAllCells', 'on'}; %

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

figHandles = [];

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if (nargin < 1 || ...
        ~isstruct(words) || ...
        ~isempty(setdiff(fields(CreateWordsStructure(0,0)), fields(words)) ))
    error('matlabFunctions:invalidArguments', 'Invalid words structure.');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Cell By Cell Report
% -------------------------------------------------------------------------
reportID = find(strcmp('numOnBitsHistByCell', parameters.reportsToGenerate(:,1)));
if ~isempty(reportID)
    cellIDs = unique([words.cellID]);
    for i=1:length(cellIDs)
        figHandles(i) = figure('Name',['numOnBitsHist_cell', num2str(cellIDs(i))], ...
            'visible', parameters.reportsToGenerate{reportID, 2});
        localWords = words([words.cellID] == cellIDs(i));
        hist([localWords.numOnBits], 0:localWords(1).numHyb);
        xlabel('Number of On Bits');
        xlim([0 localWords(1).numHyb + 1]);
        ylabel('Counts');
        title(['Cell ' num2str(cellIDs(i))]);
        PresentationPlot();

        if parameters.saveAndClose
            if parameters.useSubFolderForCellReport
                SaveFigure(figHandles(i),'overwrite',parameters.overwrite,...
                    'formats',parameters.figFormats, ...
                    'subFolder', ['Cell_' num2str(cellIDs(i))]);
            else
                SaveFigure(figHandles(i),'overwrite',parameters.overwrite,...
                    'formats',parameters.figFormats);
            end
            close(figHandles(i));
        end
    end
end

% -------------------------------------------------------------------------
% All Cell Reports
% -------------------------------------------------------------------------
reportID = find(strcmp('numOnBitsHistAllCells', parameters.reportsToGenerate(:,1)));
if ~isempty(reportID)
    figHandles(end+1) = figure('Name','numOnBitsHistAllCells', ...
        'visible', parameters.reportsToGenerate{reportID, 2});
    hist([words.numOnBits], 0:words(1).numHyb);
    xlabel('Number of On Bits');
    xlim([0 words(1).numHyb + 1]);
    ylabel('Counts');
    PresentationPlot();
    
    if parameters.saveAndClose
        SaveFigure(figHandles(end),'overwrite',parameters.overwrite,...
            'formats',parameters.figFormats);

        close(figHandles(end));
    end
end
