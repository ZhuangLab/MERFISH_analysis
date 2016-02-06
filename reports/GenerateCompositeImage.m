function [figHandles, parameters] = GenerateCompositeImage(words, imageData, varargin)
% ------------------------------------------------------------------------
% [figHandles, parameters] = GenerateCompositeImage(words, imageData, varargin)
% This function creates a composite image
%--------------------------------------------------------------------------
% Necessary Inputs
% imageData/A structure array with an element for each image used to create
%   the elements in words.  See CreateWordsStructure for information on
%   field names. 
%--------------------------------------------------------------------------
% Outputs
% hybReport/A structure containing information from the report
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Jeffrey Moffitt 
% lmoffitt@mcb.harvard.edu
% October 2, 2014
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.

% -------------------------------------------------------------------------
% Default Reports to Generate
% -------------------------------------------------------------------------
defaultReports = cell(0,2);
defaultReports(end+1,:) = {'cellWithWordsImage', 'on'}; % {report name, 'off'/'on' do not/do display figure}

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'reportsToGenerate','cell', defaultReports};
defaults(end+1,:) = {'printedUpdates', 'boolean', true};
defaults(end+1,:) = {'showNames', 'boolean', false};
defaults(end+1,:) = {'embedNames', 'boolean', true};
defaults(end+1,:) = {'saveAndClose', 'boolean', false};
defaults(end+1,:) = {'useSubFolderForCellReport', 'boolean', true};
defaults(end+1,:) = {'overwrite', 'boolean', true};
defaults(end+1,:) = {'figFormats', 'cell', {'png', 'fig'}};
defaults(end+1,:) = {'numImageColumns', 'positive', 4};
defaults(end+1,:) = {'displayHybLabel', 'boolean', true};
%defaults(end+1,:) = {'wordOverlayStyle', {'matchType', 'wordIdentity'}, 'matchType'};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if (nargin < 2 || ...
        ~isstruct(words) || ...
        ~isempty(setdiff(fields(StripWords(CreateWordsStructure(0,0))), fields(words))) || ...
        ~isstruct(imageData) || ...
        ~isempty(setdiff(fields(CreateImageDataStructure(0)), fields(imageData)) ))
    error('matlabFunctions:invalidArguments', 'Invalid words or imageData structures.');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Generate total report
% -------------------------------------------------------------------------
reportID = find(strcmp('cellWithWordsImage', parameters.reportsToGenerate(:,1)));
if ~isempty(reportID)
    cellIDs = unique([imageData.cellNum]);
    
    for i=1:length(cellIDs)
        % Index image Data and sort by hyb number
        localImageData = imageData([imageData.cellNum] == cellIDs(i));
        hybNums = [localImageData.hybNum];
        [~, sind] = sort(hybNums, 'Ascend');
        localImageData = localImageData(sind);
        
        % Create figure handle. 
        figHandles(i) = figure('Name', ['CellWithWords_Cell_' num2str(cellIDs(i))], ...
            'visible', parameters.reportsToGenerate{reportID,2});
        
        % Create aligned image
        alignedIm = zeros(localImageData(1).imageH, localImageData(1).imageW, length(localImageData));
        for h=1:length(localImageData)
            dax = max(ReadDax(localImageData(h).infFilePath,'endFrame',1, 'verbose', false),[],3); % max project the first frame
            [H,W] = size(dax); 
            tformInv = fliptform(localImageData(h).tform); 
            alignedDax = imtransform(dax,tformInv,...
                            'XYScale',1,'XData',[1 W],'YData',[1 H]);
            alignedIm(:,:,h) = alignedDax;
        end
        
        % Add color and plot image
        Io = Ncolor(uint16(alignedIm),gray(length(localImageData)));
        set(0, 'CurrentFigure', figHandles(i));
        imshow(Io, []); hold on;
        
        % Find local words
        localWords = words([words.cellID] == cellIDs(i));

        % Plot Exact matches
        inds = [localWords.isExactMatch];
        plot([localWords(inds).wordCentroidX], [localWords(inds).wordCentroidY], 'go'); hold on;
        
        % Plot corrected matches
        inds = [localWords.isCorrectedMatch];
        plot([localWords(inds).wordCentroidX], [localWords(inds).wordCentroidY], 'rx'); hold on;
        
        % Plot unmatched words
        inds = ~ ([localWords.isExactMatch] | [localWords.isCorrectedMatch]);
        plot([localWords(inds).wordCentroidX], [localWords(inds).wordCentroidY], 'b.'); hold on;
        
        PresentationPlot('MarkerWidth', 6, 'LineWidth', .5);
        
        if parameters.saveAndClose
                if parameters.useSubFolderForCellReport
                    SaveFigure(figHandles(i), 'overwrite', parameters.overwrite, ...
                        'formats', parameters.figFormats, ...
                        'subFolder', ['Cell_' num2str(cellIDs(i))]);
                else
                    SaveFigure(figHandles(i), 'overwrite', parameters.overwrite, ...
                        'formats', parameters.figFormats);
                end
            close(figHandles(i));
        end
    end
end
