function [fiducialData, parameters] = AlignFiducials(fiducialData,varargin)
% ------------------------------------------------------------------------
% [fiducialData] = AlignFiducials(fiducialData,varargin)
% This function analyzes a series of raw conventional images in the 
%   specified directory and creates a words structure which represents all
%   of the identified words in the data.
%--------------------------------------------------------------------------
% Necessary Inputs
% fiducialData/A structure array with elements equal to the number of
%   images to align. This structure must contain an mList field. 
%--------------------------------------------------------------------------
% Outputs
% fiducialData/The same structure array provided but with two new fields
%   added.
%   --transform. The transform that will bring each image into the
%     reference frame of the first image. 
%   --warpUncertainty. A matrix that describes the quality of the generated
%     transform. 
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Alistair Boettiger
% boettiger@fas.harvard.edu
% Jeffrey Moffitt 
% lmoffitt@mcb.harvard.edu
% September 6, 2014
%--------------------------------------------------------------------------
% Based on previous functions by Alistair Boettiger
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'maxD', 'nonnegative', 8}; % Maximum distance for fiducial tracking
defaults(end+1,:) = {'fiducialFrame', 'nonnegative', 1}; % Frame to use for tracking
defaults(end+1,:) = {'fiducialWarp2Hyb1', 'boolean', false}; 
defaults(end+1,:) = {'verbose', 'boolean', true};
defaults(end+1,:) = {'printedUpdates', 'boolean', true};

defaults(end+1,:) = {'reportsToGenerate','cell',cell(0,1)};
defaults(end+1,:) = {'useSubFolderForCellReport', 'boolean', true};
defaults(end+1,:) = {'overwrite', 'boolean', true};
defaults(end+1,:) = {'figFormats', 'cell', {'png', 'fig'}};
defaults(end+1,:) = {'saveAndClose', 'boolean', false};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if (nargin < 1 || ...
        ~(strcmp(class(fiducialData), 'struct')) || ...
        (length(fiducialData) < 2) || ...
        ~ismember('mList', fields(fiducialData)) )
    
    error('matlabFunctions:AlignFiducials', 'Invalid fiducial data.');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
if parameters.printedUpdates & parameters.verbose
    display('--------------------------------------------------------------');
    display('Analyzing fiducials');
end

% -------------------------------------------------------------------------
% Handling optional variables
% -------------------------------------------------------------------------
if ~isfield(parameters,'numHybs')
    parameters.numHybs = length(fiducialData); 
end

% -------------------------------------------------------------------------
% Optional plotting
% -------------------------------------------------------------------------
reportID = find(strcmp('fiducialReport1', parameters.reportsToGenerate(:,1)));
if ~isempty(reportID)
    figHandle1 = figure('Name',['fiducialReport1_cell', num2str(fiducialData(1).cellNum)], ...
        'visible', parameters.reportsToGenerate{reportID, 2});
    clrmap = jet(length(fiducialData)); 
    r = 5; % radius for color wheels of aligned hybes
end

reportID = find(strcmp('fiducialReport2', parameters.reportsToGenerate(:,1)));
if ~isempty(reportID)
    figHandle2 = figure('Name',['fiducialReport2_cell', num2str(fiducialData(1).cellNum)], ...
        'visible', parameters.reportsToGenerate{reportID, 2});
    clrmap = jet(length(fiducialData)); 
end

% -------------------------------------------------------------------------
% Build cell array of bead x,y positions
% -------------------------------------------------------------------------
fedPos = {};
w = {};
for i=1:length(fiducialData)
    try
        mList = FilterMList(fiducialData(i).mList, 'frame', [1 2]);
        tempPos = [];
        tempPos(:,1) = mList.xc; %Kludge to force orientation of array
        tempPos(:,2) = mList.yc;
        w{i} = mList.w;
        fedPos{i} = double(tempPos);
    catch er
        warning(er.getReport);
    end
end
  
% -------------------------------------------------------------------------
% Build transforms
% -------------------------------------------------------------------------
tformNN = cell(parameters.numHybs,1);
for i=1:length(fiducialData)
    try
        fiducialData(i).warpErrors = NaN(1,5);
        if parameters.fiducialWarp2Hyb1 % warp all beads to first hybe
            [fiducialData(i).tform, tempErrors] = ...
                Warp2BestPair(fedPos{1},fedPos{i}, ...
                    'parameters', parameters, ...
                    'showPlots', false); 
        else  % warp to previous hybe
            [tformNN{i}, tempErrors] = ...
                Warp2BestPair(fedPos{max(i-1, 1)}, fedPos{i}, ...
                    'parameters', parameters, ...
                    'showPlots', false); 
                fiducialData(i).tform = maketform('composite',tformNN{1:i});
        end
        fiducialData(i).hasFiducialError = false;
        fiducialData(i).fiducialErrorMessage = [];
        inds = 1:(min(length(tempErrors), 5));
        fiducialData(i).warpErrors(inds) = tempErrors(inds);
    catch er
        disp(er.getReport);
        warning(['failed to warp data for hybe ',num2str(i),'. Using previous positions']);
        if i==1
            fiducialData(i).tform = maketform('affine',[1 0 0; 0 1 0; 0 0 1]); % don't warp if you can't warp
            tformNN{i}= maketform('affine',[1 0 0; 0 1 0; 0 0 1]); % don't warp if you can't warp
            fiducialData(i).warpErrors = NaN(1,5);
        else
            fiducialData(i).tform = fiducialData(i-1).tform; % don't warp if you can't warp
            tformNN1{i}= maketform('affine',[1 0 0; 0 1 0; 0 0 1]); % don't warp if you can't warp
            fiducialData(i).warpErrors = NaN(1,5);
        end
        fiducialData(i).hasFiducialError = true;
        fiducialData(i).fiducialErrorMessage = er;
    end
    if parameters.verbose
        disp(['Alignment error = ',num2str(fiducialData(i).warpErrors)]);
    end
    
    % -------------------------------------------------------------------------
    % Plot Fiducial Report 1 Data: Bead Positions for Each Hyb
    % -------------------------------------------------------------------------
    if exist('figHandle1')
        set(0, 'CurrentFigure', figHandle1);
        [xw,yw] = tforminv(fiducialData(i).tform,fedPos{i}(:,1),fedPos{i}(:,2));
        plot(xw,yw,'.','color',clrmap(i,:)); hold on;  
        theta = pi/length(fiducialData)*i;
        xl = reshape([xw-r*cos(theta),xw+r*cos(theta),NaN*ones(length(xw),1)]',length(xw)*3,1);
        yl = reshape([yw-r*sin(theta),yw+r*sin(theta),NaN*ones(length(yw),1)]',length(xw)*3,1);
        plot(xl,yl,'color',clrmap(i,:),'LineWidth',2); hold on;  
    end
    
    % -------------------------------------------------------------------------
    % Plot Fiducial Report 2 Data: Bead Widths for Each Hyb
    % -------------------------------------------------------------------------
    if exist('figHandle2')
        set(0, 'CurrentFigure', figHandle2);
        [xw,yw] = tforminv(fiducialData(i).tform,fedPos{i}(:,1),fedPos{i}(:,2));
        theta = pi/length(fiducialData)*i;
        r = w{i}'/80;
        xl = reshape([xw-r.*cos(theta),xw+r.*cos(theta),NaN*ones(length(xw),1)]',length(xw)*3,1);
        yl = reshape([yw-r.*sin(theta),yw+r.*sin(theta),NaN*ones(length(yw),1)]',length(xw)*3,1);
        plot(xl,yl,'color',clrmap(i,:),'LineWidth',2); hold on;  
    end
    
end
  
% -------------------------------------------------------------------------
% Finalize report figures and save
% -------------------------------------------------------------------------    
if exist('figHandle1')
    xlim([0,256]); ylim([0,256]);
    xlabel('pixels'); ylabel('pixels'); 
    
    if parameters.saveAndClose
        if parameters.useSubFolderForCellReport
            SaveFigure(figHandle1, 'overwrite', parameters.overwrite, ...
                'formats', parameters.figFormats, ...
                'subFolder', ['Cell_' num2str(fiducialData(1).cellNum)]);
        else
            SaveFigure(figHandle1, 'overwrite', parameters.overwrite, ...
                'formats', parameters.figFormats);
        end

        close(figHandle1); 
    end
end

if exist('figHandle2')
    xlim([0,256]); ylim([0,256]);
    xlabel('pixels'); ylabel('pixels'); 
    
    if parameters.saveAndClose
        if parameters.useSubFolderForCellReport
            SaveFigure(figHandle2, 'overwrite', parameters.overwrite, ...
                'formats', parameters.figFormats, ...
                'subFolder', ['Cell_' num2str(fiducialData(1).cellNum)]);
        else
            SaveFigure(figHandle2, 'overwrite', parameters.overwrite, ...
                'formats', parameters.figFormats);
        end

        close(figHandle2);
    end
end
        
