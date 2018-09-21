function [report, figHandles, parameters] = GenerateGeoTransformReport(tforms, residuals, varargin)
% ------------------------------------------------------------------------
% [report, figHandles, parameters] = GenerateGeoTransformReport(tforms, residuals)
% This function generates a report on the generated geometric
% transformations
%--------------------------------------------------------------------------
% Necessary Inputs: 
%   tforms -- A cell array corresponding to geometric transform objects
%       (affine2d()) 
%   residuals -- A cell array of the distance vectors between control
%       points
%   inds -- A cell array of index pairs 
%
%   All inputs can be either 1D cell arrays or 2D cell arrays, in the later
%   case they should correspond to bits x fov. 
%--------------------------------------------------------------------------
% Outputs: 
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
% 
%--------------------------------------------------------------------------
% Example Calls
%
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 21, 2017
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2018.
%--------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Default Reports to Generate
% -------------------------------------------------------------------------
defaultReports = cell(0,3);
defaultReports(end+1,:) = {'residualTransformError', 'on', [1 1 1700 400]}; %
defaultReports(end+1,:) = {'residualTransformErrorByPosition', 'on', [1 1 1700 400]}; %
defaultReports(end+1,:) = {'transformSummary', 'on', [1 1 1230 400]}; %

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'reportsToGenerate', 'cell', defaultReports};
defaults(end+1,:) = {'edges', 'array', {1:25:512, 1:25:512}};

% Generic display parameters
defaults(end+1,:) = {'saveAndClose', 'boolean', false};
defaults(end+1,:) = {'overwrite', 'boolean', true};
defaults(end+1,:) = {'formats', 'cell', {'png', 'fig'}};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 2
    error('matlabFunctions:invalidArguments', ...
        'Both a cell array of transforms and residuals must be provided.');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Determine display type
% -------------------------------------------------------------------------
if any(size(tforms)==1)
    displayMethod = '1D';
else
    displayMethod = '2D';
end


% -------------------------------------------------------------------------
% Calculate error
% -------------------------------------------------------------------------
% Determine if error properties will be calculated
if ~isempty(find(strcmp('residualTransformError', parameters.reportsToGenerate(:,1))));

    % Allocate error properties
    % Mean residuals
    muX = nan(size(tforms));
    muY = nan(size(tforms));
    % STD in residuals
    stdX = nan(size(tforms));
    stdY = nan(size(tforms));
    % Number of control points
    numCP = nan(size(tforms));

    % Calculate errors
    numElements = numel(tforms);
    for i=1:numElements
        localResiduals = residuals{i};
        if ~isempty(localResiduals)
            muX(i) = mean(localResiduals(:,1));
            muY(i) = mean(localResiduals(:,2));
            stdX(i) = std(localResiduals(:,1));
            stdY(i) = std(localResiduals(:,2));
            numCP(i) = size(localResiduals,1);
        end
    end

    % Archive in report
    report.muX = muX;
    report.muY = muY;
    report.stdX = stdX;
    report.stdY = stdY;
    report.numCP = numCP;
end

% -------------------------------------------------------------------------
% Calculate position dependence
% -------------------------------------------------------------------------
if ~isempty(find(strcmp('residualTransformErrorByPosition', parameters.reportsToGenerate(:,1))));

    % Calculate magnitude and angle for all residuals
    allResiduals = cat(1,residuals{:});
    if iscell(allResiduals) % Handle 2D cells
        allResiduals = cat(1,allResiduals{:});
    end
    errMag = sqrt(sum(allResiduals(:,1:2).*allResiduals(:,1:2),2));
    errAngle = asin( allResiduals(:,1)./errMag);

    % Define edges for binning
    edges = parameters.edges;

    % Allocate memory
    errX = zeros(length(edges{1})-1, length(edges{2})-1);
    errY = errX;
    numValues = errX;

    % Fill
    for i=1:(length(edges{1})-1)
        for j=1:(length(edges{2})-1)
            inds = allResiduals(:,3) >= edges{1}(i) & ...
                allResiduals(:,3) < edges{1}(i+1) & ...
                allResiduals(:,4) >= edges{2}(j) & ...
                allResiduals(:,4) < edges{2}(j+1);
            errX(i,j) = mean(allResiduals(inds,1));
            errY(i,j) = mean(allResiduals(inds,2));
            numValues(i,j) = sum(inds);
        end
    end

    % Archive in report
    report.errX = errX;
    report.errY = errY;
    report.numValues = numValues;
end

% -------------------------------------------------------------------------
% Calculate position dependence
% -------------------------------------------------------------------------
if ~isempty(find(strcmp('transformSummary', parameters.reportsToGenerate(:,1))))
    
    % Allocate variables
    offsetX = nan(size(tforms));
    offsetY = nan(size(tforms));
    angle = nan(size(tforms));
    scale = nan(size(tforms));
    
    % Unit vector in x direction
    unitVector(1,:) = [0 0];
    unitVector(2,:) = [1 0];
    
    % Calculate errors
    numElements = numel(tforms);
    for i=1:numElements
        if ~isempty(tforms(i))
            movedPoints = transformPointsForward(tforms(i), unitVector);
            
            offsetX(i) = movedPoints(1,1) - unitVector(1,1);
            offsetY(i) = movedPoints(1,2) - unitVector(1,2);
            
            diffVector = diff(movedPoints);
            
            scale(i) = sqrt(sum(diffVector.^2));
            
            angle(i) = acosd(sum(diffVector.*unitVector(2,:))/scale(i));
        end
    end
    
    % Archive values
    report.offsetX = offsetX;
    report.offsetY = offsetY;
    report.scale = scale;
    report.angle = angle;
end

% -------------------------------------------------------------------------
% Display Transform Error Report
% -------------------------------------------------------------------------
figHandles = [];
reportID = find(strcmp('residualTransformError', parameters.reportsToGenerate(:,1)));

if ~isempty(reportID)
    figHandles(end+1) = figure('Name', 'Geometric Transformation Error Report', ...
        'Color', 'w', 'visible', parameters.reportsToGenerate{reportID, 2},...
        'Position', parameters.reportsToGenerate{reportID, 3});
    data = {muX, muY, numCP, stdX, stdY};
    titles = {'Residual X', 'Residual Y', 'Number of Points', 'STD X', 'STD Y'};
    quantileRanges = {[.1 .9], [.1 .9], [0 1], [0 .9], [0 .9]};
    switch displayMethod
        case '1D'
            for i=1:length(data)
                subplot(2,3,i);
                plot(data{i});
                xlabel('FOV');
                ylabel(titles{i});
                localData = data{i};
                limits = quantile(localData(:), quantileRanges{i});
                if diff(limits) > 0 && ~any(isinf(limits))
                    ylim(limits);
                end
            end
        case '2D'
            for i=1:length(data)
                subplot(2,3,i);
                imagesc(data{i});
                colorbar;
                xlabel('FOV');
                ylabel('Imaging Rounds');
                title(titles{i});
                localData = data{i};
                limits = quantile(localData(:), quantileRanges{i});
                if diff(limits) > 0 && ~any(isinf(limits))
                    caxis(limits);
                end
            end
    end
    
    if parameters.saveAndClose
        SaveFigure(figHandles(end),'parameters', parameters);
        close(figHandles(end));
    end
end

% -------------------------------------------------------------------------
% Display Position Dependence of Transform Error Report
% -------------------------------------------------------------------------
reportID = find(strcmp('residualTransformErrorByPosition', parameters.reportsToGenerate(:,1)));
if ~isempty(reportID)
    figHandles(end+1) = figure('Name', 'Geometric Transformation Error Position Dependence Report', ...
        'Color', 'w', 'visible', parameters.reportsToGenerate{reportID, 2}, ...
        'Position', parameters.reportsToGenerate{reportID, 3});
    
    data = {errX, errY, numValues};
    titles = {'X', 'Y', 'Number'};
    quantileRanges = {[.1 .9], [.1 .9], [0 1]};

    for i=1:length(data)
        subplot(1,3,i);
        imagesc(edges{1}(1:(end-1))+diff(edges{1}), ...
            edges{2}(1:(end-1))+diff(edges{2}), data{i});
        colorbar;
        xlabel('Position X (Pixels)');
        ylabel('Position Y (Pixels)')
        localData = data{i};
        limits = quantile(localData(:), quantileRanges{i});
        if diff(limits) > 0 && ~any(isinf(limits))
            caxis(limits);
        end
    end
    
    if parameters.saveAndClose
        SaveFigure(figHandles(end),'parameters', parameters);
        close(figHandles(end));
    end
end

% -------------------------------------------------------------------------
% Display properties of the transforms
% -------------------------------------------------------------------------
reportID = find(strcmp('transformSummary', parameters.reportsToGenerate(:,1)));
if ~isempty(reportID)
    figHandles(end+1) = figure('Name', 'Geometric Transformation Summary', ...
        'Color', 'w', 'visible', parameters.reportsToGenerate{reportID, 2},...
        'Position', parameters.reportsToGenerate{reportID, 3});
    data = {offsetX, offsetY, scale, angle};
    titles = {'Offset X', 'Offset Y', 'Scale', 'Angle'};
    zlabels = {'Pixels', 'Pixels', 'Fraction', 'Degrees'};
    quantileRanges = {[.1 .9], [.1 .9], [.1 .9], [.1 .9]};
    switch displayMethod
        case '1D'
            for i=1:length(data)
                subplot(2,2,i);
                plot(data{i});
                xlabel('FOV');
                ylabel(titles{i});
                localData = data{i};
                limits = quantile(localData(:), quantileRanges{i});
                if diff(limits) > 0 && ~any(isinf(limits))
                    ylim(limits);
                end
                ylabel(zlabels{i});
            end
        case '2D'
            for i=1:length(data)
                subplot(2,2,i);
                imagesc(data{i});
                colorbar;
                xlabel('FOV');
                ylabel('Imaging Rounds');
                title(titles{i});
                localData = data{i};
                limits = quantile(localData(:), quantileRanges{i});
                if diff(limits) > 0 && ~any(isinf(limits))
                    caxis(limits);
                end
                zlabel(zlabels{i});
            end
    end
    
    if parameters.saveAndClose
        SaveFigure(figHandles(end),'parameters', parameters);
        close(figHandles(end));
    end
end


