function [tforms, mList, residuals, inds, parameters] = MLists2Transform(refList, mList, varargin)
% ------------------------------------------------------------------------
% [tforms, mList, residuals, inds, parameters] = MLists2Transform(refList, mList, varargin)
% This function returns a set of geometric transforms based on a set of
% control points found in each frame of both the reference mList, refList,
% and the mList to be warped.
%--------------------------------------------------------------------------
% Necessary Inputs: 
%   refList -- A molecule list with the following fields: x, y, and frame.
%   mList -- A molecule list with the following fields: x, y, and frame.
%       See ReadMasterMoleculeList for information on molecule lists.
%       Molecule lists must be in the compact form.
%--------------------------------------------------------------------------
% Outputs: 
%   tforms -- A cell array of geometric transform objects for each frame in
%       the specified molecule lists
%   mList -- The original moving mList but with xc and yc values updated to
%       reflect the transformation. Note these are updated only if the 
%       'applyTransform' flag is true. 
%   residuals -- A cell array of the residual distances between control
%       points after the transformation is applied. Each entry is a Nx4 set
%       of points where the first two columns represent the residual error
%       vectors and the last two the position of the original points.
%       Note that these are updated only if the 'applyTransform' flag is 
%       true.
%   inds -- A cell array of the points in each frame of the mList that
%       that were assigned to points in the refList
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
% 
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% jeffrey.moffitt@childrens.harvard.edu
% September 21, 2017
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2018.
%--------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);

% Parameters for the geometric transformation
defaults(end+1,:) = {'transformationType', {'nonreflectivesimilarity', 'similarity', ...
    'affine', 'projective', 'polynomial'}, 'nonreflectivesimilarity'}; % Type of geometric transformation
defaults(end+1,:) = {'polynomialOrder', 'positive', 3}; % Order of the polynomial fit

% Parameters for control point association
defaults(end+1,:) = {'controlPointMethod', {'nearestNeighbor', 'kNNDistanceHistogram'}, 'nearestNeighbor'}; % The method used to associate points in the refList and the mList
defaults(end+1,:) = {'distanceWeight', 'positive', 0.25};   % The weight applied to the STD of the distances for valid control point selection
defaults(end+1,:) = {'thetaWeight', 'positive', 0.25};      % The weight applied to the STD of the angles for valid control point selection
defaults(end+1,:) = {'histogramEdges', 'array', -128:1:128};% The bin edges for the point different histogram
defaults(end+1,:) = {'numNN', 'positive', 10};                 % The number of nearest neighbors to compute
defaults(end+1,:) = {'pairDistTolerance','positive', 1};    % The distance threshold to find paired points after crude shift (multiples of the histogramEdges step)

% Apply transformation
defaults(end+1,:) = {'applyTransform', 'boolean', true};

% How to handle molecules in different frames
defaults(end+1,:) = {'ignoreFrames', 'boolean', false}; % If false, combine all molecules in different frames
defaults(end+1,:) = {'transpose', 'boolean', false}; % Whether to transpose data in X/Y or not

% Reporting/Debugging: WILL BE REMOVED IN THE FUTURE
defaults(end+1,:) = {'debug', 'boolean', false};
defaults(end+1,:) = {'displayFraction', 'positive', 0.05};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 2 || ...
        ~isempty(setdiff({'x', 'y', 'frame'}, fields(refList))) || ...
        ~isempty(setdiff({'x', 'y', 'frame'}, fields(mList))) || ...
        length(refList)> 1 || length(mList) > 1
    error('matlabFunctions:invalidArguments', 'Two valid molecule lists must be provided.');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Define default values
% -------------------------------------------------------------------------
if parameters.ignoreFrames
    numFrames = 1;
else
    numFrames = max(refList.frame);
end
residuals = {};
inds = {};
tforms = repmat({affine2d()}, [1 numFrames]); % Default is no transformation

% -------------------------------------------------------------------------
% Loop over individual frames
% -------------------------------------------------------------------------
for f=1:numFrames
    if parameters.ignoreFrames
        if parameters.transpose
            refPoints = [refList.x' refList.y'];
            movPoints = [mList.x' mList.y'];
        else
            refPoints = [refList.x refList.y];
            movPoints = [mList.x mList.y]; % Moving points: to warp
        end
            
    else
        if parameters.transpose
            % Extract points points for each frame
            refPoints = [refList.x(refList.frame == f)' refList.y(refList.frame == f)'];
            movPoints = [mList.x(mList.frame == f)' mList.y(mList.frame == f)']; % Moving points: to warp
        else
            % Extract points points for each frame
            refPoints = [refList.x(refList.frame == f) refList.y(refList.frame == f)];
            movPoints = [mList.x(mList.frame == f) mList.y(mList.frame == f)]; % Moving points: to warp
        end
    end
    
    % Check for empty mLists
    if isempty(refPoints)
        warning('matlabFunctions:emptyMList', 'No molecules in reference list. Using default transformation.');
        residuals{f} = zeros(0,4);
        continue;
    elseif isempty(movPoints)
        warning('matlabFunctions:emptyMList', 'No molecules in moving list. Using default transformation.');
        residuals{f} = zeros(0,4);
        continue;
    end

    % Use different methods to identify control point pairs
    switch parameters.controlPointMethod
        case 'nearestNeighbor'
            % Find the nearest neighbors and distances
            [idx, D] = knnsearch(movPoints,refPoints);

            % Find angle of vector between nearest neighbors
            theta = asin( (movPoints(idx,2)-refPoints(:,2))./D );

            if ~all(D==0) 
                % Define thresholds for excluding outliers
                DBounds = median(D) + std(D)*parameters.distanceWeight*[-1 1];
                thetaBounds = median(theta)+ std(theta)*parameters.thetaWeight*[-1 1];

                % Define the control points to keep
                pointsToKeep = find(D>=DBounds(1) & D<=DBounds(2) & ...
                    theta >= thetaBounds(1) & theta <= thetaBounds(2));
            else % If the lists match exactly
                pointsToKeep = (1:length(idx))';
            end
            
        case 'kNNDistanceHistogram'
            % Calculate differences in X and Y for all neighborhoods of
            % Find the nearest neighbors and distances
            [idx, ~] = knnsearch(movPoints,refPoints, 'K', parameters.numNN);
            
            % Flatten indices
            idx1 = idx(:);
            idx2 = repmat((1:size(refPoints,1))', [parameters.numNN 1]);
            
            % Handle the case that there are less than 'K' points
            idx2 = idx2(1:length(idx1));
            
            % Define distances
            diffX = movPoints(idx1,1) - refPoints(idx2,1);
            diffY = movPoints(idx1,2) - refPoints(idx2,2);
            
            % Create data matrix
            data = [diffX diffY];
            
            % Calculate 3D histogram
            N = hist3(data, 'edges', {parameters.histogramEdges, parameters.histogramEdges});
            
            % Find peak position
            [maxValue, ind] = max(N(:));
            [xpeak,ypeak] = ind2sub(size(N), ind);
                        
            % Determine offsets
            offset = [parameters.histogramEdges(xpeak) parameters.histogramEdges(ypeak)];
            
            % Assign control points from nearest neighbors
            [idx, D] = knnsearch(movPoints - repmat(offset, [size(movPoints,1) 1]), ...
                refPoints);
            
            % Keep all points separated by less than the pixelation of the
            % crude shift
            pointsToKeep = find(D <= parameters.pairDistTolerance*mean(diff(parameters.histogramEdges)));

            % Issue warnings
            if maxValue < mean(N(:)) + 2*std(N(:))
                warning('matlabFunctions:controlPoints', ['The maximum kNN offset is <2 times the std from the mean. ' ...
                    'Control point identification may be inaccurate']);
            end
    end

    % Build transform
    try
        if strcmp(parameters.transformationType, 'polynomial')
            tforms{f} = fitgeotrans(movPoints(idx(pointsToKeep),:), refPoints(pointsToKeep, :), ...
                parameters.transformationType, ...
                parameters.polynomialOrder);
        else
            tforms{f} = fitgeotrans(movPoints(idx(pointsToKeep),:), refPoints(pointsToKeep, :), ...
                parameters.transformationType);
        end
    catch
        warning('Did not find sufficient control points. Using default transformation');
    end

    % Archive control point indices
    inds{f} = [pointsToKeep idx(pointsToKeep)]; % Indices of reference points; indices of moving points
    
    % Apply transformation
    if parameters.applyTransform
        % Move points
        movedPoints = transformPointsForward(tforms{f}, movPoints);

        % Store in xc and yc if they exist
        if isfield(mList, 'xc') && isfield(mList, 'yc');
            if parameters.ignoreFrames
                mList.xc = movedPoints(:,1);
                mList.yc = movedPoints(:,2);
            else
                mList.xc(mList.frame==f) = movedPoints(:,1);
                mList.yc(mList.frame==f) = movedPoints(:,2);
            end
        else
            warning('matlabFunctions:missingFields', 'Transformed points were not stored because the xc and yc fields were missing from the provided mList');
        end

        residuals{f} = [movedPoints(idx(pointsToKeep),:) - refPoints(pointsToKeep,:), movPoints(idx(pointsToKeep),:)];
        % Handle empty residuals
        if isempty(residuals{f})
            residuals{f} = zeros(0,4);
        end
    end
    
    % Display progress if in debug mode
    if parameters.debug
        if rand(1) < parameters.displayFraction
            figHandle = figure('Name', ['FOV ' num2str(f-1)]);
            plot(refPoints(:,1), refPoints(:,2), 'xr'); hold on;
            plot(movPoints(:,1), movPoints(:,2), 'go'); 

            movedPoints = transformPointsForward(tforms{f}, movPoints);
            plot(movedPoints(:,1), movedPoints(:,2), 'bs');
            drawnow;
        end
    end
    
end

% -------------------------------------------------------------------------
% Flatten output if requested
% -------------------------------------------------------------------------
if parameters.ignoreFrames
    if ~isempty(residuals) % Handle case of no molecules in one of the lists
        residuals = residuals{1};
    end
    tforms = tforms{1};
end
