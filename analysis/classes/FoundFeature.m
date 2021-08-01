classdef FoundFeature < handle
% ------------------------------------------------------------------------
% [fFeature, parameters] = FoundFeature(..., varargin)
% This class is a container for segmented features in MERFISH data. 
%--------------------------------------------------------------------------
% Necessary Inputs
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% jeffrey.moffitt@childrens.harvard.edu
% September 21, 2017
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2018.
%--------------------------------------------------------------------------
% This class is a wrapper around a found feature, i.e. a segmented cell, and
% contains all the basic functionality associated with features

% -------------------------------------------------------------------------
% Define properties
% -------------------------------------------------------------------------
properties
    name = ''                   % A string name for this feature
    type = ''                   % A string that can be used to distinguish different feature types
    verbose = false             % A boolean that determines whether or not the class displays progress
end

properties (SetAccess=protected)
    % Class definition properties
    version = '1.0'
    
    % Properties associated with the image
    image_size      % The size of the original image (WxH)
    pixel_size      % The size of the pixel 
    
    % Properties to allow identification of feature
    uID             % A unique ID to this feature
    joinedUIDs      % A cell array of all unique IDs joined to create this cell (if any)
    fovID           % The ID of fov in which the cell appears
    feature_label   % The element in the label matrix corresponding to this feature
    
    % Properties of the boundaries in fov coordinates
    num_zPos        % The number of z slices
    boundaries      % A cell array of boundaries: A multidimensional cell if multiple fov have been grouped together
    
    % Properties of the boundaries in real coordinates
    abs_zPos            % An array of z positions
    abs_boundaries      % A cell array of boundaries
    
    % Meta data associated with the features
    volume = 0              % The number of voxels in the feature
    abs_volume = 0          % The total sample volume in the feature
    boundary_area           % An array of the number of pixels in each z stack
    abs_boundary_area       % An array of the area in each z stack
    is_broken = false       % Is the boundary broken?
    num_joined_features = 0 % The number of features joined together
    
    % Useful features for quickly indexing the feature
    feature_id = -1
    
    % Links to features contained or associated with this feature
    children        % References to features that this feature 'owns'
    parents         % References to features that own this feature
    
    % Metadata associated with the cell
    metaData        % Reserved for misc meta data associated with the feature
end

% -------------------------------------------------------------------------
% Define methods
% -------------------------------------------------------------------------
methods
    
    % -------------------------------------------------------------------------
    % Define constructor
    % -------------------------------------------------------------------------
    function obj = FoundFeature(labelMat, fovID, fovCenterPos, pixelSize, stageOrientation, boundingBox, zPos, featureLabel, varargin)
        % This class allows the construction of a defined feature based on
        % the pixels that are true within a 2D or 3D label matrix image
        % foundFeature = FoundFeature(labelMat, fovID, fovCenterPos,
        % pixelSize, stageOrientation, boundingBox, zPos, featureLabel)
        
        % -------------------------------------------------------------------------
        % Handle optional arguments
        % -------------------------------------------------------------------------
        % Define defaults
        defaults = cell(0,3);

        % Parameters for parsing file names
        defaults(end+1,:) = {'verbose', 'boolean', false};      % Display progress of construction
        defaults(end+1,:) = {'name', 'string', ''};             % A name for the OTTable instance
        
        % Parse varaible arguments with defaults
        parameters = ParseVariableArguments(varargin, defaults, mfilename);
       
        % Transfer to object
        defaultFields = fields(parameters);
        for f=1:length(defaultFields)
            obj.(defaultFields{f}) = parameters.(defaultFields{f});
        end
        
        % -------------------------------------------------------------------------
        % Handle empty class request
        % -------------------------------------------------------------------------
        if nargin < 1
            return;
        end

        % -------------------------------------------------------------------------
        % Transfer common information
        % -------------------------------------------------------------------------
        obj.fovID = fovID;
        obj.abs_zPos = zPos;
        obj.feature_label = featureLabel;
        obj.image_size = size(labelMat);
        obj.image_size = obj.image_size(1:2);
        obj.pixel_size = pixelSize;
        obj.num_zPos = size(labelMat,3);

        % -------------------------------------------------------------------------
        % Prepare the object
        % -------------------------------------------------------------------------
        % Assign a uID to this found feature
        obj.uID = char(java.util.UUID.randomUUID);

        % Prepare boundaries
        obj.boundaries = repmat({zeros([0,2])}, [2 obj.num_zPos]); % Two sets of boundaries to allow for joined features
        obj.abs_boundaries = repmat({zeros([0,2])}, [1 obj.num_zPos]);
        
        % -------------------------------------------------------------------------
        % Parse label matrix and build unique features
        % -------------------------------------------------------------------------
        % Loop through each z stack
        for z=1:obj.num_zPos

            % Find boundaries without holes
            boundaries = bwboundaries(labelMat(:,:,z), 'noholes');

            % Determine if this boundary is empty
            if isempty(boundaries)
                continue;
            end

            % Concatenate boundaries for multiple objects (i.e. connected in a different z plane)
            lBoundary = boundaries{1};
            for l=2:length(boundaries)
                lBoundary = cat(1, lBoundary, nan(1,2), boundaries{l});
            end

            % Flip X/Y coordinates to match with barcodes
            lBoundary = lBoundary(:, [2 1]);     

            % Add boundary to boundary list
            obj.boundaries{1,z} = lBoundary;

            % Convert boundary to absolute scale                
            abs_lBoundary = lBoundary - ...                    % Transform coordinate system to middle of image
                repmat([size(labelMat,1) size(labelMat,2)]/2, [size(lBoundary,1) 1]); 
            abs_lBoundary = abs_lBoundary.*...                 % Convert pixels to microns
                repmat(pixelSize/1000*stageOrientation, [size(abs_lBoundary,1),1]);

            % Apply crop (bounding box is in microns)
            inInds = true(1, size(abs_lBoundary,1));
            inInds(~isnan(abs_lBoundary(:,1))) = abs_lBoundary(~isnan(abs_lBoundary(:,1)),1) >= boundingBox(1) & ...
                abs_lBoundary(~isnan(abs_lBoundary(:,1)),1) <= (boundingBox(1) + boundingBox(3)) & ...
                abs_lBoundary(~isnan(abs_lBoundary(:,1)),2) >= boundingBox(2) & ...
                abs_lBoundary(~isnan(abs_lBoundary(:,1)),2) <= (boundingBox(2) + boundingBox(4));

            % Check to see if the only remaining true in ind is due to a
            % nan flag
            if sum(inInds) == sum(isnan(abs_lBoundary(:,1)))
                inInds = false(1, size(abs_lBoundary,1));
            end
            
            % Determine start and stop of line segments
            startInds = find(inInds & ~circshift(inInds,-1))+1; % The index of the start of the break
            numBreaks = length(startInds);

            % Handle case that the boundary has no breaks
            if isempty(startInds)
                startInds = 0;
            end

            % Shift the coordinate system to the location of the fov
            abs_lBoundary = abs_lBoundary + ...
                repmat(fovCenterPos, [size(abs_lBoundary,1) 1]);

            % Skip construction of boundaries if nothing survived the crop
            if all(~inInds) || numBreaks > 1 % Discard anything that has been cropped away
                obj.is_broken = true;
                continue;
            end

            % Shift the absolute boundaries and index out the remaining
            % boundaries
            abs_lBoundary = circshift(abs_lBoundary, -startInds(1)+1, 1);
            abs_lBoundary = abs_lBoundary((sum(~inInds)+1:end),:);
            
            % Add absolute/cropped boundary to boundary list
            obj.abs_boundaries{z} = abs_lBoundary;
            
            % Mark the feature as broken or not
            obj.is_broken = obj.is_broken | numBreaks > 0;

        end
        
        % -------------------------------------------------------------------------
        % Calculate the metadata associated with the feature
        % -------------------------------------------------------------------------
        obj.CalculateProperties();
        
    end
    
    % -------------------------------------------------------------------------
    % Calculate Feature Centroid
    % -------------------------------------------------------------------------
    function centroid = CalculateCentroid(obj)
        % This method calculates the centroid of the feature boundary in
        % absolute coordinates
        %
        % centroid = obj.CalculateCentroid(obj)
        
        % Calculate the x,y centroid for each boundary
        muCentroid = zeros(2, obj.num_zPos);
        for z=1:obj.num_zPos
            lBoundary = obj.abs_boundaries{z};
            muCentroid(:,z) = nanmean(lBoundary,1);
        end
        
        % Weight the x/y centroids by relative boundary area
        centroid(1:2) = nansum(muCentroid.*repmat(obj.abs_boundary_area, [2 1]),2)/nansum(obj.abs_boundary_area);
        
        % Weight the z centroid by the area of each boundary
        centroid(3) = nansum(obj.abs_zPos.*obj.abs_boundary_area)/nansum(obj.abs_boundary_area);
        
    end
    
    % -------------------------------------------------------------------------
    % Plot found feature in absolute coordinates on provide axis handle
    % -------------------------------------------------------------------------
    function lineHandles = Plot(obj, zInd, axisHandle)
        % Plot the found feature in the provided zInd into the provided
        % axes handle
        %
        % lineHandles = obj.Plot(zInd, axisHandle);
        
        % Create a figure and axis if not provided
        if ~exist('axisHandle', 'var')
            figHandle = figure();
            axisHandle = axis;
        end

        % Define the z ind if not provided
        if ~exist('zInd', 'var')
            zInd = 1;
        end
        
        % If multiple objects are provided (i.e. this method called from a
        % class array), then loop over then and plot
        lineHandles = [];
        for i=1:length(obj)
            localBoundaries = obj(i).abs_boundaries{zInd};
            if ~isempty(localBoundaries)
                lineHandles(end+1) = plot(localBoundaries(:,1),localBoundaries(:,2)); hold on;
            end
        end
        
    end
    
    % -------------------------------------------------------------------------
    % Calculate morphological properties of the feature
    % -------------------------------------------------------------------------
    function [eccentricity, hullRatio, numRegions] = CalculateMorphology(obj)
        % This method calculates morphological properties of the feature
        %
        % [eccentricity, hullRatio, numRegions] = obj.CalculateMorphology(obj)
        
        % Allocate memory
        eccentricity = zeros(1, obj.num_zPos);
        hullRatio = zeros(1, obj.num_zPos);
        numRegions = zeros(1, obj.num_zPos);
        
        % Loop over all z plane
        for z=1:obj.num_zPos
            % Extract boundary
            lBoundary = obj.abs_boundaries{z};
            
            % Break if is empty
            if isempty(lBoundary)
                continue;
            end
            
            % Convert to pixel values
            lBoundary = round(lBoundary/(obj.pixel_size/1000));
            
            % Remove nan
            isBadIndex = isnan(lBoundary(:,1));
            lBoundary = lBoundary(~isBadIndex,:);
            
            % Calculate the values to convert to a binary image
            edges{1} = unique(lBoundary(:,1));
            edges{2} = unique(lBoundary(:,2));
            
            % Confirm a valid boundary
            numUniquePixels = cellfun(@length, edges);
            if any(numUniquePixels <= 1)
                continue;
            end
            
            % Compute the 2D image of boundary
            image = hist3(lBoundary, 'Edges', edges)>0;
            
            % Fill the holes in the image
            image = imfill(image, 'holes');
            
            % Compute the properties
            props = regionprops(image, 'Eccentricity', 'Area', 'ConvexArea');
            
            % Keep the eccentricity of the largest region
            eccentricity(z) = props(1).Eccentricity;

            % Keep the ratio of the area to the convex area for the largest
            % region
            hullRatio(z) = props(1).Area/props(1).ConvexArea;
            
            % Record the number of regions
            numRegions(z) = length(props);
                        
        end
                
    end

    % -------------------------------------------------------------------------
    % Calculate the distance to the feature
    % -------------------------------------------------------------------------
    function dist = DistanceToFeature(obj, position)
        % This method calculates the distance between a position and the
        % boundaries of the local feature (in absolute coordinates)
        %
        % dist = obj.DistanceToFeature(position);
        
        % Check necessary input
        if nargin < 2 || ~all(size(position) == [1 3])
            error('matlabFunctions:invalidArguments', 'A 1 x 3 position vector must be provided.');
        end
        
        % Determine the proper z index
        zInd = position(1,3) >= [obj.abs_zPos] & position(1,3) < [obj.abs_zPos(2:end) Inf];

        % Calculate distance to boundary
        dist = min( sqrt(...
            sum(...
            (obj.abs_boundaries{zInd} - repmat(position(1,1:2),[size(obj.abs_boundaries{zInd},1) 1])).^2 )));
        
    end
    
    % -------------------------------------------------------------------------
    % Determine if a point is within the feature
    % -------------------------------------------------------------------------
    function inFeature = IsInFeature(obj, position)
        % This method returns a boolean specifying if the point or position
        % is within the feature. 
        %
        % inFeature = obj.IsInFeature(position);
        
        % Check necessary input
        if nargin < 2 || ~all(size(position) == [1 3])
            error('matlabFunctions:invalidArguments', 'A 1 x 3 position vector must be provided.');
        end
        
        % Determine the proper z index
        zInd = position(1,3) >= [obj.abs_zPos] & position(1,3) < [obj.abs_zPos(2:end) Inf];

        % Calculate distance to boundary
        inFeature = inpolygon(position(1,1), position(1,2), ...
            obj.abs_boundaries{zInd}(:,1), obj.abs_boundaries{zInd}(:,2));
    end
    % -------------------------------------------------------------------------
    % Assign feature id
    % -------------------------------------------------------------------------
    function AssignFeatureID(obj, featureID)
        % Assign a feature id to the object, must be a positive integer
        %
        % obj.AssignFeatureID(featureID)
        
        % Check required input
        if nargin < 2 || ~(featureID > 0) || ~(round(featureID) == featureID)
            error('matlabFunctions:invalidArguments', 'A positive integer feature ID must be provided.');
        end
        
        % Assign the feature id
        obj.feature_id = featureID;
    end
            
    % -------------------------------------------------------------------------
    % Check overlap between this feature and another feature
    % -------------------------------------------------------------------------
    function doesOverlap = DoesFeatureOverlap(obj, fFeature)
        % This method checks to see if any portion of the provided fFeature
        % overlaps with any portion of the current feature. This comparison
        % is done using boundaries in the absolute (real world) coordinates
        %
        % doesOverlap = obj.DoesFeatureOverlap(fFeature)
        
        % Check to confirm that the provided fFeature is a a foundFeature
        if nargin < 1 || ~isa(fFeature, 'FoundFeature')
            error('matlabFunctions:invalidArgument', 'A valid foundFeature object must be provided');
        end
        
        % Check to confirm that the zPositions are equivalent between these
        % two features
        if ~all(obj.abs_zPos == fFeature.abs_zPos)
            warning('matlabFunctions:differentCoordinates', 'The two FoundFeature objects do not share the same z coordinate system.');
        end
        
        % Initialize doesOverlap
        doesOverlap = false;

        % Loop through the z positions
        for z=1:length(obj.abs_zPos)
            % Identify corresponding z postion in the other feature
            zInd = find(fFeature.abs_zPos == obj.abs_zPos(z));
            
            % Handle the case that this position does not exist in both
            % features
            if isempty(zInd)
                continue;
            end
            
            % Extract boundaries for these z planes
            absBoundary1 = obj.abs_boundaries{z};
            absBoundary2 = fFeature.abs_boundaries{zInd};
            
            % Handle the position that either boundary is empty
            if isempty(absBoundary1) || isempty(absBoundary2)
                continue;
            end
            
            % Determine the points that are within or overlapping
            [in, on] = inpolygon(absBoundary2(:,1), absBoundary2(:,2), ...
                absBoundary1(:,1), absBoundary1(:,2));
            
            % Find if any points are in (but not shared between the
            % boundaries)
            doesOverlap = any(in(~on)); 
        
            % Exit as soon as any overlap is found
            if doesOverlap
                break;
            end
        end

    end
    
    % -------------------------------------------------------------------------
    % Return the primary fov associated with this feature
    % -------------------------------------------------------------------------
    function fovID = ReturnPrimaryFovID(obj)
        % Return the primary fov id associated with this feature
        % 
        % obj.ReturnPrimaryFovID()
        %
        % Also works on arrays
        % [fovIDs] = objectArray.ReturnPrimaryFovID();
        
        % Return the first fov for each object
        fovID = zeros(1, length(obj));
        for f=1:length(obj)
            fovID(f) = obj(f).fovID(1);
        end
    end
    
    % -------------------------------------------------------------------------
    % Calculate area/volume properties of the feature
    % -------------------------------------------------------------------------
    function CalculateProperties(obj)
        % Calculate the area, the volume, and the bounding boxes for both
        % pixel and absolute coordinate values
        %
        % obj.CalculateProperties();
        
        % Prepare the object
        obj.boundary_area = zeros(1, obj.num_zPos);
        obj.abs_boundary_area = zeros(1, obj.num_zPos);
        
        % Determine the thickness of the optical slice (if appropriate)
        if obj.num_zPos > 1
            zSliceThickness = mean(diff(obj.abs_zPos));
        else
            zSliceThickness = 1;
        end

        % Loop over all z positions
        for z=1:obj.num_zPos
            % Loop over all fov associated with each pixel coordinate
            % system boundary
            for f=1:length(obj.fovID)
                % Extract pixel coordinate boundary
                lBoundary = obj.boundaries{f,z};
                
                % Handle case of empty boundary
                if isempty(lBoundary)
                    continue;
                end

                % Compute metadata on this boundary (in pixel coordinates)
                obj.boundary_area(z) = obj.boundary_area(z) + polyarea(lBoundary(~isnan(lBoundary(:,1)),1), lBoundary(~isnan(lBoundary(:,1)),2));
                obj.volume = obj.volume + obj.boundary_area(z);
            end
            
            % Extract absolute coordinate boundary
            abs_lBoundary = obj.abs_boundaries{z};
            
            % Handle case of empty boundary
            if isempty(abs_lBoundary)
                continue;
            end
            
            % Determine meta data for absolute boundary
            obj.abs_boundary_area(z) = polyarea(abs_lBoundary(~isnan(abs_lBoundary(:,1)),1), abs_lBoundary(~isnan(abs_lBoundary(:,1)),2));
            obj.abs_volume = obj.abs_volume + obj.abs_boundary_area(z)*zSliceThickness;
        end
    end

    % -------------------------------------------------------------------------
    % Implement/overload a copy command to make a deep copy of the object
    % -------------------------------------------------------------------------
    function cpObj = copy(obj)
        % Make a deep copy of the object (ignoring handles)
        
        % Make a new found feature object
        cpObj = FoundFeature();
        
        % Find all properties associated with the old object
        fields = fieldnames(obj);
        
        % Remove the uID
        fields = setdiff(fields, 'uID');
        
        % Transfer fields
        for f=1:length(fields)
            cpObj.(fields{f}) = obj.(fields{f});
        end
        
        % Create new uID
        cpObj.uID = char(java.util.UUID.randomUUID);
        
    end
    % -------------------------------------------------------------------------
    % Determine the penalty associated with joining two found features
    % -------------------------------------------------------------------------
    function penalty = CalculateJoinPenalty(obj, featureToJoin)
        % Calculate the penalty associated with joining two feature objects
        % or of joining an object with itself. The penalty is defined as
        % the average distance between the broken ends of all of the
        % corresponding boundaries
                
        % Examine required input
        if ~exist('featureToJoin', 'var')
            joiningSelf = true;
        elseif ~isa(featureToJoin, 'FoundFeature')
            error('matlabFunctions:invalidArguments', 'An instance of FoundFeature must be provided');
        else
            joiningSelf = false;
        end
        
        % Examine the features to see if they can be joined
        if ~joiningSelf
            % Confirm that the features are broken
            if ~(obj.is_broken && featureToJoin.is_broken)
                error('matlabFunctions:invalidArguments', 'Found Features must have a broken boundary to be capable of being joined');
            end
            
            % Check that the two objects have the same z coordinate system
            if ~all(obj.abs_zPos == featureToJoin.abs_zPos)
                error('matlabFunctions:invalidArguments', 'To be joined, the two FoundFeature objects must have the same z coordinate system');
            end
        else
            if ~obj.is_broken
                error('matlabFunctions:invalidArguments', 'Found Features must have a broken boundary to be capable of being joined');
            end
        end

        % Handle the case that the penalty is being calculated for the
        % object itself
        if joiningSelf
            penalty = 0;
            numZSlices = 0;
            for z=1:obj.num_zPos
                % Extract absolute coordinate system boundaries
                abs_boundary = obj.abs_boundaries{z};
                
                % Handle empty boundaries for this slice
                if isempty(abs_boundary)
                    continue;
                end
                
                % Compute start to end distance
                DES = pdist2(abs_boundary(end,:), abs_boundary(1,:));
                
                % Accumulate penalty
                penalty = DES + penalty;
                
                % Increment the number of z slices
                numZSlices = numZSlices + 1;
            end
            
            % Normalize by the number of z slices
            penalty = penalty/numZSlices;
            
            % Return the penalty value
            return;
        end
        
        % Continue with the case that there are two features to compare
        penalty = 0;
        numZSlices = 0;
        for z=1:obj.num_zPos
            
            % Extract boundaries
            boundary1 = obj.abs_boundaries{z};
            boundary2 = featureToJoin.abs_boundaries{z};

            % If either is empty, skip this feature
            if isempty(boundary1) || isempty(boundary2)
                continue;
            end
            
            % Compute the distance 
            DES = pdist2(boundary1(end,:), boundary2(1,:)) + ...
                pdist2(boundary1(1,:), boundary2(end,:));
            DEE = pdist2(boundary1(end,:), boundary2(end,:)) + ...
                pdist2(boundary1(1,:), boundary2(1,:));
        
            % Determine the minimum of the two
            minD = min([DES DEE]);
            
            % Accumulate penalty
            penalty = penalty + minD;
            
            % Increment number of z slices
            numZSlices = numZSlices + 1;
        end
        
        % Normalize the penalty
        penalty = penalty/numZSlices;
        
        % Handle the case that neither boundary set existed in the same z
        % plane
        if numZSlices == 0
            penalty = inf;
        end
    end

    % -------------------------------------------------------------------------
    % Join two FoundFeature objects
    % -------------------------------------------------------------------------
    function joinedFeature = JoinFeature(obj, featureToJoin)
        % Join two found feature objects (or one found feature object
        % with itself
        %
        % joinedFeature = obj.JoinFeature(); % Join the feature with itself
        % joinedFeature = obj.JoinFeature(featureToJoin); % Join the
        % feature with another feature
        
        % Examine required input
        if ~exist('featureToJoin', 'var')
            joiningSelf = true;
        elseif ~isa(featureToJoin, 'FoundFeature')
            error('matlabFunctions:invalidArguments', 'An instance of FoundFeature must be provided');
        else
            joiningSelf = false;
        end
        
        % Examine the features to see if they can be joined
        if ~joiningSelf
            if ~(obj.is_broken && featureToJoin.is_broken)
                error('matlabFunctions:invalidArguments', 'Found Features must have a broken boundary to be capable of being joined');
            end
            % Check that the two objects have the same z coordinate system
            if ~all(obj.abs_zPos == featureToJoin.abs_zPos)
                error('matlabFunctions:invalidArguments', 'To be joined, the two FoundFeature objects must have the same z coordinate system');
            end

        else
            if ~obj.is_broken
                error('matlabFunctions:invalidArguments', 'Found Features must have a broken boundary to be capable of being joined');
            end
        end
        
        % Make copy of the current object
        joinedFeature = copy(obj);
        
        % Handle the case that we are joining this feature
        if joiningSelf
            % Transfer the uID of the previous feature
            joinedFeature.joinedUIDs{1} = obj.uID;
            
            % Loop over z positions
            for z=1:obj.num_zPos

                % Extract boundary
                boundary1 = joinedFeature.abs_boundaries{z};
                
                % Escape if the boundary is empty
                if isempty(boundary1)
                    continue;
                end
                
                % Fill edges
                endDist = sqrt(sum( (boundary1(end,:) - boundary1(1,:)).^2,2));
                numPoints = round(endDist/(obj.pixel_size/1000));
                boundary1ToBoundary1 = cat(2, ...
                    linspace(boundary1(end,1), boundary1(1,1), numPoints+2)', ... % X positions
                    linspace(boundary1(end,2), boundary1(1,2), numPoints+2)');    % Y positions
                boundary1ToBoundary1 = boundary1ToBoundary1(2:(end-1),:); % Trim off replicate points

                joinedFeature.abs_boundaries{z} = cat(1, boundary1, boundary1ToBoundary1);
            end
            
        else % Handle the case that we are joining two features
            
            % Transfer the uID of the previous feature
            joinedFeature.joinedUIDs{1} = obj.uID;
            joinedFeature.joinedUIDs{2} = featureToJoin.uID;
            
            % Remove pixel boundaries and fov information from joined
            % feature
            joinedFeature.boundaries = repmat({zeros(0,2)}, [2,obj.num_zPos]);
            joinedFeature.fovID = [obj.fovID featureToJoin.fovID];
            
            % Loop over z positions
            for z=1:obj.num_zPos
                % Extract boundaries
                boundary1 = obj.abs_boundaries{z};
                boundary2 = featureToJoin.abs_boundaries{z};
        
                % Exit if either is empty
                if isempty(boundary1) || isempty(boundary2)
                    continue;
                end
                
                % Compute the distance 
                DES = pdist2(boundary1(end,:), boundary2(1,:)) + ...
                    pdist2(boundary1(1,:), boundary2(end,:));
                DEE = pdist2(boundary1(end,:), boundary2(end,:)) + ...
                    pdist2(boundary1(1,:), boundary2(1,:));
                
                % If the end better matches to the end, then one boundary
                % needs to be inverted
                if DEE < DES
                    boundary2 = flipud(boundary2);
                end
                
                % Fill gaps: boundary 2 to boundary 1
                endDist = sqrt(sum( (boundary2(end,:) - boundary1(1,:)).^2,2));
                numPoints = round(endDist/(obj.pixel_size/1000));
                boundary2ToBoundary1 = cat(2, ...
                    linspace(boundary2(end,1), boundary1(1,1), numPoints+2)', ... % X positions
                    linspace(boundary2(end,2), boundary1(1,2), numPoints+2)');    % Y positions
                boundary2ToBoundary1 = boundary2ToBoundary1(2:(end-1),:); % Trim off replicate points
                
                % Fill gaps: boundary 1 to boundary 2
                endDist = sqrt(sum( (boundary1(end,:) - boundary2(1,:)).^2,2));
                numPoints = round(endDist/(obj.pixel_size/1000));
                boundary1ToBoundary2 = cat(2, ...
                    linspace(boundary1(end,1), boundary2(1,1), numPoints+2)', ... % X positions
                    linspace(boundary1(end,2), boundary2(1,2), numPoints+2)');    % Y positions
                boundary1ToBoundary2 = boundary1ToBoundary2(2:(end-1),:); % Trim off replicate points

                % Update the absolute boundary in the joined object
                joinedFeature.abs_boundaries{z} = cat(1, boundary1, boundary1ToBoundary2, ...
                    boundary2, boundary2ToBoundary1);
                
                % Update the pixel coordinate boundaries 
                joinedFeature.boundaries{1,z} = obj.boundaries{1,z};
                joinedFeature.boundaries{2,z} = featureToJoin.boundaries{1,z};
            end
            
            % Update the meta data associated with the joined feature
            joinedFeature.CalculateProperties();
        end
        
        % Update the number of joined features
        joinedFeature.num_joined_features = length(joinedFeature.fovID);
        joinedFeature.is_broken = false;

    end
    
    % -------------------------------------------------------------------------
    % Determine if this feature falls within a fov
    % -------------------------------------------------------------------------
    function isInFov = InFov(obj, fovIDs)
        % Determine if the feature falls within any of the specified fovIDs
        %
        % isInFov = obj.InFov(fovIDs)
        
        % Check required input
        if ~exist('fovIDs', 'var')
            error('matlabFunctions:invalidArguments', 'A set of fov IDs must be provided');
        end
        
        % Determine if the fovIDs for this feature are within the specified
        % fovIDs
        isInFov = any(ismember(obj.fovID, fovIDs));
        
    end
    
    % -------------------------------------------------------------------------
    % Dilate an absolute boundary
    % -------------------------------------------------------------------------
    function dBoundary = DilateBoundary(obj, zIndex, dilationSize)
        % Produce a dilated boundary useful for determining if an RNA is
        % within a cell or not
        %
        % boundary = obj.DilateBoundary(zIndex, dilationSize)
        
        % Extract the absolute boundary
        oBoundary = obj.abs_boundaries{zIndex};
        
        % Compute the arc length
        s = sqrt(sum((oBoundary - circshift(oBoundary,-1)).^2,2));

        % Compute the derivatives
        dx = gradient(oBoundary(:,1))./s;
        ddx = gradient(dx)./s;
        dy = gradient(oBoundary(:,2))./s;
        ddy = gradient(dy)./s;

        % Compute the curvature
        num = dx .* ddy - ddx .* dy;
        denom = dx .* dx + dy .* dy;
        denom = sqrt(denom);
        denom = denom .* denom .* denom;
        kappa = num ./ denom;
        kappa(denom < 0) = NaN;

        % Handle the case of differential directions to the curves
        signValues = [1 -1];
        
        % Use the area to determine if we have selected the correct sign
        area = polyarea(oBoundary(:,1), oBoundary(:,2));
        
        % Loop over sign values
        for i=1:length(signValues)
            % Create the normal vector
            normal = cat(2, ddx, ddy);
            normal = normal./repmat(sqrt(sum(normal.^2, 2)), [1 2]);
            normal = normal.*repmat(sign(kappa), [1 2]); % Position it correctly

            % Handle when the vector is not defined
            noVectorInds = all(normal==0,2);
            normal(noVectorInds,:) = nan;

            % Fill in all nan values with linear interpolation
            normal = fillmissing(normal, 'linear');

            % Provide a slight filter
            normal = 1/3*(circshift(normal,-1) + normal + circshift(normal,1));

            % Renormalize
            normal = normal./repmat(sqrt(sum(normal.^2, 2)), [1 2]);
            
            % Provide the new boundary
            dBoundary = oBoundary + dilationSize*signValues(i)*normal;
            
            % Calculate new area
            newArea = polyarea(dBoundary(:,1), dBoundary(:,2));
            
            % If it is larger the sign is correct, break
            if newArea > area
                break;
            end
        end
        
    end
    
    % -------------------------------------------------------------------------
    % Generate a pixel mask for the specified z indices
    % -------------------------------------------------------------------------
    function mask = GeneratePixelMask(obj, fovID, zIndices)
        % Produce a binary mask (2D or 3D) for the boundaries in the given
        % fov and with the specified zIndices
        %
        % mask = obj.GeneratePixelMask(fovID, zIndices)
        
        % Confirm that the provided values are valid
        if ~ismember(fovID, obj.fovID)
            error('matlabFunctions:invalidArgument', 'The requested fov does not contain the desired feature');
        end
        if ~all(ismember(zIndices, 1:obj.num_zPos))
            error('matlabFunctions:invalidArgument', 'The requested z indices are not present in this feature');
        end
        
        % Allocate mask memory
        mask = false(obj.image_size(1), obj.image_size(2), length(zIndices));
        
        % Determine the local fovID
        localFovID = find(obj.fovID == fovID);
        
        % Build this mask for the specified z indices
        for z=1:length(zIndices)
            % Build the outline mask
            localMask = false(obj.image_size);
            
            % Extract boundary
            lBoundary = obj.boundaries{localFovID, zIndices(z)};
            
            % Remove nan values
            lBoundary = lBoundary(~isnan(lBoundary(:,1)),:);

            % Handle empty case
            if isempty(lBoundary)
                continue;
            end
            
            % Add boundaries to mask
            localMask(sub2ind(obj.image_size, lBoundary(:,2), lBoundary(:,1))) = true;
            
            % Build the full mask by filling in the holes in this image
            mask(:,:,z) = imfill(localMask, 'holes');
        end
    end

    % -------------------------------------------------------------------------
    % Export the found feature to a table object 
    % -------------------------------------------------------------------------
    function oTable = Feature2Table(obj)
	% Output all properties of a the feature as a table. This output functionality will eventually allow features to be saved/read from a csv format.
	% oTable = obj.Feature2Table();
    
        % Determine the number of zPos
        numZPos = length(obj(1).abs_zPos);

        % Create the cell arrays for these boundaries
        boundariesX = repmat({''}, length(obj), numZPos);
        boundariesY = repmat({''}, length(obj), numZPos);

        % Fill these items
        uIDs = {obj.uID}';
        featureIDs = cat(1,obj.feature_id);
        isBrokenFlags = cat(1, obj.is_broken);
        numJoinedFeatures = cat(1, obj.num_joined_features);
        absVolumes = cat(1, obj.abs_volume);

        fovIDs = obj.ReturnPrimaryFovID()';

        % Fill the boundary cells
        for z=1:numZPos
            for o=1:length(obj)
                localBoundary = obj(o).abs_boundaries{z};
                boundaryX = cell(1, 2*size(localBoundary,1));
%                 boundaryX(1:2:end) = arrayfun(@num2str, localBoundary(:,1), 'UniformOutput', false);
                boundaryX(1:2:end) = cellstr(num2str(localBoundary(:,1)));
                boundaryX(2:2:end) = repmat({';'}, [1 size(localBoundary,1)]);
                boundariesX{o,z} = cat(2,boundaryX{:});
                boundaryY = cell(1, 2*size(localBoundary,1));
%                 boundaryY(1:2:end) = arrayfun(@num2str, localBoundary(:,2), 'UniformOutput', false);
                boundaryY(1:2:end) = cellstr(num2str(localBoundary(:,2)));
                boundaryY(2:2:end) = repmat({';'}, [1 size(localBoundary,1)]);
                boundariesY{o,z} = cat(2,boundaryY{:});
            end
        end
        
        % Create the table
        oTable = table();
        
        % Add metadata
        oTable.feature_uID = uIDs;
        oTable.feature_ID = featureIDs;
        oTable.fovID = fovIDs;
        oTable.is_broken = isBrokenFlags;
        oTable.num_joined_features = numJoinedFeatures;
        oTable.abs_volume = absVolumes;
        
        % Add boundaries
        for z=1:numZPos
            oTable.(['abs_x_boundary_' num2str(z)]) = boundariesX(:,z);
            oTable.(['abs_y_boundary_' num2str(z)]) = boundariesY(:,z);
        end
    end
    
    
%         % Create an empty table
%         oTable = table();
%         
%         % Add the uID
%         oTable.feature_uID = obj.uID;
%         
%         % Add the feature id
%         oTable.feature_iD = obj.feature_id;
%         
%         % Add feature metadata
%         oTable.fovID = obj.fovID(1);        % Primary fov ID
%         oTable.is_broken = obj.is_broken;
%         oTable.num_joined_features = obj.num_joined_features;
%         oTable.abs_volume = obj.abs_volume;
%         
%         % Add properties of the boundaries: 
%         zPos = obj.abs_zPos;
%         
%         for z=1:length(zPos)
%             localBoundary = obj.abs_boundaries{z};
%             boundaryX = cell(1, 2*size(localBoundary,1));
%             boundaryX(1:2:end) = arrayfun(@num2str, localBoundary(:,1), 'UniformOutput', false);
%             boundaryX(2:2:end) = repmat({';'}, [1 size(localBoundary,1)]);
%             oTable.(['abs_x_boundary_' num2str(z)]) = {cat(2, boundaryX{:})};
%             boundaryY = cell(1, 2*size(localBoundary,1));
%             boundaryY(1:2:end) = arrayfun(@num2str, localBoundary(:,2), 'UniformOutput', false);
%             boundaryY(2:2:end) = repmat({';'}, [1 size(localBoundary,1)]);
%             oTable.(['abs_y_boundary_' num2str(z)]) = {cat(2, boundaryY{:})};
%         end
        
end


end


