classdef MERFISHDecoder < handle
% ------------------------------------------------------------------------
% [mdObject] = MERFISHDecoder(varargin)
% This class provides a wrapper for a MERFISH data set, the basic functions
% for data normalization, data analysis, and data visualization. 
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% jeffrey.moffitt@childrens.harvard.edu
% September 21, 2017
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2018.
%--------------------------------------------------------------------------
% This class is a wrapper around a MERFISH data set. It coordinates all 
% aspects of analysis of MERFISH data, including decoding, segmentation, and
% parsing. It also provides basic functionality to access different aspects 
% of a MERFISH data set. 

% -------------------------------------------------------------------------
% Define properties
% -------------------------------------------------------------------------
properties 
    verbose                 % Display progress in various methods
    overwrite               % Determine if existing files will be overwrite
    sliceIDs                % A cell array containing fovIDs that should be associated within individual slices
    name                    % A string that can be used to distinguish different decoder objects
    hal_version = 'hal1';   % A string that determines the version of hal info/xml files to use for image metadata
end
 
properties (SetAccess=protected)
    % Version 
    version = '0.6';        % Version number--used to identify how to load saved classes

    % Metadata associated with raw data
    rawDataPath             % The path to the raw data
    dataOrganization        % The data organization file
    dataOrganizationPath    % Path to the data organization file
    
    % Properties for parallel processing
    parallel                % A parallel.Pool object
    numPar                  % The number of parallel workers
    
    % Metadata associated with codebook
    codebookPath            % Path to the codebook
    codebook                % The codebook
    codebookHeader          % Metadata associated with the codebook
    
    % Metadata associated with normalized data
    normalizedDataPath                                  % Base path to the data set
    mDecoderPath = ['mDecoder' filesep];                % Relative path to the MERFISHDecoder instance
    fiducialDataPath = ['fiducial_data' filesep];       % Relative path to information on the fiducial warping process
    warpedDataPath = ['warped_data' filesep];           % Relative path to warped tiff stacks
    processedDataPath = ['processed_data' filesep];     % Relative path to pre-processed tiff stacks
    barcodePath = ['barcodes' filesep];                 % Relative path to barcode data
    reportPath = ['reports' filesep];                   % Relative path to various reports/perfomance metrics    
    segmentationPath = ['segmentation' filesep];        % Relative path to segmentation results
    summationPath = ['summation' filesep];              % Relative path to summation results
    mosaicPath = ['mosaics' filesep];                   % Relative path to low resolution mosaic images
    
    % Metadata associated with the MERFISH run
    numDataChannels         % Number of data channels, i.e. specific image channels
    
    numBits                 % The number of bits
    bitNames                % The names of the individual readout probes for each bit
    numBarcodes             % The number of barcodes
    
    numFov                  % The number of FOV
    fovIDs                  % The fov ids
    fovPos                  % The position of each fov in microns (Nx2)
    
    numZPos                 % The number of z stacks for each data frame
    zPos                    % The position of each z-stack in microns 

    imageSize               % The number of pixels (HxW)
    imageExt                % Extension of the raw image data
    pixelSize               % The size of the pixel
    
    cameraIDs = {''}        % The ids associated with cameras for imaging
    numCameraIDs = 1        % The number of cameras used to collect data 

    numImagingRounds        % The number of imaging rounds
    imageRoundIDs           % The ids of the different imaging rounds
    
    % Storage for parameters associated with all aspects of decoder behavior
    parameters              % Storage for all parameters
    
    % Parameters associated with decoding
    scaleFactors            % The relative scaling for each imaging round
    initScaleFactors        % The initial scale factors determined from all fov
    pixelHistograms         % The pixel histograms for the processed data
    optFovIDs               % The FOV ids used for optimization
    
    % Parameters associated with the raw data and data IO
    rawDataFiles            % List of raw data files
    fov2str                 % A function for converting fov ids to strings
    
    % Properties associated with warping
    affineTransforms        % All affine transformations
    residuals               % All residuals
    geoTransformReport      % The geometric transform report
    
    % Properties associated with sliced decoders
    originalMaxFovID = []   % The original maximum fov id number from downsampled data
        
end

% -------------------------------------------------------------------------
% Public Methods
% -------------------------------------------------------------------------
methods
    
    % -------------------------------------------------------------------------
    % Define constructor
    % -------------------------------------------------------------------------
    function obj = MERFISHDecoder(rawDataPath, normalizedDataPath, varargin)
        % This method is a wrapper for a MERFISH dataset. 
        % It provides an internal mapping of the raw data, normalized data, 
        % and methods for normalizing, decoding, and visualizing data
        %
        % Examples
        % mDecoder = MERFISHDecoder(rawDataPath) % Create a MERFISH decoder object
        % mDecoder = MERFISHDecoder() % Create an empty decoder containing all default values for parameters
        % mDecoder = MERFISHDecoder(rawDataPath, 'variableName', value) %
        %    Pass a named variable to the decoder to overwrite a default.
        %    See property summary by typing doc MERFISHDecoder
        
        % -------------------------------------------------------------------------
        % Define default values for main class parameters
        % -------------------------------------------------------------------------
        % Create defaults cell
        defaults = cell(0,3); 
        
        % Paths to various metadata
        defaults(end+1,:) = {'codebookPath', ...                % Path to the codebook
            'filePath', []};
        defaults(end+1,:) = {'dataOrganizationPath', ...        % Path to a data organization file
            'filePath', []};
        
        % Basic class functionality parameters
        defaults(end+1,:) = {'verbose', 'boolean', true};       % Verbosity of class
        defaults(end+1,:) = {'overwrite', 'boolean', false};    % Determine if some files will be overwritten
        
        % Properties of the data 
        defaults(end+1,:) = {'imageExt', ...                    % Type of the raw data
            {'dax', 'tif', 'tiff'}, 'dax'}; 
        defaults(end+1,:) = {'pixelSize', ...                   % The pixel size in nm
            'nonnegative', 109};
        defaults(end+1,:) = {'hal_version', ...                     % The microscope control software version
            {'hal1', 'hal2'}, 'hal1'};
        
        % Metadata associated with the dataset
        defaults(end+1,:) = {'sliceIDs', 'cell', {}};          % Information regarding the organization of slices
        
        % Properties associated with parallel computing
        defaults(end+1,:) = {'parallel', 'parallel', []};       % A parallel pool object
        
        % Properties associated with the internal organization of the
        % decoder data and the location of generated files
        defaults(end+1,:) = {'fiducialDataPath', 'dir', []};    % Path where fiducial data will be saved
        defaults(end+1,:) = {'warpedDataPath', 'dir', []};      % Path to location of warped tiff stacks
        defaults(end+1,:) = {'processedDataPath', 'dir', []};   % Path to the location of pre-processed tiff stacks
        defaults(end+1,:) = {'barcodePath', 'dir', []};         % Path to location of barcodes
        defaults(end+1,:) = {'mDecoderPath', 'dir', []};        % Path to the MERFISH decoder
        
        % -------------------------------------------------------------------------
        % Parse the optional 'parameters' input
        % -------------------------------------------------------------------------
        % Determine if parameters was provided, and extract this field
        inputParameterNames = varargin(1:2:end);
        inputParameterValues = varargin(2:2:end);
        parameterStructNameInd = find(ismember('parameters', inputParameterNames));
        
        % Extract these extra values
        additionalInputParameterNames = {};
        additionalInputParameterValues = {};
        if ~isempty(parameterStructNameInd)
            % Extract names of the parameters structure fields
            additionalInputParameterNames = fieldnames(inputParameterValues{parameterStructNameInd(1)});
            % Extract their values
            for f=1:length(additionalInputParameterNames)
                additionalInputParameterValues{f} = inputParameterValues{parameterStructNameInd(1)}.(additionalInputParameterNames{f});
            end
        end
        
        % Remove the parameters entry
        inputParameterNames = inputParameterNames(setdiff(1:length(inputParameterValues), parameterStructNameInd));
        inputParameterValues = inputParameterValues(setdiff(1:length(inputParameterValues), parameterStructNameInd));
        
        % Add the extra values to the provided input
        inputParameterNames(end+1:end+length(additionalInputParameterNames)) = additionalInputParameterNames;
        inputParameterValues(end+1:end+length(additionalInputParameterNames)) = additionalInputParameterValues;
            
        % -------------------------------------------------------------------------
        % Load the values
        % -------------------------------------------------------------------------

        % Identify parameters that correspond to properties of the decoder
        % object
        matchInd = find(ismember(inputParameterNames, setdiff(defaults(:,1), 'parallel')));
        combinedParameterInfo = cell(1, 2*length(matchInd));
        combinedParameterInfo(1:2:end) = inputParameterNames(matchInd);
        combinedParameterInfo(2:2:end) = inputParameterValues(matchInd);

        parameters = ParseVariableArguments(combinedParameterInfo, defaults, mfilename);
        fieldsToTransfer = fields(parameters);

        % Assign main class parameters
        for f=1:length(fieldsToTransfer)
            if ~isempty(parameters.(fieldsToTransfer{f})) % Check for empty parameters to not overwrite defaults
                obj.(fieldsToTransfer{f}) = parameters.(fieldsToTransfer{f});
            end
        end
        
        % Set parallel pool object
        SetParallel(obj, parameters.parallel); 

        % -------------------------------------------------------------------------
        % Parse variable input for all other methods
        % -------------------------------------------------------------------------
        matchInd = find(~ismember(inputParameterNames, defaults(:,1)));
        combinedParameterInfo = cell(1, 2*length(matchInd));
        combinedParameterInfo(1:2:end) = inputParameterNames(matchInd);
        combinedParameterInfo(2:2:end) = inputParameterValues(matchInd);

        obj.parameters = obj.InitializeParameters(combinedParameterInfo{:}); % Pass variable arguments to be parsed for parameters
        
        % -------------------------------------------------------------------------
        % Return empty decoder if no arguments were provided
        % -------------------------------------------------------------------------
        if nargin < 1
            return;
        end
        
        % -------------------------------------------------------------------------
        % Parse required input
        % -------------------------------------------------------------------------
        if nargin < 2 || ~exist(rawDataPath, 'dir')
            error('matlabFunctions:invalidArguments', 'Both a valid raw data path and normalized data path must be provided');
        end
        
        % -------------------------------------------------------------------------
        % Coerce final character of provided paths to filesep
        % -------------------------------------------------------------------------
        if ~strcmp(rawDataPath(end), filesep)
            rawDataPath(end+1) = filesep;
        end
        if ~strcmp(normalizedDataPath(end), filesep)
            normalizedDataPath(end+1) = filesep;
        end
        
        % -------------------------------------------------------------------------
        % Set internal paths
        % -------------------------------------------------------------------------
        obj.rawDataPath = rawDataPath;
        if obj.verbose
            PageBreak();
            display(['Generating MERFISHDecoder for data in ' obj.rawDataPath]);
        end
        
        % -------------------------------------------------------------------------
        % Check validity of normalized data path
        % -------------------------------------------------------------------------
        if ~exist(normalizedDataPath, 'dir')
            status = mkdir(normalizedDataPath);
            if ~status
                error('matlabFunctions:invalidArguments', 'Could not create normalizedDataPath');
            end
        end
        display(['Saving results in ' normalizedDataPath]);
        obj.normalizedDataPath = normalizedDataPath;
        
        % -------------------------------------------------------------------------
        % Find data organization file
        % -------------------------------------------------------------------------
        if isempty(obj.dataOrganizationPath) % Look for file if no path was provided
            placesToLook = {'', ['settings' filesep]}; % Define relative paths to check
            
            % Look for the data organization file in these paths, stop at the first one found
            for p=1:length(placesToLook) 
                if exist([obj.rawDataPath placesToLook{p} 'data_organization.csv'], 'file')
                    obj.dataOrganizationPath = [obj.rawDataPath placesToLook{p} 'data_organization.csv'];
                    break; % Stop if a valid file has been found
                end
            end
        end
        
        disp(obj.dataOrganizationPath);

        % Confirm that a valid file was found
        if ~exist(obj.dataOrganizationPath, 'file')
            error('matlabFunctions:invalidArguments', ...
                'Could not find the data organization file.');
        end

        % Copy data organization file and redefine path to copied version
        copyfile([obj.dataOrganizationPath], [obj.normalizedDataPath  'data_organization.csv']);
        obj.dataOrganizationPath = 'data_organization.csv';
        
        if obj.verbose
            display(['Copied data organization file to ' normalizedDataPath filesep 'data_organization.csv']);
        end

        % -------------------------------------------------------------------------
        % Find codebook file
        % -------------------------------------------------------------------------
        if isempty(obj.codebookPath) % Look for codebook if no path was provided
            placesToLook = {'', ['settings' filesep]}; % Define places to look
            for p=1:length(placesToLook) % Look for a valid codebook
                localBasePath = [obj.rawDataPath placesToLook{p}];
                fileStruct = dir([localBasePath '*codebook.csv']);
                if ~isempty(fileStruct) || length(fileStruct) == 1
                    obj.codebookPath = [localBasePath fileStruct.name];
                    break; % Stop when the first codebook is found
                end
            end
        end
        
        % Confirm that a valid file was found
        if ~exist(obj.codebookPath, 'file')
            warning('matlabFunctions:missingFile', 'Could not find codebook');
        end
        
        % Copy codebook and update codebook path to the copied version
        [~, codebookName, codebookExt] = fileparts(obj.codebookPath);
        copyfile([obj.codebookPath], [obj.normalizedDataPath codebookName codebookExt]);
        obj.codebookPath = [codebookName codebookExt];
        
        if obj.verbose
            display(['Copied codebook file to ' normalizedDataPath obj.codebookPath]);
        end

        % -------------------------------------------------------------------------
        % Load data organization file
        % -------------------------------------------------------------------------
        [obj.dataOrganization, metaData] = LoadDataOrganization(...
            [obj.normalizedDataPath obj.dataOrganizationPath], ...
            'verbose', obj.verbose);

        % -------------------------------------------------------------------------
        % Load codebook
        % -------------------------------------------------------------------------
        [codebook, codebookHeader] = LoadCodebook([obj.normalizedDataPath obj.codebookPath], ...
            'verbose', obj.verbose, ...
            'barcodeConvFunc', @(x)bi2de(x=='1'));
        obj.bitNames = codebookHeader.bit_names;
        obj.numBits = length(codebookHeader.bit_names);
        obj.numBarcodes = length(codebook);
        obj.codebook = codebook;
        obj.codebookHeader = codebookHeader;
                
        % -------------------------------------------------------------------------
        % Extract properties of data from data organization file
        % -------------------------------------------------------------------------
        obj.numDataChannels = metaData.numDataChannels;
        obj.numZPos = metaData.numZPos;
        obj.zPos = metaData.zPos;
        
        % -------------------------------------------------------------------------
        % Sort data organization file to match the order of the bits in the
        % codebook
        % -------------------------------------------------------------------------
        [isInCodebook, sind] = ismember({obj.dataOrganization.bitName}, obj.bitNames); % Sort data organization to the order of bit names
        
        tempDataOrg = cat(1,obj.dataOrganization(sind(isInCodebook)), ...              % Order data organization: bits in codebook, then all other data channels in order of the data org file
            obj.dataOrganization(~isInCodebook));
        
        obj.dataOrganization = tempDataOrg; % Assign temporary data organization to the object property to set the proper order
        
        % Cross check
        if sum(isInCodebook) ~= length(obj.bitNames)
            error('matlabFunctions:invalidArguments', 'Not all bits in the codebook are present in the data organization file');
        end
        
        % -------------------------------------------------------------------------
        % Organize the raw data files
        % -------------------------------------------------------------------------
        obj.rawDataFiles = obj.MapRawData(); % Map the raw data 
        
        % -------------------------------------------------------------------------
        % Load image metadata
        % -------------------------------------------------------------------------
        obj.LoadImageMetaData();
                
    end
    
    % -------------------------------------------------------------------------
    % Calculate barcode counts per feature
    % -------------------------------------------------------------------------
    function T = FoundFeaturesToCSV(obj, varargin)
        % Export a single z plane of boundaries to a csv file with
        % additional found features metadata
        %
        % obj.FoundFeaturesToCSV();
        % obj.FoundFeaturesToCSV('downSampleFactor',N); % Down sample
        % boundaries by N-fold
        % obj.FoundFeaturesToCSV('zIndex', M); Export the boundaries found
        % in the M-th z plane
        
        % Handle the variable input parameters
        defaults = cell(0,3);
        
        defaults(end+1,:) = {'downSampleFactor', ...    % The degree to which the boundaries are downsampled
            'nonnegative', 10};
        defaults(end+1,:) = {'zIndex', ...              % The index of the z plane for the downsample boundaries
            'nonnegative', 4};
        
        parameters = ParseVariableArguments(varargin, defaults, mfilename);
        
        % Display progress
        if obj.verbose
            PageBreak();
            disp(['Exporting found features']);
            disp(['   Down sampling: ' num2str(parameters.downSampleFactor)]);
            disp(['   z Index: ' num2str(parameters.zIndex)]);
            localTimer = tic;
        end
        
        % Load the found features
        foundFeatures = obj.GetFoundFeatures();
        
        % Extract metadata for the table
        featureUIDs = {foundFeatures.uID};
        featureIDs = [foundFeatures.feature_id];
        primaryFovIDs = [foundFeatures.ReturnPrimaryFovID()];
        featureVolume = [foundFeatures.abs_volume];
        
        % Handle the case of a improper requested z bondary
        if parameters.zIndex > foundFeatures(1).num_zPos
            warning(['The requested z index is not within the found features. Defaulting to the first plane']);
            parameters.zIndex = 1;
        end
        
        % Extract the downsample boundaries
        boundaryX = cell(1, length(foundFeatures));
        boundaryY = cell(1, length(foundFeatures));
        centroids = zeros(length(foundFeatures),3);
        delimiter = ';';      
        
        % Loop over the found features to extract feature-specific info
        for f=1:length(foundFeatures)
            % Extract centroid
            centroids(f,:) = foundFeatures(f).CalculateCentroid();
            % Extract Downsampled boundary
            lBoundary = foundFeatures(f).abs_boundaries{parameters.zIndex};
            if isempty(lBoundary)
                boundaryX{f} = '';
                boundaryY{f} = '';
                continue;
            end
            % Downsample the boundary
            lBoundary = lBoundary(1:parameters.downSampleFactor:end,:);
            lBoundary(end+1,:) = lBoundary(1,:); % Close the boundary for display purposes
            
            % Convert to delimited string for simple export
            boundaryString = cell(2, size(lBoundary,1)*2 -1);
            boundaryString(:,2:2:end) = repmat({delimiter}, [2 size(lBoundary,1)-1]);
            boundaryString(:,1:2:end) = arrayfun(@num2str, lBoundary', 'UniformOutput', false);
            
            boundaryX{f} = cat(2, boundaryString{1,:});
            boundaryY{f} = cat(2, boundaryString{2,:});
            
            % Display progress
            if obj.verbose && ~mod(f, 1000)
                disp(['...completed parsing ' num2str(f) ' of ' num2str(length(foundFeatures)) ' features']);
            end
        end
        
        % Create a table to which the information is going to be added
        T = table(featureUIDs', featureIDs', primaryFovIDs',featureVolume', ...
            centroids, boundaryX', boundaryY', ...
            'VariableNames', {'feature_uID', 'feature_ID', 'primary_fovID', ...
            'abs_volume', 'centroid', 'boundaryX', 'boundaryY'});
        
        % Write table
        if ~exist([obj.normalizedDataPath obj.reportPath])
            mkdir([obj.normalizedDataPath obj.reportPath]);
        end
        tablePath = [obj.normalizedDataPath obj.reportPath 'feature_metadata.csv'];
        writetable(T, tablePath);
        
        % Display progress
        if obj.verbose
            disp(['...wrote to ' tablePath]);
            disp(['...completed export in ' num2str(toc(localTimer)) ' s']);
        end

    end
    
    % -------------------------------------------------------------------------
    % Calculate barcode counts per feature
    % -------------------------------------------------------------------------
    function BarcodesToCSV(obj, varargin)
        % Load, parse, and save properties of barcodes
        %
        % obj.BarcodesToCSV();
        % obj.BarcodesToCSV('fieldsToExport',
        % {'fieldName1',....,'fieldNameN'}); % Export only the provided
        % fields
        % obj.BarcodesToCSV('parsedBarcodes', true); Export the parsed
        % barcodes
        
        % -------------------------------------------------------------------------
        % Handle defaults for varariable arguments
        % -------------------------------------------------------------------------
        % Create defaults cell
        defaults = cell(0,3); 
        
        % Paths to various metadata
        defaults(end+1,:) = {'fieldsToExport', ...                % Default barcode metadata to export
            'cell', {'barcode_id', 'fov_id', 'total_magnitude', ...
            'area', 'abs_position', 'error_bit', 'error_dir', 'feature_id', 'in_feature'}};
        defaults(end+1,:) = {'parsedBarcodes', ...                % Export parsed or not parsed barcodes?
            'boolean', true};
        
        % Parse variable input
        parameters = ParseVariableArguments(varargin, defaults, mfilename);

        % -------------------------------------------------------------------------
        % Prepare for loading barcode list files
        % -------------------------------------------------------------------------
        % Display progress
        if obj.verbose
            PageBreak();
            disp(['Loading barcodes to export metadata to csv file']);
            totalTimer = tic;
        end
        
        % Determine the requested barcode type
        if parameters.parsedBarcodes
            barcodePath = [obj.normalizedDataPath obj.barcodePath filesep 'parsed_fov' filesep];
            if ~all(obj.CheckStatus([], 'p'))
                error('matlabFunctions:missingData', 'Some fov have not yet been parsed!');
            end
        else
            barcodePath = [obj.normalizedDataPath obj.barcodePath filesep 'barcode_fov' filesep];
            if ~all(obj.CheckStatus([], 'd'))
                error('matlabFunctions:missingData', 'Some fov have not yet been decoded!');
            end
        end
        
        % Transfer some information to local variables
        fovIDs = obj.fovIDs;
        numFov = length(fovIDs);
        
        % Setup parallel processes
        spmd (obj.numPar)
            % Initialize storage for barcodes
            combinedList = [];
            
            % Loop over a set of blocks per worker
            for b=labindex:numlabs:numFov
                % Determine local fov
                localFovID = fovIDs(b);
                
                % Display progress
                if obj.verbose
                    displayStrings = {};
                    displayStrings{end+1} = PageBreak('nodisplay');
                    displayStrings{end+1} = ['Exporting metadata for fov ' obj.fov2str(localFovID)];
                    loadTimer = tic;
                end
                
                % Load barcodes
                aList = ReadBinaryFile([barcodePath 'fov_' obj.fov2str(localFovID) '_blist.bin']);
                
                % Check for empty barcodes
                if isempty(aList)
                    continue;
                end

                % Cut List
                aList = aList([aList.area] >= obj.parameters.quantification.minimumBarcodeArea & ...
                    double([aList.total_magnitude])./double([aList.area]) >= obj.parameters.quantification.minimumBarcodeBrightness);
                
                % Display progress
                if obj.verbose
                    displayStrings{end+1} = (['Loaded and cut ' num2str(length(aList)) ' barcodes']); 
                end
                
                % Check for empty barcodes
                if isempty(aList)
                    continue;
                end
                
                % Cut list to only the requested fields
                fieldsToRemove = setdiff(fieldnames(aList), parameters.fieldsToExport);
                aList = rmfield(aList, fieldsToRemove);
                
                % Append list
                combinedList = cat(1, combinedList, aList);

                % Display progress
                if obj.verbose
                    displayStrings{end+1} = (['..completed in ' num2str(toc(loadTimer)) ' s']);
                    display(char(displayStrings)); % Flush buffer
                end

            end
        end
        
        % -------------------------------------------------------------------------
        % Compile composite objects
        % -------------------------------------------------------------------------
        % Display progress
        if obj.verbose
            PageBreak();
            disp(['Combining barcode metadata from all parallel objects']);
            localTimer = tic;
        end

        % Combine composite objects
        temp = combinedList{1};
        for i=2:length(combinedList)
            temp = cat(1, temp, combinedList{i});
            combinedList{i} = []; % Destroy the extra copy of the bList to save space
        end
        combinedList = temp;

        % Display progress
        if obj.verbose
            disp(['...completed in ' num2str(toc(localTimer)) ' s']);
            localTimer = tic;
        end

        % ------------------------------------------------------------------------
        % Save the calculated data
        %-------------------------------------------------------------------------
        % Define reports path
        reportsPath = [obj.normalizedDataPath 'reports' filesep];
        if ~exist(reportsPath)
            mkdir(reportsPath);
        end
        
        % Define path to barcode_metadata
        barcodeMetadataPath = [reportsPath 'barcode_metadata.csv'];

        % Display progress
        if obj.verbose
            PageBreak();
            disp(['Converting to table']);
            conversionTimer = tic;
        end
        % Convert to table
        T = struct2table(combinedList);
        
        % Display progress
        if obj.verbose
            disp(['...completed in ' num2str(toc(conversionTimer)) ' s']);
            disp(['Writing barcode metadata: ' barcodeMetadataPath]);
            writeTimer = tic;
        end
        
        % Write table to file
        writetable(T, barcodeMetadataPath);
        
        % Display progress
        if obj.verbose
            disp(['..completed write in ' num2str(toc(writeTimer)) ' s']);
            PageBreak();
            disp(['...completed feature counts in ' num2str(toc(totalTimer)) ' s']);
        end
    end

    % -------------------------------------------------------------------------
    % Calculate barcode counts per feature
    % -------------------------------------------------------------------------
    function CalculateFeatureCounts(obj, varargin)
        % Calculate the number of barcodes per feature
        % obj.CalculateFeatureCounts();
        
        % Display progress
        if obj.verbose
            PageBreak();
            disp(['Compute barcode counts for all features']);
            localTimer = tic;
        end
        
        % -------------------------------------------------------------------------
        % Confirm all files are present
        % -------------------------------------------------------------------------
        parsedBarcodePath = [obj.normalizedDataPath obj.barcodePath filesep 'parsed_fov' filesep];
        if ~all(obj.CheckStatus([], 'p'))
            error('matlabFunctions:missingData', 'Some fov have not yet been parsed!');
        end
   
        % Confirm that the final found features are present
        if ~all(obj.CheckStatus([], 'c'))
            error('matlabFunctions:missingData', 'Combined found features appear to be missing.');
        end
        
        % -------------------------------------------------------------------------
        % Load feature boundaries
        % -------------------------------------------------------------------------
        if obj.verbose
            disp(['Loading found features']);
        end
        foundFeatures = obj.GetFoundFeatures();
        
        % Determine number of unique feature ids
        featureIDs = unique([foundFeatures.feature_id]);
        numFeatures = length(featureIDs);
        
        if obj.verbose
            disp(['...found ' num2str(length(featureIDs)) ' features']);
        end
        
        % Transfer some information to local variables
        numBarcodes = obj.numBarcodes;
        numFov = obj.numFov;
        fovIDs = obj.fovIDs;
        minDist = obj.parameters.quantification.minimumDistanceToFeature;
        
        % Setup parallel processes
        spmd (obj.numPar)
            % Define local variables for accumulation per worker
            countsPerCellExactIn = zeros(numBarcodes, numFeatures);
            countsPerCellCorrectedIn = zeros(numBarcodes, numFeatures);
            countsPerCellExactOut = zeros(numBarcodes, numFeatures);
            countsPerCellCorrectedOut = zeros(numBarcodes, numFeatures);
            
            % Loop over a set of blocks per worker
            for b=labindex:numlabs:numFov
                % Determine local fov
                localFovID = fovIDs(b);
                
                % Display progress
                if obj.verbose
                    displayStrings = {};
                    displayStrings{end+1} = PageBreak('nodisplay');
                    displayStrings{end+1} = ['Counting for fov ' obj.fov2str(localFovID)];
                    loadTimer = tic;
                end
                
                % Load barcodes
                aList = ReadBinaryFile([parsedBarcodePath 'fov_' obj.fov2str(localFovID) '_blist.bin']);
                
                % Check for empty barcodes
                if isempty(aList)
                    continue;
                end

                % Cut List
                aList = aList([aList.area] >= obj.parameters.quantification.minimumBarcodeArea & ...
                    double([aList.total_magnitude])./double([aList.area]) >= obj.parameters.quantification.minimumBarcodeBrightness);
                
                % Cut list to a specific Z range if requested
                if ~isempty(obj.parameters.quantification.zSliceRange)
                    pos = cat(1,[aList.abs_position]);
                    zPos = pos(:,3);
                    aList = aList(zPos >= obj.parameters.quantification.zSliceRange(1) && ...
                        zPos <= obj.parameters.quantification.zSliceRange(2));
                end

                if obj.verbose
                    displayStrings{end+1} = (['Loaded and cut ' num2str(length(aList)) ' barcodes in ' num2str(toc(loadTimer)) ' s']); 
                    countTimer = tic;
                end
                
                % Check for empty barcodes
                if isempty(aList)
                    continue;
                end

                % Extract needed information for exact information
                isExact = [aList.is_exact]==1;
                barcodeID = [aList.barcode_id];
                inFeature = [aList.in_feature];
                featureDist = [aList.feature_dist];
                
                % Convert featureID to the corresponding index in the count
                % matrix
                featureInd = discretize([aList.feature_id], [featureIDs inf]); % The lower edge is <= the upper is <
                
%--------------------------------------------------------------------------
                featureInd = uint32(featureInd);
                barcodeID = uint32(barcodeID);
%--------------------------------------------------------------------------
                
                % Compute the various quantities
                % Exact/in feature
                data = cat(2, barcodeID(isExact & inFeature)', featureInd(isExact & inFeature)'); % Exact/in features
                if ~isempty(data)
                    countsPerCellExactIn = countsPerCellExactIn + hist3(data, ...
                        'Edges', {1:numBarcodes, 1:numFeatures});
                end
                
                % Corrected/in feature
                data = cat(2, barcodeID(~isExact & inFeature)', featureInd(~isExact & inFeature)'); % Corrected/in features
                if ~isempty(data)
                    countsPerCellCorrectedIn = countsPerCellCorrectedIn + hist3(data, ...
                        'Edges', {1:numBarcodes, 1:numFeatures});
                end

                % Exact/out of feature
                data = cat(2, barcodeID(isExact & ~inFeature & featureDist <= minDist)', featureInd(isExact & ~inFeature & featureDist <= minDist)'); % Corrected/in features
                if ~isempty(data)
                    countsPerCellExactOut = countsPerCellExactOut + hist3(data, ...
                        'Edges', {1:numBarcodes, 1:numFeatures});
                end

                % Exact/out of feature
                data = cat(2, barcodeID(~isExact & ~inFeature & featureDist <= minDist)', featureInd(~isExact & ~inFeature & featureDist <= minDist)'); % Corrected/in features
                if ~isempty(data)
                    countsPerCellCorrectedOut = countsPerCellCorrectedOut + hist3(data, ...
                        'Edges', {1:numBarcodes, 1:numFeatures});
                end

                if obj.verbose
                    displayStrings{end+1} = (['..completed in ' num2str(toc(countTimer)) ' s']);
                    display(char(displayStrings)); % Flush buffer
                end
            end
        end
        
        % -------------------------------------------------------------------------
        % Compile composite objects
        % -------------------------------------------------------------------------
        % Display progress
        if obj.verbose
            PageBreak();
            disp(['Combining counts from all parallel objects']);
            localTimer = tic;
        end

        temp = countsPerCellExactIn{1};
        for i=2:length(countsPerCellExactIn)
            temp = temp + countsPerCellExactIn{i};
        end
        countsPerCellExactIn = temp;

        temp = countsPerCellCorrectedIn{1};
        for i=2:length(countsPerCellCorrectedIn)
            temp = temp + countsPerCellCorrectedIn{i};
        end
        countsPerCellCorrectedIn = temp;

        temp = countsPerCellExactOut{1};
        for i=2:length(countsPerCellExactOut)
            temp = temp + countsPerCellExactOut{i};
        end
        countsPerCellExactOut = temp;

        temp = countsPerCellCorrectedOut{1};
        for i=2:length(countsPerCellCorrectedOut)
            temp = temp + countsPerCellCorrectedOut{i};
        end
        countsPerCellCorrectedOut = temp;

        if obj.verbose
            disp(['...completed in ' num2str(toc(localTimer)) ' s']);
            localTimer = tic;
        end

        % ------------------------------------------------------------------------
        % Save the calculated data
        %-------------------------------------------------------------------------
        reportsPath = [obj.normalizedDataPath 'reports' filesep];
        if ~exist(reportsPath)
            mkdir(reportsPath);
        end

        if obj.verbose
            PageBreak();
            disp(['Writing data']);
        end
        
        % Write
        csvwrite([reportsPath 'countsPerCellExactIn.csv'], countsPerCellExactIn);
        csvwrite([reportsPath 'countsPerCellCorrectedIn.csv'], countsPerCellCorrectedIn);
        csvwrite([reportsPath 'countsPerCellExactOut.csv'], countsPerCellExactOut);
        csvwrite([reportsPath 'countsPerCellCorrectedOut.csv'], countsPerCellCorrectedOut);

        if obj.verbose
            disp(['...wrote ' reportsPath 'countsPerCellExactIn.csv']);
            disp(['...wrote ' reportsPath 'countsPerCellCorrectedIn.csv']);
            disp(['...wrote ' reportsPath 'countsPerCellExactOut.csv']);
            disp(['...wrote ' reportsPath 'countsPerCellCorrectedOut.csv']);
            disp(['Completed compiling performance statistics at ' datestr(now)]);
        end 
        
        % ------------------------------------------------------------------------
        % Write metadata
        %-------------------------------------------------------------------------
        % Define the feature name file
        featuresNameFilePath = [reportsPath 'featureNames.csv'];
        
        % Generate the contents of the feature names file
        featureUIDs = {foundFeatures.uID}; % Feature ids
        
        % Open and write file
        fid = fopen(featuresNameFilePath, 'W');
        fprintf(fid, '%s\n', featureUIDs{:});
        fclose(fid);
        
        % Display progress
        if obj.verbose
            disp(['...wrote ' featuresNameFilePath]);
        end
        
        % Define the gene name file
        geneNamesFilePath = [reportsPath 'geneNames.csv'];
        geneNames = {obj.codebook.name};

        % Open and write file
        fid = fopen(geneNamesFilePath, 'W');
        fprintf(fid, '%s\n', geneNames{:});
        fclose(fid);

        % Display progress
        if obj.verbose
            disp(['...wrote ' featuresNameFilePath]);
        end

        % Display progress
        if obj.verbose
            PageBreak();
            disp(['...completed feature counts in ' num2str(toc(localTimer)) ' s']);
        end
    end
    
    % -------------------------------------------------------------------------
    % Create found features report
    % -------------------------------------------------------------------------
    function GenerateFoundFeaturesReport(obj)
        % Generate a report on the segmented and joined features
        % 
        % obj.GenerateFoundFeaturesReport();
        
        % Display progress
        if obj.verbose
            PageBreak();
            disp(['Generating the found features reports']);
        end
        
        % Check status
        if ~all(obj.CheckStatus([], 'combine'))
            error('matlabFunctions:missingData', 'Feature segmentation and combination has not yet been completed!');
        end
        
        % Define paths to found features
        allFoundFeaturesPath = [obj.normalizedDataPath obj.segmentationPath 'all_found_features.matb'];
        finalFoundFeaturesPath = [obj.normalizedDataPath obj.segmentationPath 'final_found_features.matb'];

        % Load the original found features to determine some basic stats        
        if obj.verbose
            disp(['...loading original found features: ' allFoundFeaturesPath]);
            loadTimer = tic;
        end

        allFoundFeatures = LoadSplitByteStream(allFoundFeaturesPath, 'verbose', obj.verbose);
        
        if obj.verbose
            disp(['...loaded ' num2str(length(allFoundFeatures)) ' features in ' num2str(toc(loadTimer)) ' s']);
            loadTimer = tic;
        end

        % Save some basic statistics
        numOriginalFeatures = length(allFoundFeatures);
        numBroken = sum([allFoundFeatures.is_broken] == 1);
        numUnbroken = sum([allFoundFeatures.is_broken] == 0);
        zPos = allFoundFeatures(1).abs_zPos;
        
        % Extract area
        area = cat(1, allFoundFeatures.boundary_area);
        abs_area = cat(1, allFoundFeatures.abs_boundary_area);
        
        % Extract volume
        volume = [allFoundFeatures.volume];
        
        % Create figure handle
        figHandle = figure('Name', 'Found feature statistics', ...
            'Color', 'w', 'Position', [1 1 1000 500]);
        
        % Generate subplot on feature numbers
        subplot(2,3,1);
        bar([numOriginalFeatures numBroken numUnbroken]);
        set(gca, 'XTick', 1:3, 'XTickLabels', {'All', 'Broken', 'Unbroken'}, ...
            'XTickLabelRotation', 90);
        ylabel('Number of features');
        
        % Generate subplot on feature area by z position
        subplot(2,3,2);
        avArea = mean(area,1);
        errArea = std(area,0,1)./sqrt(length(avArea));
        bar(zPos, avArea, 'b'); hold on;
        errorbar(zPos, avArea, errArea, 'k.');
        xlabel('Z Position (micron)');
        ylabel('Average number of voxels');
        title([num2str(mean(area(:)),2) '+/-' num2str(std(area(:))/sqrt(length(area(:))),2) ' (SEM)']);
        if length(zPos) > 1
            zDiff = mean(diff(zPos));
            xlim([min(zPos) max(zPos)] + zDiff*[-1 1]);
        end
        
        % Plot distribution of volume
        subplot(2,3,3);
        [n,x] = hist(volume, 1000);
        bar(x, n);
        xlabel('Volume (voxels)');
        ylabel('Count');
        title([num2str(mean(volume),2) '+/-' num2str(std(volume)/sqrt(length(volume)),2) ' (SEM)']);

        % Clear memory for the original found features
        allFoundFeatures = [];
        
        % Load the final found features        
        if obj.verbose
            disp(['...loading final found features: ' finalFoundFeaturesPath]);
            loadTimer = tic;
        end

        finalFoundFeatures = LoadSplitByteStream(finalFoundFeaturesPath, 'verbose', obj.verbose);
        
        if obj.verbose
            disp(['...loaded ' num2str(length(finalFoundFeatures)) ' features in ' num2str(toc(loadTimer)) ' s']);
            loadTimer = tic;
        end

        % Save some basic statistics
        numOriginalFeatures = length(finalFoundFeatures);
        numSelfJoin = sum([finalFoundFeatures.num_joined_features] == 1);
        numDoubleJoin = sum([finalFoundFeatures.num_joined_features] == 2);
        
        % Extract area
        area = cat(1, finalFoundFeatures.abs_boundary_area);
        
        % Extract volume
        volume = [finalFoundFeatures.abs_volume];
                
        % Generate subplot on feature numbers
        subplot(2,3,4);
        bar([numOriginalFeatures numSelfJoin numDoubleJoin]);
        set(gca, 'XTick', 1:3, 'XTickLabels', {'All', 'Self-Join', 'Double-Join'}, ...
            'XTickLabelRotation', 90);
        ylabel('Number of features');
        
        % Generate subplot on feature area by z position
        subplot(2,3,5);
        avArea = mean(area,1);
        errArea = std(area,0,1)./sqrt(length(avArea));
        bar(zPos, avArea, 'b'); hold on;
        errorbar(zPos, avArea, errArea, 'k.');
        xlabel('Z Position (micron)');
        ylabel('Average area (microns^2)');
        title([num2str(mean(area(:)),2) '+/-' num2str(std(area(:))/sqrt(length(area(:))),2) ' (SEM)']);
        if length(zPos) > 1
            zDiff = mean(diff(zPos));
            xlim([min(zPos) max(zPos)] + zDiff*[-1 1]);
        end
        
        % Plot distribution of volume
        subplot(2,3,6);
        [n,x] = hist(volume, 1000);
        bar(x, n);
        xlabel('Volume (microns^3)');
        ylabel('Count');
        title([num2str(mean(volume),2) '+/-' num2str(std(volume)/sqrt(length(volume)),2) ' (SEM)']);
        
        % Save report figure
        SaveFigure(figHandle, 'parameters', obj.parameters.display, ...
            'savePath', [obj.normalizedDataPath obj.reportPath]);
        
    end
    
    % -------------------------------------------------------------------------
    % Generate Summation Report
    % -------------------------------------------------------------------------
    function GenerateSummationReport(obj)
        % Generate a report on the summation of features
        % 
        % obj.GenerateSummationReport();
        
        % Display progress
        if obj.verbose
            PageBreak();
            disp(['Generating a report on the summation of raw data channels']);
        end
        
        % Check status: combination of raw sum
        if ~obj.CheckStatus(obj.fovIDs(1), 'i')
            error('matlabFunctions:missingData', 'Raw signal summation has not yet been completed!');
        end
        
        % Check status for mosaics
        useMosaic = false;
        if obj.CheckStatus(obj.fovIDs(1), 'l')
           useMosaic = true;
           if obj.verbose
               disp(['Use mosaic: ' num2str(useMosaic)]);
           end
        end
        
        % Load found features
        foundFeatures = obj.GetFoundFeatures();
        foundFeatureFovIDs = foundFeatures.ReturnPrimaryFovID();
        
        % Load the summation data
        [normalizedSignal, ~,~,dataChannelNames] = obj.GetSummedSignal();
                
        % Create figure handle
        figHandle = figure('Name', 'Summation statistics', ...
            'Color', 'w');
        
        % Create a master approach
        for d=1:length(dataChannelNames)
            % Extract local data
            localData = normalizedSignal(d,:);
            
            % Create kernel estimate
            [p,x] = ksdensity(localData, linspace(0, 2*quantile(localData,0.95), 250), ...
                'support', [0 Inf]);
       
            % Create scatter plot
            plot(d + 0.75*(rand(1, length(localData))-0.5), ...
                localData, '.', 'Color', [.9 .9, .9]); hold on;
        
            % Create violin
            plot(d + 0.75*[p -fliplr(p)]/max(p)/2, [x fliplr(x)]); hold on;
        end
        
        % Set axis labels
        set(gca, 'XTick', 1:length(dataChannelNames), 'XTickLabel', dataChannelNames, ...
            'XTickLabelRotation', 90);
        ylabel('Normalized Signal');
        
        % Save report figure
        SaveFigure(figHandle, 'parameters', obj.parameters.display, ...
            'savePath', [obj.normalizedDataPath obj.reportPath]);
        % Close figure handle
        close(figHandle);
        
        % Generate distribution reports for all cell types
        numColors = 25;
        cMap = jet(numColors);

        % Determine properties of slices
        numSlices = length(obj.sliceIDs);
        
        % Loop over slices
        for s=1:numSlices
            % Check to confirm that the dataset contains the fovIDs for
            % this slice
            if ~isempty(setdiff(obj.sliceIDs{s}, obj.fovIDs))
                continue; % Not all fovIDs are present, skip analysis
            end
            
            % Extract the proper features
            goodFeatureInds = ismember(foundFeatureFovIDs, obj.sliceIDs{s});
            goodFeatures = foundFeatures(goodFeatureInds);
            
            % Load the slice mosaic (if requested)
            if useMosaic
                [mosaicImageStack, coordinates] = obj.GetMosaic(s);
            end
            
            % Loop over data channels
            for d=1:length(dataChannelNames)
                % Extract local data
                localData = normalizedSignal(d,goodFeatureInds);

                % Discretize the data into specific bins
                edges = [-Inf linspace(quantile(localData, 0.05), quantile(localData, 0.95), numColors-1) Inf];

                colorID = discretize(localData, edges);

                % Create figure
                figHandle = figure('Name', ['Cell distributions for ' dataChannelNames{d}], ...
                    'Color', 'w');

                % Add the mosaic
                if useMosaic                    
                    imshow(mosaicImageStack(:,:,obj.parameters.summation.dcIndsForSummation(d)), ...
                        [], ...
                        'XData', coordinates.xLimits, ...
                        'YData', coordinates.yLimits);
                    hold on;
                end
                
                % Loop over the colorIDs
                for c=1:numColors
                    localFeatures = goodFeatures(colorID == c);

                    pos = zeros(0,2);

                    for g=1:length(localFeatures)
                        localPos = localFeatures(g).abs_boundaries{ceil(localFeatures(g).num_zPos/2)};
                        pos = cat(1, pos, localPos(1:obj.parameters.display.downSample:end,:), nan(1,2));
                    end

                    plot(pos(:,1), pos(:,2), 'Color', cMap(c,:)); hold on;

                    drawnow;
                end

                % Format figure
                xlabel('X Position (microns)');
                ylabel('Y Position (microns)');

                % Save report figure
                SaveFigure(figHandle, 'parameters', obj.parameters.display, ...
                    'savePath', [obj.normalizedDataPath obj.reportPath filesep 'summation_distributions' filesep 'slice_' num2str(s) filesep]);
                close(figHandle);
            end
        end
        
    end
    
    % -------------------------------------------------------------------------
    % Create feature counts report
    % -------------------------------------------------------------------------
    function GenerateFeatureCountsReport(obj)
        % Generate a report on the counts per feature
        % 
        % obj.GenerateFeatureCountsReport();
        
        % Display progress
        if obj.verbose
            PageBreak();
            disp(['Generating the feature count reports']);
        end
        
        % Check for the existance of the feature reports
        reportsPath = [obj.normalizedDataPath 'reports' filesep];
        if (~exist([reportsPath 'countsPerCellExactIn.csv'], 'file') || ...
            ~exist([reportsPath 'countsPerCellExactIn.csv'], 'file') || ...
            ~exist([reportsPath 'countsPerCellExactIn.csv'], 'file') || ...
            ~exist([reportsPath 'countsPerCellExactIn.csv'], 'file'))
                error('matlabFunctions:missingData', 'The feature counts must be calculated first.');
        end
            
        % Load feature counts
        cInE = csvread([reportsPath 'countsPerCellExactIn.csv']);
        cInC = csvread([reportsPath 'countsPerCellCorrectedIn.csv']);
        cOutE = csvread([reportsPath 'countsPerCellExactOut.csv']);
        cOutC = csvread([reportsPath 'countsPerCellCorrectedOut.csv']);

        % Load features
        finalFoundFeaturesPath = [obj.normalizedDataPath obj.segmentationPath 'final_found_features.matb'];

        if obj.verbose
            disp(['...loading final found features: ' finalFoundFeaturesPath]);
            loadTimer = tic;
        end

        finalFoundFeatures = LoadSplitByteStream(finalFoundFeaturesPath, 'verbose', obj.verbose);
        
        if obj.verbose
            disp(['...loaded ' num2str(length(finalFoundFeatures)) ' features in ' num2str(toc(loadTimer)) ' s']);
            loadTimer = tic;
        end
        
        % Save basic counts/feature statistics
        figHandle = figure('Name', 'Counts per feature', 'Color', 'w');
    
        subplot(2,3,1);
        localData = sum(cInC,1);
        hist(localData, 0:1:500);
        xlabel('Counts per cell');
        ylabel('Number of cells');
        title(['CI: ' num2str(mean(localData),2) '\pm' num2str(std(localData),2)]);
        xlim([0 500]);
        set(gca, 'YScale', 'log');

        subplot(2,3,2);
        localData = sum(cInE,1);
        hist(localData,  0:1:1000);
        xlabel('Counts per cell');
        ylabel('Number of cells');
        title(['EI: ' num2str(mean(localData),2) '\pm' num2str(std(localData),2)]);
        xlim([0 500]);
        set(gca, 'YScale', 'log');

        subplot(2,3,3);
        localData = sum(cInE+cInC,1);
        hist(localData,  0:1:1000);
        xlabel('Counts per cell');
        ylabel('Number of cells');
        title(['All in: ' num2str(mean(localData),2) '\pm' num2str(std(localData),2)]);
        xlim([0 500]);
        set(gca, 'YScale', 'log');

        subplot(2,3,4);
        localData = sum(cOutC,1);
        hist(localData,  0:1:1000);
        xlabel('Counts per cell');
        ylabel('Number of cells');
        title(['CO: ' num2str(mean(localData),2) '\pm' num2str(std(localData),2)]);
        xlim([0 500]);
        set(gca, 'YScale', 'log');

        subplot(2,3,5);
        localData = sum(cOutE,1);
        hist(localData,  0:1:1000);
        xlabel('Counts per cell');
        ylabel('Number of cells');
        title(['EO: ' num2str(mean(localData),2) '\pm' num2str(std(localData),2)]);
        xlim([0 1000]);
        set(gca, 'YScale', 'log');

        subplot(2,3,6);
        localData = sum(cOutE,1);
        hist(localData,  0:1:2000);
        xlabel('Counts per cell');
        ylabel('Number of cells');
        title(['All out: ' num2str(mean(localData),2) '\pm' num2str(std(localData),2)]);
        xlim([0 1000]);
        set(gca, 'YScale', 'log');

        % Save report figure
        SaveFigure(figHandle, 'parameters', obj.parameters.display, ...
            'savePath', [obj.normalizedDataPath obj.reportPath]);
        
        % Save basic counts/volume stats
        figHandle = figure('Name', 'Counts per volume', 'Color', 'w');

        % Extract the volume values for each feature
        volume = [finalFoundFeatures.volume];
        
        % Volume versus all counts
        subplot(1,2,1);
        localData = sum(cInC+cInE,1);
        plot(volume, localData, '.');
        xlabel('Volume (microns^3)');
        ylabel('Total counts within feature');
        
        subplot(1,2,2);
        localData = sum(cInC+cInE,1)./volume;
        localData = localData(~isnan(localData));
        [n, x] = hist(localData, linspace(0, 0.025, 100));
        bar(x, n);
        xlabel('Density (counts/microns^3)');
        ylabel('Counts');
        xlim([-0.0001 0.03]);
        title(['Median: ' num2str(median(localData),2) '+/-' num2str(iqr(localData),2) ' (iqr)']);
        
        % Save report figure
        SaveFigure(figHandle, 'parameters', obj.parameters.display, ...
            'savePath', [obj.normalizedDataPath obj.reportPath]);

    end

    
    % -------------------------------------------------------------------------
    % Load image metadata
    % -------------------------------------------------------------------------
    function LoadImageMetaData(obj)
        % fovPos = MapFOVPositions()
        % Extract fov positions from raw data meta data and set image size
        
        % -------------------------------------------------------------------------
        % Check that raw file metadata has been loaded
        % -------------------------------------------------------------------------
        if isempty(obj.rawDataFiles)
            error('matlabFunctions:nonexistantData', 'Raw data files have not yet been loaded');
        end
        
        % -------------------------------------------------------------------------
        % Loop over fov in order of fovIDs and load stage pos
        % -------------------------------------------------------------------------
        for f=1:obj.numFov
            % Select the first data channel image file to determine
            % location of meta data
            
            % Handle camera
            if ~isfield(obj.dataOrganization(1), 'imagingCameraID')
                fileInd = find(...
                    strcmp({obj.rawDataFiles.imageType}, obj.dataOrganization(1).imageType) & ...
                    [obj.rawDataFiles.imagingRound] == obj.dataOrganization(1).imagingRound & ...
                    [obj.rawDataFiles.fov] == obj.fovIDs(f));
            else
                fileInd = find(...
                    strcmp({obj.rawDataFiles.imageType}, obj.dataOrganization(1).imageType) & ...
                    [obj.rawDataFiles.imagingRound] == obj.dataOrganization(1).imagingRound & ...
                    strcmp({obj.rawDataFiles.cameraID}, obj.dataOrganization(1).imagingCameraID) & ...
                    [obj.rawDataFiles.fov] == obj.fovIDs(f));
            end
            
            % Check for existance
            if isempty(fileInd)
                display(['Error: missing fov id: ' obj.fov2str(obj.fovIDs(f))]);
                error('matlabFunctions:missingFile', 'Could not find a raw data file');
            end
           
            % Switch on raw data type to determine how load metadata
            switch obj.imageExt
                case {'dax', 'tiff', 'tif'} 
                    
                    % Switch on the version of hal
                    switch obj.hal_version
                        case 'hal1'
                            if ismember(obj.imageExt, {'dax'})
                                % Load image meta data
                                infoFile = ReadInfoFile([obj.rawDataPath obj.rawDataFiles(fileInd).name], 'verbose', false);
                            else
                                % THIS CASE IS RESERVED FOR FUTURE USE
                                error('matlabFunctions:unsupportedExt', 'Tiff data are not yet for hal1 info file structure');
                            end

                        case 'hal2'
                            % UNDER CONSTRUCTION: THIS NEEDS TO BE UDPDATED
                            % TO HAVE A ROBUST FUNCTION FOR LOADING THE XML
                            % FILE ASSOCIATED WITH THE IMAGE
                            
                            % Determine the filename and strip off the
                            % camera tag if necessary
                            baseFile = [obj.rawDataPath obj.rawDataFiles(fileInd).name];
                            baseFile = baseFile(1:(end-4)); % Strip off the dax extension
                            if ismember(baseFile((end-2):end), {'_c1', '_c2'})
                                baseFile = baseFile(1:(end-3));
                            end
                                                        
                            % Load the xml file
                            xDoc = xmlread([baseFile '.xml']);
                            
                            % Find the stage position node and get its data
                            allPosItems = xDoc.getElementsByTagName('stage_position');
                            thisListItem = allPosItems.item(0);
                            posString = thisListItem.getFirstChild.getData;
                            
                            pos = strsplit(char(posString), ',');
                            
                            infoFile.Stage_X = str2num(pos{1});
                            infoFile.Stage_Y = str2num(pos{2});
                            
                            % Find the height and width of the image to
                            % compute image size
                            xSizeList = xDoc.getElementsByTagName('x_pixels');
                            xSizeItem = xSizeList.item(0);
                            xSizeString = xSizeItem.getFirstChild.getData;
                            
                            infoFile.frame_dimensions(1) = str2num(char(xSizeString));
                            
                            ySizeList = xDoc.getElementsByTagName('y_pixels');
                            ySizeItem = ySizeList.item(0);
                            ySizeString = ySizeItem.getFirstChild.getData;
                            
                            infoFile.frame_dimensions(2) = str2num(char(ySizeString));

                        otherwise
                            error('matlabFunctions:unsupportedExt', 'An unrecognized hal version was provided');
                    end

                    % Archive stage position
                    obj.fovPos(f,:) = [infoFile.Stage_X infoFile.Stage_Y];  % Record stage x, y
                    
                    % Archive image size (only the first fov)
                    if f==1
                        obj.imageSize = infoFile.frame_dimensions;
                    end
                otherwise
                    error('matlabFunctions:unsupportedExt', 'This file extension is not supported');
            end
            
        end
    end
    
    % -------------------------------------------------------------------------
    % Map raw data file
    % -------------------------------------------------------------------------
    function foundFiles = MapRawData(obj)
        % foundFiles = mDecoder.MapRawData()
        % Map the organization of the raw data files in the raw data path
        % according to the data organization file
    
        % Display progress
        if obj.verbose
            PageBreak();
            display(['Finding all ' obj.imageExt ' files']);
        end

        % Extract unique image types
        [uniqueImageType, ia] = unique({obj.dataOrganization.imageType});
        display(['Found ' num2str(length(uniqueImageType)) ' image types']);
        for i=1:length(uniqueImageType)
            display(['...Parsing ' uniqueImageType{i} ' with: ' obj.dataOrganization(ia(i)).imageRegExp]);
        end

        % Loop over all file name patterns
        foundFiles = [];
        parseTimer = tic;
        for i=1:length(uniqueImageType)

            % Build file data structure
            newFiles = BuildFileStructure(obj.rawDataPath, ...
                'fileExt', obj.imageExt, ...
                'fieldNames', {'imageType', 'fov', 'imagingRound', 'cameraID'}, ... % Required fields
                'fieldConv', {@char, @str2num, @str2num, @char}, ...
                'regExp', obj.dataOrganization(ia(i)).imageRegExp, ...
                'requireFlag', uniqueImageType{i});

            % Coerce empty imageRound fields
            if any(arrayfun(@(x)isempty(x.imagingRound), newFiles))
                [newFiles(:).imagingRound] = deal(-1); % Flag that no image round was specified
            end
            
            % Coerce empty cameraID fields
            if any(arrayfun(@(x)isempty(x.cameraID), newFiles))
                [newFiles(:).cameraID] = deal(''); % Flag that no image round was specified
            end

            % Combine file structures
            foundFiles = [foundFiles newFiles];
        end
        disp(['...completed parse in ' num2str(toc(parseTimer)) ' s']);
        
        % Check to see if any files were found
        if isempty(foundFiles)
            error('matlabFunctions:noFoundFiles', 'No files matching the patterns in the data organization file were found!');
        end

        % Handle the occasional parsing error
        indsToKeep = ismember({foundFiles.imageType}, uniqueImageType);
        if any(~indsToKeep) && obj.verbose
            display(['...found parsing errors']);
            display(['...removing ' num2str(sum(~indsToKeep)) ' files']);
        end
        foundFiles = foundFiles(indsToKeep);

		% Remove replicates
		[~, ia] = unique({foundFiles.name}); % Find all unique file names
		if obj.verbose
			display(['...removing ' num2str(length(foundFiles) - length(ia)) ' replicated files']);
		end
		foundFiles = foundFiles(ia);
		
        % Compile properties
        obj.fovIDs = unique([foundFiles.fov]);
        obj.imageRoundIDs = unique([foundFiles.imagingRound]);
        obj.numFov = length(obj.fovIDs);
        obj.numImagingRounds = length(obj.imageRoundIDs);
        obj.cameraIDs = unique({foundFiles.cameraID});
        obj.numCameraIDs = length(obj.cameraIDs);
        
        % Allocate memory for stage positions
        obj.fovPos = nan(obj.numFov, 2); % These entries are in the same order as the fovIDs.
        
        % Display progress
        if obj.verbose
            display(['Found ' num2str(length(foundFiles)) ' files']);
            display(['...' num2str(obj.numFov) ' fov']);
            display(['...' num2str(obj.numImagingRounds) ' imaging rounds']);
            PageBreak();
            display(['Checking consistency of raw data files and expected data organization']);
            localTimer = tic;
        end        
        
        % Define useful fov to string conversion function (uniform padding)
        padNum2str = @(x,y)num2str(x,['%0',num2str(ceil(log10(y+1))),'d']);
        obj.fov2str = @(x)padNum2str(x, max(obj.fovIDs));
        
        %Run a cross check
        % Loop over each data channcel for each fov
        for c=1:obj.numDataChannels
            % Extract local properties of the channel
            localImageType = obj.dataOrganization(c).imageType;
            localImagingRound = obj.dataOrganization(c).imagingRound;
            if isfield(obj.dataOrganization(c), 'imagingCameraID') % Handle backwards compatibility or single camera systems
                localCameraID = obj.dataOrganization(c).imagingCameraID;
            else
                localCameraID = '';
            end

            % Extract the fovIDs for files with these parameters
            foundFovIDs = [foundFiles(strcmp({foundFiles.imageType}, localImageType) & ...
                [foundFiles.imagingRound] == localImagingRound & ...
                strcmp({foundFiles.cameraID}, localCameraID)).fov];

            % Find the fovIDs that are missing
            missingFovIDs = setdiff(foundFovIDs, obj.fovIDs);

            for f=1:length(missingFovIDs)
                disp(['An error has been found with the following file!']);
                disp(['   FOV: ' num2str(missingFovIDs)]);
                disp(['   imageType: ' obj.dataOrganization(c).imageType]);
                disp(['   imagingRound: ' obj.dataOrganization(c).imagingRound]);
                disp(['   imagingCameraID: ' localCameraID]);
                error('matlabFunctions:invalidFileInformation', ...
                    'Either a file is missing or there are multiple files that match an expected pattern.');
            end
            
            %Display Progress
            disp(['...completed ' num2str(c) ' channel of ' num2str(obj.numDataChannels)]);
            
        end
        
        % Display progress
        if obj.verbose
            disp(['...completed in ' num2str(toc(localTimer)) ' s']);
            disp(['...no problems were found.']);
        end        

    end
    
    % -------------------------------------------------------------------------
    % Segment features (cells) within individual fov
    % -------------------------------------------------------------------------
    function SegmentFOV(obj, fovIDs)
        % Segment features, e.g. cells, in the specified fov 
        % SegmentFOV([]);       % Segment all fov
        % SegmentFOV(fovIDs);   % Segment the fov that match the specified fovids
        
        % -------------------------------------------------------------------------
        % Determine properties of the requested fov ids
        % -------------------------------------------------------------------------
        if isempty(fovIDs)
            fovIDs = obj.fovIDs;
        elseif ~all(ismember(fovIDs, obj.fovIDs))
            error('matlabFunctions:invalidArguments', 'An invalid fov id has been requested');
        end
        
        % -------------------------------------------------------------------------
        % Create segmentation directory if needed
        % -------------------------------------------------------------------------
        if ~exist([obj.normalizedDataPath obj.segmentationPath], 'dir')
            mkdir([obj.normalizedDataPath obj.segmentationPath]);
        end
        
        % -------------------------------------------------------------------------
        % Extract local copy of parameters
        % -------------------------------------------------------------------------
        parameters = obj.parameters.segmentation;
        
        % -------------------------------------------------------------------------
        % Extract necessary parameters for segmentation
        % -------------------------------------------------------------------------
        switch parameters.segmentationMethod
            case 'seededWatershed'
                % Identify seed and watershed frames based on data
                % organization file
                dataChannelNames = {obj.dataOrganization.bitName};
                
                [~, ID] = ismember(parameters.watershedSeedChannel, dataChannelNames);
                if length(ID) > 1 || any(ID == 0)
                    error('matlabFunctions:invalidDataType', 'The requested seed channel has problems');
                end
                seedFrames = ((ID-1)*obj.numZPos +1):(ID*obj.numZPos);
                
                [~, ID] = ismember(parameters.watershedChannel, dataChannelNames);
                if length(ID) > 1 || any(ID == 0)
                    error('matlabFunctions:invalidDataType', 'The requested watershed channel has problems');
                end
                watershedFrames = ((ID-1)*obj.numZPos +1):(ID*obj.numZPos);
                                
            otherwise
                error('matlabFunctions:invalidArguments', 'The provided segmentation method is not currently supported.');
        end
        
        % Handle the request to ignore z
        if parameters.ignoreZ
            seedFrames = seedFrames(1);
            watershedFrames = watershedFrames(1);
        end
        
        % Extract parameters
        parameters = obj.parameters.segmentation;
        
        % -------------------------------------------------------------------------
        % Run processing on individual fov in parallel (if requested)
        % -------------------------------------------------------------------------
        spmd (obj.numPar) % Run in parallel
            % Loop over requested fov
            for f=labindex:numlabs:length(fovIDs)
				% Determine local fov id
                localFovID = fovIDs(f);
								
                % Create display strings
                if obj.verbose
                    displayStrings = {};
                    displayStrings{end+1} = PageBreak('nodisplay'); 
                    displayStrings{end+1} = ['Started segmentation for fov ' obj.fov2str(localFovID)];
                    fovTimer = tic;
                    localTimer = tic;
                end
                
                % Create file path and check for existance
                foundFeaturesPath = [obj.normalizedDataPath obj.segmentationPath ...
                    'found_features_fov_'  obj.fov2str(localFovID) '.matb'];

                if exist(foundFeaturesPath, 'file')
                    if obj.overwrite % If overwrite, then delete and overwrite
                        delete(foundFeaturesPath);
                        if obj.verbose
                            displayStrings{end+1} = ['...overwriting existing analysis'];
                        end
                    else
                        if obj.verbose
                            displayStrings{end+1} = ['...found existing analysis. Skipping.'];
                        end
                        continue;
                    end
                end                
                                    
                % Display progress
                if obj.verbose
                    displayStrings{end+1} = ['...using the ' parameters.segmentationMethod ' method'];
                end
                
                
                % Define tiff file name
                tiffFileName = [obj.normalizedDataPath obj.warpedDataPath 'fov_' obj.fov2str(localFovID) '.tif'];
                         
                % Clear any previous figure handles
                figHandles = [];
                
                % Allocate memory for the necessary frames
                localWatershedFrames = zeros([obj.imageSize length(watershedFrames)], 'uint16');
                localSeedFrames = zeros([obj.imageSize length(watershedFrames)], 'uint16');
                
                % Load and preprocess the necessary frames
                for z=1:length(watershedFrames)
                    % Load frames for watershed and seed generation
                    localWatershedFrames(:,:,z) = imread(tiffFileName, watershedFrames(z));
                    localSeedFrames(:,:,z) = imread(tiffFileName, seedFrames(z));
                    
                    
                    % Filter seed frame (if requested)
                    if ~isempty(parameters.seedFrameFilterSize)
                        localSeedFrames(:,:,z) = imgaussfilt(localSeedFrames(:,:,z), parameters.seedFrameFilterSize);
                    end

                    % Filter watershed image
                    if ~isempty(parameters.watershedFrameFilterSize)
                        localWatershedFrames(:,:,z) = imgaussfilt(localWatershedFrames(:,:,z), parameters.watershedFrameFilterSize);
                    end
                    
                    % Generate report
                    if parameters.saveSegmentationReports
                        % Generate a figure handle for this z slice
                        figHandles(z) = figure('Name', ['Segmentation report fov_' obj.fov2str(localFovID) ' z_' num2str(z)], ...
                            'Color', 'w');

                        % Update segmentation report
                        subplot(2,3,1);
                        imshow(localSeedFrames(:,:,z), []); hold on;
                        title('Filtered seed frame');
                    end
                end
                
                % Erode seed frame (if requested)
                if ~isempty(parameters.seedFrameErosionKernel)
                    localSeedFrames = imerode(localSeedFrames, parameters.seedFrameErosionKernel);
                end
                
                % Create a mask image
                if isempty(parameters.seedThreshold) % Use Otsu's method to determine threshold
                    mask = (localSeedFrames-min(localSeedFrames(:)))/(max(localSeedFrames(:)) - min(localSeedFrames(:))) >= ...
                        graythresh(localSeedFrames);
                elseif ischar(parameters.seedThreshold) && strcmp(parameters.seedThreshold, 'adaptive') % Use an adaptive method
                    for z=1:length(watershedFrames)
                        mask(:,:,z) = imbinarize(localSeedFrames(:,:,z), 'adaptive');
                    end
                else % Use a constant user specified value
                    mask = localSeedFrames >= parameters.seedThreshold;
                end
                                                
                % Create a 3D mask
                localSeedFrames(~mask) = 0;

                % Find regional max
                seeds = imregionalmax(localSeedFrames);
                    
                % Connect seeds that are close
                if ~isempty(parameters.seedConnectionKernel)
                    seeds = imdilate(seeds, parameters.seedConnectionKernel);
                end
                
                % Update the segmentation reports
                if parameters.saveSegmentationReports

                    for z=1:length(watershedFrames)
                        % Set the active figure
                        figure(figHandles(z));

                        % Update debug display
                        subplot(2,3,2);
                        %imshow(imoverlay(imadjust(localSeedFrames(:,:,z)), ~mask(:,:,z), 'red'));
                        imshow(localSeedFrames(:,:,z), []);
                        title('Eroded, masked seed frame');

                        ax(3) = subplot(2,3,3);
                        imshow(imoverlay(imadjust(localSeedFrames(:,:,z)), seeds(:,:,z), 'red'));
                        title('Seed frame with seed locations');

                        % Update debug display
                        ax(4) = subplot(2,3,4);
                        imshow(localWatershedFrames(:,:,z), []); hold on;
                        title('Watershed frame'); hold on;

                    end
                end
                
                % Parse the seeds to find seeds that are independent in one
                % frame and joined in another
                CC = bwconncomp(seeds);
                
                % Clear the previous seeds (they will be added back below)
                seeds = false(size(seeds));
                
                % Loop over all found 3D seeds and examine them in 2D
                for c=1:CC.NumObjects
                    % Prepare an image of just this seed
                    seedImage = false(size(seeds));
                    seedImage(CC.PixelIdxList{c}) = 1;
                    
                    % Determine the properties of this seed in each slice
                    localProps = {};
                    for z=1:length(watershedFrames)
                        localProps{z} = regionprops(seedImage(:,:,z),'Centroid');
                    end
                    
                    % Determine how many unique regions are in each slice
                    numSeeds = cellfun(@length, localProps);
                    
                    % And build a consensus centroid for each of them
                    if all(numSeeds < 2) % A single seed in all frames (or no seed)
                        % Extract seed centroids: the median position in
                        % the consensus position
                        goodFrames = numSeeds==1;
                        allGoodProps = cat(1,localProps{goodFrames});
                        seedPos = round(median(cat(1, allGoodProps.Centroid),1));
                    else % In at least one frame there is a double seed that has been combined in another frame
                        % Find the mediod using kmediods and the maximum
                        % number of observed seeds in any frame
                        goodFrames = numSeeds>1;
                        allGoodProps = cat(1,localProps{goodFrames});
                        [~, seedPos] = kmedoids(cat(1, allGoodProps.Centroid), max(numSeeds));
                        seedPos = round(seedPos);
                    end
                    
                    % Add back the seedPos in all z slices
                    for s=1:size(seedPos,1)
                        seeds(seedPos(s,2), seedPos(s,1),goodFrames) = 1;
                    end
                end
                
                % Dilate seed centroids
                if ~isempty(parameters.seedDilationKernel)
                    seeds = imdilate(seeds, parameters.seedDilationKernel); 
                end
                
                % Create in cell mask
                if isempty(parameters.watershedFrameThreshold)
                    % Calculate threshold using Otsu's method
                    level = graythresh(localWatershedFrames);
                    
                    % Calculate mask (recognizing that graythreshold
                    % preoduces value between 0 and 1)
                    inCellMask = (localWatershedFrames - min(localWatershedFrames(:)))/...
                        (max(localWatershedFrames(:)) - min(localWatershedFrames(:))) >= level;
                elseif ischar(parameters.watershedFrameThreshold) && strcmp(parameters.watershedFrameThreshold, 'adaptive')
                    for z=1:length(watershedFrames)
                        inCellMask(:,:,z) = imbinarize(localWatershedFrames(:,:,z), 'adaptive');
                    end
                else
                    inCellMask = localWatershedFrames >= parameters.watershedFrameThreshold;
                end
                
                % Fill in holes in the cell mask
                for z=1:size(inCellMask,3)
                    inCellMask(:,:,z) = imfill(inCellMask(:,:,z), 'holes');
                end
                    
                % Convert to double
                localWatershedFrames = double(localWatershedFrames);
                    
                % Set uniform range to watershed image
                localWatershedFrames = (localWatershedFrames - min(localWatershedFrames(:)))/...
                    (max(localWatershedFrames(:)) - min(localWatershedFrames(:)));
                    
                % Create image for watershed
                localWatershedFrames = imcomplement(localWatershedFrames); % Invert image
                localWatershedFrames(~inCellMask) = 0; % Mask with in cell mask
                localWatershedFrames(seeds) = 0; % Insert seeds for watershed

                % Remove local minima, forcing to only define watersheds via
                % imposed catch basins
                localWatershedFrames = imimposemin(localWatershedFrames, ~inCellMask | seeds);

                % Watershed
                L = watershed(localWatershedFrames);

                % Fill the empty, unassigned pixels
                for iter=1:100
                    % Check the stop early condition
                    if ~any(L(:)==0)
                        break;
                    end
                    
                    % First remove the regions not within cells
                    dilatedL = L;
                    dilatedL(~inCellMask) = 0;
                    
                    % Dilate the remaining regions by 1 pixel
                    dilatedL = imdilate(dilatedL, strel('disk',1));
                    
                    % Assign the values from the dilation process (uint8
                    % format forces selection one one label when two
                    % collide)
                    L(L==0) = dilatedL(L==0);
                end
                
                % Compute label matrix color map
                cMap = jet(double(max(L(:))));
                
                % Update the segmentation reports
                if parameters.saveSegmentationReports
                    for z=1:length(watershedFrames)
                        % Set the active figure
                        figure(figHandles(z));

                        subplot(2,3,5);
                        imshow(imoverlay(imadjust(localWatershedFrames(:,:,z)), ~inCellMask(:,:,z), 'red'));
                        title('Watershed frame');

                        subplot(2,3,6);
                        imshow(label2rgb(L(:,:,z), cMap, 'k'), []); hold on;
                        title('Label matrix');
                    end

                    % Link all children axes
                    allChildren = get(figHandles, 'Children');
                    if iscell(allChildren) % Handle the output for a array of figure handles
                        allChildren = cat(1,allChildren{:});
                    end
                    linkaxes(allChildren, 'xy');
                end

                % Display progress
                if obj.verbose
                    displayStrings{end+1} = ['...completed watershed in ' num2str(toc(localTimer)) ' s'];
                    localTimer = tic;
                end

                % Find the label(s) that corresponds to the regions not in
                % cells
                uniqueLabels = unique(L(:));
                isBadLabel = false(1, length(uniqueLabels));
                for J=1:length(isBadLabel)
                    isBadLabel(J) = any(L(:) == uniqueLabels(J) & ~inCellMask(:));
                end
                badLabels = uniqueLabels(isBadLabel);
                
                % Remove this label from feature construction
                L(ismember(L(:), badLabels)) = 0; % Label for no feature
                            
                % Determine the number of features
                uniqueLabels = unique(L(:));
                uniqueLabels = setdiff(uniqueLabels, 0); % Remove a marker of no features
                numFeatures = length(uniqueLabels);

                % Loop over all features creating a found feature
                foundFeatures = repmat(FoundFeature(), [0 1]);
                for n=1:numFeatures
                    % Convert label matrix into FoundFeatures
                    localFeature = FoundFeature(L == uniqueLabels(n), ...   % label matrix
                        localFovID, ...                                     % fovID
                        obj.fovPos(obj.fovIDs == localFovID,:), ...         % fov center position
                        obj.pixelSize, ...                                  % pixelSize
                        obj.parameters.decoding.stageOrientation, ...       % stage orientation 
                        obj.parameters.segmentation.boundingBox, ...        % bounding box
                        obj.zPos, ...                                       % z positions for each stack
                        uniqueLabels(n));                                   % the label associated with this feature                                          
                    
                    % Only keep features if they were not completely
                    % removed upon crop
                    if localFeature.abs_volume > 0
                        foundFeatures = cat(1, foundFeatures, localFeature);
                    end
                end
                
                % Display progress
                if obj.verbose 
                    displayStrings{end+1} = ['...found ' num2str(numFeatures) ' features in ' num2str(toc(localTimer)) ' s'];
                    displayStrings{end+1} = ['...keeping ' num2str(length(foundFeatures)) 'features after cuts by bounding box'];
                    localTimer = tic;
                end

                % Identify features that are contained within another feature
                % (These arise due to failure of a seed to drive the
                % creation of a watershed region larger than itself)
                isGoodFeature = true(1, length(foundFeatures));
                for i=1:length(foundFeatures)
                    doesOverlap = false(1, length(foundFeatures));
                    for j=1:length(foundFeatures)
                        % Skip the identity case
                        if i==j
                            continue;
                        end
                        % Check if feature j contains feature i
                        doesOverlap(j) = foundFeatures(j).DoesFeatureOverlap(foundFeatures(i));
                    end
                    isGoodFeature(i) = ~any(doesOverlap); % If no feature j contains feature i, then feature i is good
                end
                
                % Slice away features that are contained (at least
                % partially) within another feature
                foundFeatures = foundFeatures(isGoodFeature);

                % Display progress
                if obj.verbose 
                    displayStrings{end+1} = ['...completed cross checks for overlapping features in ' num2str(toc(localTimer)) ' s'];
                    displayStrings{end+1} = ['...keeping ' num2str(length(foundFeatures)) ' non-overlapping features'];
                    localTimer = tic;
                end
                
                % Create the path for the progress figures
                localSavePath = [obj.normalizedDataPath obj.segmentationPath filesep 'fov_images' filesep];
                if ~exist(localSavePath, 'dir')
                    mkdir(localSavePath);
                    displayStrings{end+1} = ['Created ' localSavePath];
                end
                
                % Update, save, and close the segmentation reports
                if parameters.saveSegmentationReports
                    for z=1:length(watershedFrames)
                        % Set figure
                        figure(figHandles(z));

                        % Plot all features
                        subplot(2,3,4);
                        hold on;
                        for F=1:length(foundFeatures)
                            plot(foundFeatures(F).boundaries{z}(:,1), foundFeatures(F).boundaries{z}(:,2), 'r');
                        end

                        % Force any residual draw
                        drawnow;

                        % Save figure
                        SaveFigure(figHandles(z), 'overwrite', true, 'savePath', localSavePath, ...
                            'formats', {'fig'}, 'verbose', false);

                        % Close figure
                        close(figHandles(z));
                    end
                end
                
                % Save the found features
                SaveAsByteStream(foundFeaturesPath, foundFeatures, 'verbose', false);
                
                % Flush display buffer
                if obj.verbose
                    displayStrings{end+1} = ['...saved ' foundFeaturesPath ' and figures in ' num2str(toc(localTimer)) ' s'];
                    displayStrings{end+1} = ['...completed fov ' obj.fov2str(localFovID) ' in ' num2str(toc(fovTimer)) ' s'];
                    display(char(displayStrings));
                end
             end
        end % SPMD loop
    end
    
    % -------------------------------------------------------------------------
    % Generate a low resolution mosaic stacks of the image data
    % -------------------------------------------------------------------------
    function GenerateLowResolutionMosaic(obj)
        % Generate a low resolution mosaic of the data using sliceIDs to
        % define different slices (if provided)
        %
        % obj.GenerateLowResolutionMosaic()
        
        % Handle the case that zInd was not provided
        zInd = obj.parameters.display.mosaicZInd;
        if zInd < 1 || zInd > obj.numZPos
            warning('matlabFunctions:invalidParameter', 'The provided mosaic z index is not within the z range. Assuming the first z position');
            zInd = 1;
        end
        
        % Determine the slice ids
        if isempty(obj.sliceIDs)
            sliceIDs = {obj.fovIDs};
        else
            sliceIDs = obj.sliceIDs;
        end
        numSlices = length(sliceIDs);
        
        % Determine the crop properties
        boundingBox = round(obj.parameters.segmentation.boundingBox/(obj.pixelSize/1000)); % Scale bounding box from microns to pixels
        startPixels = boundingBox(1:2) + obj.imageSize/2;
        endPixels = startPixels + boundingBox(3:4);

        % Display progress
        if obj.verbose
            PageBreak();
            disp(['Downsampling and saving ' num2str(length(sliceIDs)) ' mosaics']);
            totalTimer = tic;
        end
        
        % Check to see if warping has been completed
        doneWarping = obj.CheckStatus([], 'w');
        if ~all(doneWarping)
            disp('Some fov were not warped:');
            disp(find(~doneWarping));
            error('matlabFunctions:incompleteAnalysis', 'The data must be warped prior to construction of a low resolution mosaic stack');
        end
%         if ~all(obj.CheckStatus([], 'c')) 
%             if ~all(obj.CheckStatus([], 'w')) % Check for segmentation first since this check is faster than warping
%                 error('matlabFunctions:incompleteAnalysis', 'The data must be warped prior to construction of a low resolution mosaic stack');
%             end
%         end
                 
        % Make the downsample folder if necessary
        if ~exist([obj.normalizedDataPath obj.mosaicPath])
            mkdir([obj.normalizedDataPath obj.mosaicPath]);
        end
        
        % Loop over the individual slices
        for s=1:length(sliceIDs)
            % Display progress
            if obj.verbose
                disp(['...rendering channels for slice ' num2str(s)]);
                disp(['...preparing files']);
                localTimer = tic; 
            end
            
            % Determine the local fovs
            localFovIDs = obj.fovIDs(ismember(obj.fovIDs,sliceIDs{s}));
            
            % Check for existing slices and skip if not present
            if isempty(localFovIDs)
                if obj.verbose
                    disp(['Could not find the requested fov IDs']);
                end
                continue;
            end

            % Compute locations to determine ordering of files
            fovPositions = obj.fovPos(ismember(obj.fovIDs, localFovIDs),:); % Extract positions
            [~, ~, xInd] = unique(fovPositions(:,1)); % Map positions to row/column indices (1st, 2nd, etc...)
            [~, ~, yInd] = unique(fovPositions(:,2));
            
            numXFrames = max(xInd);
            numYFrames = max(yInd);
            
            % Determine the size of the final image by a test downsample
            testImage = zeros(obj.imageSize); % Generate empty image
            testImage = testImage(startPixels(1):endPixels(1), ...
                startPixels(2):endPixels(2)); % Crop
            testImage = imresize(testImage, 1/obj.parameters.display.downSample); % Resize
            [numXPixels, numYPixels] = size(testImage); % Determine final downsampled size
            
            numFinalPixelsX = numXPixels*numXFrames;
            numFinalPixelsY = numYPixels*numYFrames;
            
            % Determine the coordinate system of the final image
            coordinates.xLimits = [min(fovPositions(:,1)) max(fovPositions(:,1))] + ...
                obj.parameters.segmentation.boundingBox(1) + ...
                [0 obj.parameters.segmentation.boundingBox(3)] + ...
                obj.parameters.display.downSample/2*obj.pixelSize/1000*[1 -1];

            coordinates.yLimits = [min(fovPositions(:,2)) max(fovPositions(:,2))] + ...
                obj.parameters.segmentation.boundingBox(2) + ...
                [0 obj.parameters.segmentation.boundingBox(4)] + ...
                obj.parameters.display.downSample/2*obj.pixelSize/1000*[1 -1];
            
            % Save the coordinate system
            SaveAsByteStream([obj.normalizedDataPath obj.mosaicPath 'coordinates_slice_' num2str(s) '.matb'], coordinates);
            
            % Create tiff file to write
            tiffFileName = [obj.normalizedDataPath obj.mosaicPath 'slice_' num2str(s) '.tif'];
            tiffFile = Tiff(tiffFileName, 'w8');

            % Create tiff tags
            tiffTagStruct.ImageLength = numFinalPixelsX;
            tiffTagStruct.ImageWidth =numFinalPixelsY;
            tiffTagStruct.Photometric = Tiff.Photometric.MinIsBlack;
            tiffTagStruct.BitsPerSample = 16;
            tiffTagStruct.SamplesPerPixel = 1;
            tiffTagStruct.RowsPerStrip = 16;
            tiffTagStruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
            tiffTagStruct.Software = 'MATLAB';
            tiffTagStruct.ImageDescription = sprintf(['ImageJ=1.47a\n' ...      % ImageJ label
                'images=' num2str(obj.numDataChannels) '\n' ...            % The total number of images
                'channels=1\n' ...                                              % The number of channels, 1 for now
                'slices=' num2str(1) '\n' ...                                   % The number of z slices
                'frames=' num2str(obj.numDataChannels) '\n' ...            % The number of frames, i.e. data org entries
                'hyperstack=true\n' ...
                'loop=false\n' ...
                ]);
            
            % Loop over data channels
            for c=1:obj.numDataChannels
                % Reset frame (unnecessary)
                dsFrame = uint16(inf(numFinalPixelsX, numFinalPixelsY));

                % Loop over fields of view
                for v=1:length(localFovIDs)

                    % Read image from this file
                    localFileName = [obj.normalizedDataPath obj.warpedDataPath 'fov_' obj.fov2str(localFovIDs(v)) '.tif'];
                    localImage = imread(localFileName, obj.numZPos*(c-1) + zInd);

                    % Crop
                    localImage = localImage(startPixels(1):endPixels(1), ...
                        startPixels(2):endPixels(2)); % Crop

                    % Downsample
                    localImage = imresize(localImage, 1/obj.parameters.display.downSample);
                    
                    [dsNumPixelsX, dsNumPixelsY] = size(localImage);

                    % Invert/rotation (kludge)
                    localImage = localImage';

                    % Place into frame
                    xPixels = (1 + (xInd(v)-1)*dsNumPixelsX):xInd(v)*dsNumPixelsX;
                    yPixels = (1 + (yInd(v)-1)*dsNumPixelsY):yInd(v)*dsNumPixelsY;
                    dsFrame(xPixels, yPixels) = localImage;
                end

                % Write tiff frame
                tiffFile.setTag(tiffTagStruct);
                tiffFile.write(dsFrame'); % The matrix transpose is kludge
                if c~= obj.numDataChannels
                    tiffFile.writeDirectory(); % Write the directory for the next frame
                end

                % Display progress
                disp(['...completed channel ' num2str(c) ' of ' num2str(obj.numDataChannels)]);
                disp(['...in ' num2str(toc(localTimer)) ' s']);
                localTimer = tic;

            end
            % Close tiff file
            tiffFile.close();
        end
        % Display progress
        if obj.verbose
            disp(['...completed all low resolution mosaics in ' num2str(toc(totalTimer)) ' s']);
        end

    end
    % -------------------------------------------------------------------------
    % Sum raw fluorescence signals within individual boudnaries for each FOV
    % -------------------------------------------------------------------------
    function SumRawSignalFOV(obj,fovIDs)
        % Sum the raw signal from each channel within boundaries 
        %
        % obj.SumRawSignalFOV(fovIDs); % Analyze the specified fov 
        % obj.SumRawSignalFOV([]);     % Analyze all barcodes
        
        % -------------------------------------------------------------------------
        % Determine properties of the requested fov ids
        % -------------------------------------------------------------------------
        if isempty(fovIDs)
            fovIDs = obj.fovIDs;
        elseif ~all(ismember(fovIDs, obj.fovIDs))
            error('matlabFunctions:invalidArguments', 'An invalid fov id has been requested');
        end
        
        % -------------------------------------------------------------------------
        % Confirm analysis status before proceeding
        % -------------------------------------------------------------------------
        allWarped = obj.CheckStatus(fovIDs, 'warp');
        segmentationComplete = obj.CheckStatus(fovIDs, 'combine');
        if ~(all(allWarped) && all(segmentationComplete)) 
            error('matlabFunctions:missingData', 'Preprocessing of the data as well as segmentation must be complete prior to this analysis.');
        end
        
        % -------------------------------------------------------------------------
        % Make directories if they do not exist
        % -------------------------------------------------------------------------
        localSummationPath = [obj.normalizedDataPath obj.summationPath filesep];
        % Directory for barcodes by fov
        if ~exist(localSummationPath, 'dir')
            mkdir(localSummationPath);
        end
         
        % -------------------------------------------------------------------------
        % Load features
        % -------------------------------------------------------------------------
        % Check for existence
        foundFeaturesPath = [obj.normalizedDataPath obj.segmentationPath filesep 'final_found_features.matb'];
        if ~exist(foundFeaturesPath,'file')
            error('matlabFunctions:missingData', 'Final segmentation boundaries could not be found.');
        end
        
        % Load boundaries
        foundFeatures = LoadSplitByteStream(foundFeaturesPath, 'verbose', obj.verbose);
        
        % Copy (to preserve a direct link between parsed barcodes and the
        % boundaries
        if ~exist([localSummationPath 'final_found_features.matb'], 'file')
            % Load split bytestream info file to determine what files to
            % copy
            sbsInfo = LoadSplitByteStreamHeader(foundFeaturesPath);

            % Copy the spilt bytestream info file
            copyfile(foundFeaturesPath, [localSummationPath 'final_found_features.matb']);
            for C=1:sbsInfo.numBlocks
                copyfile([obj.normalizedDataPath obj.segmentationPath filesep sbsInfo.fileNames{C}], ...
                    [localSummationPath sbsInfo.fileNames{C}]);
            end
        end
        
        % Make local copy of parameters for segmentation
        parameters = obj.parameters.segmentation;

        % -------------------------------------------------------------------------
        % Determine which data channels and z slices to sum
        % -------------------------------------------------------------------------
        if isempty(obj.parameters.summation.dcIndsForSummation)
            dataChannelInds = 1:obj.numDataChannels;
        else
            dataChannelInds = obj.parameters.summation.dcIndsForSummation;
        end
        
        % Check validity of provided data channel inds
        if ~isempty(setdiff(dataChannelInds, 1:obj.numDataChannels))
            error('matlabFunctions:invalidArguments', 'The data channel indices provided do not match the known data channels');
        end
        
        % -------------------------------------------------------------------------
        % Prepare mask to remove pixels outside of the bounding box
        % -------------------------------------------------------------------------
        % Create pixel coordinate system
        [X,Y] = meshgrid(1:obj.imageSize(1), 1:obj.imageSize(2)); 

        % Create real-world scale coordinate system (but not yet
        % centered on fov position)
        x = obj.pixelSize/1000*obj.parameters.decoding.stageOrientation(1)*(X-obj.imageSize(1)/2); % x in microns
        y = obj.pixelSize/1000*obj.parameters.decoding.stageOrientation(2)*(Y-obj.imageSize(2)/2); % y in microns

        % Create ROI mask
        ROIMask = x >= parameters.boundingBox(1) & ...
            x <= (parameters.boundingBox(1) + parameters.boundingBox(3)) & ...
            y >= parameters.boundingBox(2) & ...
            y <= (parameters.boundingBox(2) + parameters.boundingBox(4));
                
        % -------------------------------------------------------------------------
        % Loop over FOV and sum raw signals
        % -------------------------------------------------------------------------
        spmd (obj.numPar)
            % Clear memory for accumulation registers
            totalSignal = zeros(length(dataChannelInds), max([foundFeatures.feature_id]));
            numberOfPixels = zeros(length(dataChannelInds), max([foundFeatures.feature_id]));
            
            % Loop over individual fov
            for f=labindex:numlabs:length(fovIDs)
                % Determine local fovID
                localFovID = fovIDs(f);
                
                % Create display strings
                if obj.verbose
                    displayStrings = {};
                    displayStrings{end+1} = PageBreak('nodisplay');
                    displayStrings{end+1} = ['Started summation of raw signal in fov ' obj.fov2str(localFovID) ' at ' datestr(now)];
                    displayStrings{end+1} = ['...extracting boundaries for this fov'];
                    fovTimer = tic;
                    localTimer = tic;
                    % Flush buffer
                    if numlabs==1        
                        display(char(displayStrings));
                        displayStrings = {};
                    end
                end
                
                % Identify the features the correspond to this fov
                isInFov = false(1, length(foundFeatures));
                for F=1:length(foundFeatures)
                    isInFov(F) = foundFeatures(F).InFov(localFovID);
                end
                
                % Crop these features
                localFeatures = foundFeatures(isInFov);
                                                
                % Create display strings
                if obj.verbose
                    displayStrings{end+1} = ['...searching ' num2str(length(localFeatures)) ' features']; 
                    displayStrings{end+1} = ['...completed in ' num2str(toc(localTimer)) ' s'];
                    localTimer = tic;
                    % Flush buffer
                    if numlabs==1        
                        display(char(displayStrings));
                        displayStrings = {};
                    end
                end
                
                % Define image data
                tiffName2Read = [obj.normalizedDataPath obj.warpedDataPath ...
                    'fov_' obj.fov2str(localFovID) '.tif'];
                if ~exist(tiffName2Read, 'file')
                    error('matlabFunctions:missingFile', 'The requsted tiff stack is not present.');
                end
                                              
                % Loop over the data channels and load z stacks
                for c=1:length(dataChannelInds)
                    
                    % Define frames to load
                    allZInds = 1:obj.numZPos;
                    possibleFrames = (dataChannelInds(c)-1)*obj.numZPos + allZInds;
                    
                    % Find the z positions in the data organization file
                    usedZPos = obj.dataOrganization(dataChannelInds(c)).zPos;
                    
                    % Determine the frames to keep
                    keptZInds = allZInds(ismember(obj.zPos, usedZPos));
                    framesToLoad = possibleFrames(ismember(obj.zPos, usedZPos));
                
                    % Display progress
                    if obj.verbose
                        displayStrings{end+1} = ['... loading stack for data channel ' num2str(dataChannelInds(c))]; 
                        loadTimer = tic;
                        % Flush buffer
                        if numlabs==1        
                            display(char(displayStrings));
                            displayStrings = {};
                        end
                    end

                    % Allocate memory for the necessary frames
                    localImageStack = zeros([obj.imageSize length(framesToLoad)], 'uint16');

                    % Load these frames
                    for z=1:length(framesToLoad)
                        % Load frames for watershed and seed generation
                        localImageStack(:,:,z) = imread(tiffName2Read, framesToLoad(z));
                    end
                    
                    % Display progress
                    if obj.verbose
                        displayStrings{end+1} = ['...completed in ' num2str(toc(loadTimer)) ' s'];
                        displayStrings{end+1} = ['...parsing ' num2str(length(localFeatures)) ' features'];
                        localTimer = tic;
                        totalTimer = tic;
                        % Flush buffer
                        if numlabs==1        
                            display(char(displayStrings));
                            displayStrings = {};
                        end
                    end
                    
                    % Loop over local features and add
                    for F=1:length(localFeatures)
                        % Create feature mask
                        mask = localFeatures(F).GeneratePixelMask(localFovID, keptZInds);
                        
                        % Crop to ROI
                        mask = mask & repmat(ROIMask, [1 1 size(mask,3)]);
                        
                        % Add signal
                        totalSignal(c,localFeatures(F).feature_id) = totalSignal(c,localFeatures(F).feature_id) + ...
                            sum(localImageStack(mask(:)));
                        
                        % Add number of pixels
                        numberOfPixels(c,localFeatures(F).feature_id) = numberOfPixels(c,localFeatures(F).feature_id) + ...
                            sum(mask(:));
                        
                        % Display progress
                        if obj.verbose && ~mod(F, 100)
                            displayStrings{end+1} = ['...completed ' num2str(F) ' features in ' num2str(toc(localTimer)) ' s'];
                            localTimer = tic;
                            % Flush buffer
                            if numlabs==1        
                                display(char(displayStrings));
                                displayStrings = {};
                            end
                        end

                    end
                    
                    % Display progress
                    if obj.verbose 
                        displayStrings{end+1} = ['...completed all features in ' num2str(toc(totalTimer)) ' s'];
                        % Flush buffer
                        if numlabs==1        
                            display(char(displayStrings));
                            displayStrings = {};
                        end
                    end
                end
                                    
                % Save total signal
                totalSignalFilePath = [localSummationPath ...
                            'total_signal_fov_'  obj.fov2str(localFovID) '.matb'];
                SaveAsByteStream(totalSignalFilePath, totalSignal);
                
                % Save number of pixels
                numberOfPixelsFilePath = [localSummationPath ...
                            'total_pixels_fov_'  obj.fov2str(localFovID) '.matb'];
                SaveAsByteStream(numberOfPixelsFilePath, numberOfPixels);

                % Flush display buffer
                if obj.verbose
                    displayStrings{end+1} = ['...saved ' totalSignalFilePath];
                    displayStrings{end+1} = ['...saved ' numberOfPixelsFilePath];
                    displayStrings{end+1} = ['...completed fov ' obj.fov2str(localFovID) ' in ' num2str(toc(fovTimer)) ' s'];
                    display(char(displayStrings));
                end
                
            end % End loop over fov
        end % End spmd loop        
    end % End function
    
    % -------------------------------------------------------------------------
    % Combine raw summation of signal
    % -------------------------------------------------------------------------
    function CombineRawSum(obj)
        % Combine the sum of raw signals within features for each fov
        %
        % obj.CombineRawSum();
        
        % Display progress
        if obj.verbose
            PageBreak();
            display('Combining raw signal summation for boundaries in individual fov');
        end
        
        % Create paths to combined final values
        combinedSignalFilePath = [obj.normalizedDataPath obj.summationPath ...
                'total_signal.matb'];

        combinedPixelsFilePath = [obj.normalizedDataPath obj.summationPath ...
                'total_pixels.matb'];
        
        % Handle existing analysis with overwrite off
        if exist(combinedSignalFilePath, 'file') && exist(combinedPixelsFilePath, 'file')
            display('...found existing combined sum. Skipping analysis.');
            return;
        end
        
        % Map the existing signal and pixel files
        foundSignalFiles = BuildFileStructure([obj.normalizedDataPath obj.summationPath], ...
            'fileExt', 'matb', ...
            'regExp', 'total_signal_fov_(?<fov>[0-9]+)', ...
            'fieldNames', {'fov'}, ...
            'fieldConv', {@str2num});
        foundPixelsFiles = BuildFileStructure([obj.normalizedDataPath obj.summationPath], ...
            'fileExt', 'matb', ...
            'regExp', 'total_pixels_fov_(?<fov>[0-9]+)', ...
            'fieldNames', {'fov'}, ...
            'fieldConv', {@str2num});
        
        % Display progress
        if obj.verbose
            display(['...found ' num2str(length(foundSignalFiles)) ' signal files']);
            display(['...found ' num2str(length(foundPixelsFiles)) ' pixel files']);
        end

        % Check for missing fov
        missingFovIDs = union(setdiff(obj.fovIDs, [foundSignalFiles.fov]), ...
            setdiff(obj.fovIDs, [foundPixelsFiles.fov]));
        if ~isempty(missingFovIDs)
            display(['...the following fov are missing some files associated with raw signal summation'])
            display(['... ' num2str(missingFovIDs)]);
            error('matlabFunctions:missingData', 'Combination cannot proceed until all fov have been summed');
        end
        
        % Load and sum signals and piuxels  
        if obj.verbose
            display(['...loading total sum and pixel numbers']);
            localTimer = tic;
        end
        totalSignal = [];
        totalPixels = [];
        for f=1:obj.numFov
            % Handle the case of the first file load
            if f==1
                totalSignal = LoadByteStream(foundSignalFiles([foundSignalFiles.fov] == obj.fovIDs(f)).filePath, 'verbose', false);
                totalPixels = LoadByteStream(foundPixelsFiles([foundPixelsFiles.fov] == obj.fovIDs(f)).filePath, 'verbose', false);
            else
                totalSignal = totalSignal + LoadByteStream(foundSignalFiles([foundSignalFiles.fov] == obj.fovIDs(f)).filePath, 'verbose', false);
                totalPixels = totalPixels + LoadByteStream(foundPixelsFiles([foundPixelsFiles.fov] == obj.fovIDs(f)).filePath, 'verbose', false);
            end
        end
                
        % Save combined signals as csv
        combinedSignalFilePath = [obj.normalizedDataPath obj.summationPath ...
            'total_signal.csv'];
        combinedPixelsFilePath = [obj.normalizedDataPath obj.summationPath ...
            'total_pixels.csv'];
        csvwrite(combinedSignalFilePath, totalSignal);
        csvwrite(combinedPixelsFilePath, totalPixels);

        % Define the feature name file
        featuresNameFilePath = [obj.normalizedDataPath obj.summationPath 'featureNames.csv'];
        
        % Load features
        foundFeatures = obj.GetFoundFeatures();
        
        % Generate the contents of the feature names file
        featureUIDs = {foundFeatures.uID}; % Feature ids
        
        % Open and write file
        fid = fopen(featuresNameFilePath, 'W');
        fprintf(fid, '%s\n', featureUIDs{:});
        fclose(fid);

        % Reconstitute and return the channel names
        channelNamesFilePath = [obj.normalizedDataPath obj.summationPath 'channelNames.csv'];
        if isempty(obj.parameters.summation.dcIndsForSummation)
            indsForSum = 1:obj.numDataChannels;
        else
            indsForSum = obj.parameters.summation.dcIndsForSummation;
        end
        channelNames = {obj.dataOrganization(indsForSum).bitName};
                
        % Open and write file
        fid = fopen(channelNamesFilePath, 'W');
        fprintf(fid, '%s\n', channelNames{:});
        fclose(fid);
        
        % Delete intermediate files
        for f=1:length(foundSignalFiles)
            delete(foundSignalFiles(f).filePath);
        end
        for f=1:length(foundPixelsFiles)
            delete(foundPixelsFiles(f).filePath);
        end        
        
        % Display progress
        if obj.verbose
            display(['...completed in ' num2str(toc(localTimer)) ' s']);
        end

    end
    
    % -------------------------------------------------------------------------
    % Get counts per feature
    % -------------------------------------------------------------------------
    function [exactCountsIn, correctedCountsIn, exactCountsOut, correctedCountsOut, ...
            geneNames, featureUIDs] = GetCountsPerFeature(obj)
        % Return the counts per feature
        % 
        % [exactCountsIn, correctedCountsIn, exactCountsOut, correctedCountsOut, geneNames, featureUIDs] = obj.GetCountsPerFeature()

        % Check to see if the feature counts have been created
        if ~obj.CheckStatus(obj.fovIDs(1), 'n')
            error('matlabFunctions:incompleteAnalysis', 'The feature counts have not yet been computed.');
        end
        
        % Display progress
        if obj.verbose
            PageBreak();
            disp('Loading feature counts');
            localTimer = tic;
        end
        
        % Load the counts
        reportsPath = [obj.normalizedDataPath 'reports' filesep];

        exactCountsIn = csvread([reportsPath 'countsPerCellExactIn.csv']);
        correctedCountsIn = csvread([reportsPath 'countsPerCellCorrectedIn.csv']);
        exactCountsOut = csvread([reportsPath 'countsPerCellExactOut.csv']); 
        correctedCountsOut = csvread([reportsPath 'countsPerCellCorrectedOut.csv']);
        
        % Load the feature unique ids
        fid = fopen([reportsPath 'featureNames.csv'], 'r');
        line = fgetl(fid);
        fclose(fid);
        featureUIDs = strsplit(line, ',');
        
        % Load the gene names
        fid = fopen([reportsPath 'geneNames.csv'], 'r');
        line = fgetl(fid);
        fclose(fid);
        geneNames = strsplit(line, ',');

        % Display progress
        if obj.verbose
            disp(['...completed in ' num2str(toc(localTimer)) ' s']);
        end

    end
    
    % -------------------------------------------------------------------------
    % Get fov absolute coordinates
    % -------------------------------------------------------------------------
    function coordinates = GetFOVCoordinates(obj, fovID)
        % Return the absolute coordinates of the specified fov
        % 
        % coordinates = obj.GetFOVCoordinates(fovID)
        % coordinates = [xstart xend ystart yend] where xstart is the
        % center of the first pixel and xend is the center of the final
        % pixel in x

        % Check the input
        if nargin < 1
            error('matlabFunctions:invalidArguments', 'A fov id must be specified');
        end
        if ~ismember(fovID, obj.fovIDs)
            error('matlabFunctions:invalidArguments', 'The specified fov id is not valid');
        end
        
        % Compute the coordinates
        fovPos = obj.fovPos(obj.fovIDs == fovID,:);
        
        coordinates(1:2) = [-1 1]*obj.imageSize(1)*(obj.pixelSize/1000)/2 + fovPos(1);
        coordinates(3:4) = [-1 1]*obj.imageSize(2)*(obj.pixelSize/1000)/2 + fovPos(2);

    end
        
    
    % -------------------------------------------------------------------------
    % Return a low resolution mosaic of the sample
    % -------------------------------------------------------------------------
    function [mosaicImageStack, coordinates] = GetMosaic(obj, sliceID, framesToLoad)
        % Return the low resolution mosaic and the absolute coordinate
        % system
        % 
        % mosaicImageStack = obj.GetMosaic(sliceID)
        
        % Check to see if the mosaics have been created
        if ~obj.CheckStatus(obj.fovIDs(1), 'l')
            error('matlabFunctions:incompleteAnalysis', 'The mosaics have not yet been created.');
        end
        
        % Check to see if the requested sliceID is valid
        if sliceID > length(obj.sliceIDs) || sliceID < 1
            error('matlabFunctions:invalidArgument', 'The provided sliceID is not valid.');
        end

        % Determine mosaic stack info
        stackInfo = imfinfo([obj.normalizedDataPath obj.mosaicPath 'slice_' num2str(sliceID) '.tif']);
                
        % Determine the frames to load if not provided
        if ~exist('framesToLoad', 'var')
            framesToLoad = 1:length(stackInfo);
        end
        
        % Check the provided frames to load
        if any(~ismember(framesToLoad, 1:length(stackInfo)))
            error('matlabFunctions:invalidArguments', 'The provided frame ids to load are not valid.');
        end
        
        % Allocate memory
        mosaicImageStack = zeros(stackInfo(1).Width, stackInfo(1).Height, ...
            length(framesToLoad), 'uint16');

        % Load the requested frames
        for s=1:length(framesToLoad)
            mosaicImageStack(:,:,s) = imread([obj.normalizedDataPath obj.mosaicPath 'slice_' num2str(sliceID) '.tif'],framesToLoad(s));
        end
        
        % Load the coordinate system
        coordinates = LoadByteStream([obj.normalizedDataPath obj.mosaicPath 'coordinates_slice_' num2str(sliceID) '.matb']);
        
    end    
    
    % -------------------------------------------------------------------------
    % Return the combined signal data
    % -------------------------------------------------------------------------
    function [normalizedSignal, sumSignal, sumPixels, channelNames, featureUIDs] = GetSummedSignal(obj)
        % Return the output of the signal summation process
        % [normalizedSignal, sumSignal, sumPixels, channelNames, featureUIDs] = obj.GetSummedSignal
        
        % Check to see if the sum signals have been combined
        if ~all(obj.CheckStatus([], 'i'))
            error('matlabFunctions:missingData', 'The summation has not yet been completed.');
        end
        
        % Define paths to the different objects
        combinedSignalFilePath = [obj.normalizedDataPath obj.summationPath ...
            'total_signal.csv'];
        combinedPixelsFilePath = [obj.normalizedDataPath obj.summationPath ...
            'total_pixels.csv'];

        % Load the sum objects
        sumSignal = csvread(combinedSignalFilePath);
        sumPixels = csvread(combinedPixelsFilePath);

        % Compute the normalized sum
        % Normalized the sum signal by the feature volume
        normalizedSignal = sumSignal./sumPixels;
        
        % Load the featureUIDs
        featureUIDs = {};
        fid = fopen([obj.normalizedDataPath obj.summationPath 'featureNames.csv'], 'r');
        line = fgetl(fid);
        while line ~= -1
            featureUIDs{end+1} = line;
            line = fgetl(fid);
        end
        fclose(fid);
        
        % Load the channel names
        channelNames = {};
        fid = fopen([obj.normalizedDataPath obj.summationPath 'channelNames.csv'], 'r');
        line = fgetl(fid);
        while line ~= -1
            channelNames{end+1} = line;
            line = fgetl(fid);
        end
        fclose(fid);
        
    end
    
    % -------------------------------------------------------------------------
    % Parse barcodes into final feature boundaries for individual fov
    % -------------------------------------------------------------------------
    function ParseFOV(obj,fovIDs)
        % Parse decoded barcodes into feature boundaries within individual
        % fov 
        %
        % obj.ParseFOV(fovIDs); % Parse the specified fov barcodes
        % obj.ParseFOV([]); % Parse all barcodes
        
        % -------------------------------------------------------------------------
        % Determine properties of the requested fov ids
        % -------------------------------------------------------------------------
        if isempty(fovIDs)
            fovIDs = obj.fovIDs;
        elseif ~all(ismember(fovIDs, obj.fovIDs))
            error('matlabFunctions:invalidArguments', 'An invalid fov id has been requested');
        end
        
        % -------------------------------------------------------------------------
        % Make directories if they do not exist
        % -------------------------------------------------------------------------
        barcodeByFovPath = [obj.normalizedDataPath obj.barcodePath filesep 'barcode_fov' filesep];
        % Check for existing barcodes
        if ~exist(barcodeByFovPath, 'dir')
            error('matlabFunctions:missingData', 'No barcodes could be found.');
        end
        
        % Make directory for the parsed barcodes
        parsedBarcodePath = [obj.normalizedDataPath obj.barcodePath filesep 'parsed_fov' filesep];
        % Directory for barcodes by fov
        if ~exist(parsedBarcodePath, 'dir')
            mkdir(parsedBarcodePath);
        end
         
        % -------------------------------------------------------------------------
        % Load feature boundaries
        % -------------------------------------------------------------------------
        % Check for existence
        foundFeaturesPath = [obj.normalizedDataPath obj.segmentationPath filesep 'final_found_features.matb'];
        if ~exist(foundFeaturesPath,'file')
            error('matlabFunctions:missingData', 'Final segmentation boundaries could not be found.');
        end
        
        % Load boundaries
        foundFeatures = LoadSplitByteStream(foundFeaturesPath, 'verbose', obj.verbose);
        
        % Copy (to preserve a direct link between parsed barcodes and the
        % boundaries
        % Load split bytestream info file to determine what files to
        % copy
        sbsInfo = LoadSplitByteStreamHeader(foundFeaturesPath);

        % Copy the spilt bytestream info file
        copyfile(foundFeaturesPath, [parsedBarcodePath 'final_found_features.matb']);
        for C=1:sbsInfo.numBlocks
            copyfile([obj.normalizedDataPath obj.segmentationPath filesep sbsInfo.fileNames{C}], ...
                [parsedBarcodePath sbsInfo.fileNames{C}]);
        end
        
        % -------------------------------------------------------------------------
        % Prepare copies of variables for parsing loop
        % -------------------------------------------------------------------------
        % Make local copy of parameters for segmentation
        parameters = obj.parameters.segmentation;

        % -------------------------------------------------------------------------
        % Parse via individual fov
        % -------------------------------------------------------------------------
        spmd (obj.numPar)
            % Loop over individual fov
            for f=labindex:numlabs:length(fovIDs)
                % Determine local fovID
                localFovID = fovIDs(f);
                
                % Create display strings
                if obj.verbose
                    displayStrings = {};
                    displayStrings{end+1} = PageBreak('nodisplay');
                    displayStrings{end+1} = ['Started parsing of barcodes in fov ' num2str(localFovID) ' at ' datestr(now)];
                    displayStrings{end+1} = ['...extracting boundaries for this fov'];
                    fovTimer = tic;
                    localTimer = tic;
                    % Flush buffer
                    if numlabs==1        
                        display(char(displayStrings));
                        displayStrings = {};
                    end
                end
                
                % Find all possible adjacent fov (to cut down on the
                % boundaries that need to be parsed)
                nnIDX = knnsearch(obj.fovPos, obj.fovPos(obj.fovIDs==localFovID,:), ...
                    'K', 9); % 8 neighbors + 1 duplicate of self
                fovIDsToSearch = obj.fovIDs(nnIDX);
                
                % Reduce the finalBoundaries
                goodFeatureInds = false(1, length(foundFeatures));
                for F=1:length(foundFeatures)
                    goodFeatureInds(F) = foundFeatures(F).InFov(fovIDsToSearch);
                end
                localFeatures = foundFeatures(goodFeatureInds);
                                
                % Create display strings
                if obj.verbose
                    displayStrings{end+1} = ['...searching ' num2str(length(localFeatures)) ' boundaries']; 
                    displayStrings{end+1} = ['...completed in ' num2str(toc(localTimer)) ' s'];
                    displayStrings{end+1} = ['...loading barcodes'];
                    localTimer = tic;
                    % Flush buffer
                    if numlabs==1        
                        display(char(displayStrings));
                        displayStrings = {};
                    end
                end
                
				
                % Check for existing or corrupt barcode list
                barcodeFilePath = [barcodeByFovPath 'fov_' obj.fov2str(localFovID) '_blist.bin'];
                if ~exist(barcodeFilePath, 'file') 
                    error('matlabFunctions:missingData', 'Could not find specified barcode file.');
                end
                try 
                    header = ReadBinaryFileHeader(barcodeFilePath);
                    isCorrupt = header.isCorrupt;
                catch
                    isCorrupt = true;
                end
                if isCorrupt
                    error('matlabFunctions:missingData', 'The requested barcode file appears to be corrupt.');
                end
                    
                % Load bList
                bList = ReadBinaryFile(barcodeFilePath);
                
                % Display progress
                if obj.verbose
                    displayStrings{end+1} = ['...loaded ' num2str(length(bList)) ' barcodes'];
                    displayStrings{end+1} = ['...completed load in ' num2str(toc(localTimer)) ' s'];
                    displayStrings{end+1} = ['...trimming barcodes to bounding box'];
                    localTimer = tic;
                    % Flush buffer
                    if numlabs==1        
                        display(char(displayStrings));
                        displayStrings = {};
                    end
                end
                
                % Extract positions of loaded barcodes
                pos = cat(1, bList.abs_position);
                
                % Define bounding box (map to absolute position in the
                % sample)
                boundingBox = parameters.boundingBox; % Bounding box used to cut segmentation boundaries
                boundingBox(1:2) = boundingBox(1:2) + obj.fovPos(obj.fovIDs==localFovID,:); % Map to real world coordinates
                
                % Define positions within this bounding box
                indsToKeep = pos(:,1) >= boundingBox(1) & ...
                    pos(:,1) <= (boundingBox(1) + boundingBox(3)) & ...
                    pos(:,2) >= boundingBox(2) & ...
                    pos(:,2) <= (boundingBox(2) + boundingBox(4));
                                
                % Slice list to remove barcodes not inside the bounding box
                bList = bList(indsToKeep);
                
                % Display progress
                if obj.verbose
                    displayStrings{end+1} = ['...cut ' num2str(sum(~indsToKeep)) ' barcodes outside of segmentation bounding box'];
                    displayStrings{end+1} = ['...completed in ' num2str(toc(localTimer)) ' s'];
                    displayStrings{end+1} = ['...parsing barcodes'];
                    localTimer = tic;
                    % Flush buffer
                    if numlabs==1        
                        display(char(displayStrings));
                        displayStrings = {};
                    end
                end
                				
                % Add additional fields to barcode list
                [bList(:).feature_id] = deal(int32(-1));            % The id of the feature to which the barcode was assigned 
                [bList(:).feature_dist] = deal(zeros(1, 'single')); % The distance to the nearest feature edge
                [bList(:).in_feature] = deal(zeros(1, 'uint8'));    % A boolean indicating whether or not (true/false) the RNA falls within the feature to which it was associated
                
                % Prepare a list of all barcode positions
                barcodePos = cat(1, bList.abs_position);

                %Discretize z positions for all barcodes
                zInds = discretize(barcodePos(:,3), [obj.zPos Inf]);
                
				% Handle the case that there are no barcodes or no features
				if ~isempty(bList) & ~isempty(localFeatures)
					%Loop over z indices
					for z=1:obj.numZPos
						%Compile a composite list of boundaries in this z
						%plane
						combinedBoundaries = zeros(0,2);
						isDilatedBoundary = zeros(0,1);
						localFeatureIDs = zeros(0,1);
						for F=1:length(localFeatures)
							% Concatenate boundaries and dilated boundaries
							combinedBoundaries = cat(1, combinedBoundaries, ...
								localFeatures(F).abs_boundaries{z},...
								localFeatures(F).DilateBoundary(z, obj.pixelSize/1000*parameters.dilationSize));
							
							% Concatenate indices for features
							localFeatureIDs = cat(1, localFeatureIDs, ...
								F*ones(length(localFeatures(F).abs_boundaries{z}),1), ...
								F*ones(length(localFeatures(F).abs_boundaries{z}),1));
							
							% Concatenate flags for inside or outside feature
							isDilatedBoundary = cat(1, isDilatedBoundary, ...
								false(length(localFeatures(F).abs_boundaries{z}),1), ...
								true(length(localFeatures(F).abs_boundaries{z}),1));
								
                        end
						
                        %Handle the case of no boundaries in the z-plane
                        if isempty(combinedBoundaries)
                            continue;
                        end
                        
						%Find the barcode indices in this zPlane
						localBarcodeIndex = find(zInds == z);
						
						%Find the nearest neighbor point in the boundaries and
						%the distance
						[nnIDX, D] = knnsearch(combinedBoundaries, barcodePos(localBarcodeIndex,1:2));
						
						%Loop through barcodes to assign found values and
						%determine if they are in the appropriate features
						for b=1:length(localBarcodeIndex)
							bList(localBarcodeIndex(b)).feature_id = int32(localFeatures(localFeatureIDs(nnIDX(b))).feature_id);
							bList(localBarcodeIndex(b)).feature_dist = single(D(b));
							bList(localBarcodeIndex(b)).in_feature = uint8(~isDilatedBoundary(nnIDX(b)));
						end
					end
				else
					% Display progress
					if obj.verbose
						displayStrings{end+1} = ['...handling case of no features or no barcodes...'];
					end
				end
                                
                % Display progress
                if obj.verbose
                    displayStrings{end+1} = ['...completed in ' num2str(toc(localTimer)) ' s'];
                    displayStrings{end+1} = ['...writing barcode list'];

                    localTimer = tic;
                    % Flush buffer
                    if numlabs==1        
                        display(char(displayStrings));
                        displayStrings = {};
                    end
                end
                
                % Save the barcode list
                % Write binary file for all measured barcodes
                barcodeFile = [parsedBarcodePath 'fov_' obj.fov2str(localFovID) '_blist.bin'];
                WriteBinaryFile(barcodeFile, bList);

                % Display progress
                if obj.verbose
                    displayStrings{end+1} = ['...wrote ' barcodeFile];
                    displayStrings{end+1} = ['...complete in ' num2str(toc(localTimer)) ' s'];
                    displayStrings{end+1} = ['Completed analysis of fov ' obj.fov2str(localFovID) ' at ' datestr(now) ' in ' num2str(toc(fovTimer)) ' s'];
                    display(char(displayStrings)); % Final buffer flush
                end
            end % End loop over fov
        end % End spmd loops
        
    end % End function
        
    % -------------------------------------------------------------------------
    % Find individual molecules within a fov
    % -------------------------------------------------------------------------
    function FindMoleculesFOV(obj,fovIDs)
        % Find individual molecules (via image regional max), assign to
        % individual cells 
        %
        % obj.FindMoleculesFOV(fovIDs); % Find molecules, assign to cells,
        % for all z stacks in all data channels in the specified fov
        % obj.FindMoleculesFOV([]); % Run this analysis for all fov

        % -------------------------------------------------------------------------
        % Determine properties of the requested fov ids
        % -------------------------------------------------------------------------
        if isempty(fovIDs)
            fovIDs = obj.fovIDs;
        elseif ~all(ismember(fovIDs, obj.fovIDs))
            error('matlabFunctions:invalidArguments', 'An invalid fov id has been requested');
        end
                
        % -------------------------------------------------------------------------
        % Make directories if they do not exist
        % -------------------------------------------------------------------------
        moleculeByFovPath = [obj.normalizedDataPath filesep 'smFISH' filesep];
        % Check for existing barcodes
        if ~exist(moleculeByFovPath, 'dir')
            mkdir(moleculeByFovPath);
        end
                         
        % -------------------------------------------------------------------------
        % Prepare copies of variables for parsing loop
        % -------------------------------------------------------------------------
        % Make local copy of parameters for segmentation
        parameters = obj.parameters.molecules;
        boundingBox = obj.parameters.segmentation.boundingBox;
        
        % Handle the case that no information was provided on the
        % dataChannels and zstacks
        if isempty(parameters.molDataChannels)
            dataChannelDescription = cell(0,2);
            for i=1:size(obj.dataOrganization,1)
                dataChannelDescription(end+1,:) = {obj.dataOrganization(i).bitName, 1:obj.numZPos};
            end
        else
            dataChannelDescription = parameters.molDataChannels;
        end
        
        % -------------------------------------------------------------------------
        % Load feature boundaries
        % -------------------------------------------------------------------------
        foundFeatures = obj.GetFoundFeatures();

        % -------------------------------------------------------------------------
        % Find spots within individual fov
        % -------------------------------------------------------------------------
        spmd (obj.numPar)
            % Loop over individual fov
            for f=labindex:numlabs:length(fovIDs)
                % Determine local fovID
                localFovID = fovIDs(f);
                
                % Loop over the requested data Channels to compile a list
                % of molecules in this fov for all data channels
                mList = [];
                
                for D = 1:size(dataChannelDescription,1)

                    % Create display strings
                    if obj.verbose
                        displayStrings = {};
                        displayStrings{end+1} = PageBreak('nodisplay');
                        displayStrings{end+1} = ['Identifying molecules in channel ' dataChannelDescription{D,1} ' for fov ' num2str(localFovID) ' at ' datestr(now)];
                        fovTimer = tic;
                        localTimer = tic;
                        % Flush buffer
                        if numlabs==1        
                            display(char(displayStrings));
                            displayStrings = {};
                        end
                    end
                                
                    % Extract the requested frame ids
                    frameIDs = dataChannelDescription{D,2};
                    
                    % Loop over the requested image frames
                    for z=1:length(frameIDs)
                        
                        % Load the requested image frame
                        frame = obj.GetImage(localFovID, dataChannelDescription{D,1}, frameIDs(z));
                        
                        % Low pass filter the image to remove background
                        frame = frame - imgaussfilt(frame, parameters.molLowPassfilterSize);
                        
                        % Compute the region max
                        rMax = imregionalmax(frame);
                        
                        % Remove values that are too dim
                        rMax(frame < parameters.molIntensityThreshold) = 0;
                        
                        % Combine touching pixels
                        props = regionprops(rMax, 'Centroid');
                        
                        % Create a combined list of centroids
                        centroids = cat(1, props.Centroid);
                        
                        % Add the z index
                        centroids(:,3) = frameIDs(z);
                        
                        % Convert to real coordinates
                        abs_position = obj.Pixel2Abs(centroids, localFovID);
                        
                        % Drop molecules that fall outside of the
                        % segmentation bounding box
                        shiftedPos = abs_position(:,1:2) - repmat(obj.fovPos(obj.fovIDs == localFovID,:), [size(abs_position,1) 1]);
                                                
                        moleculesToKeep = shiftedPos(:,1) >= boundingBox(1) & shiftedPos(:,2) >= boundingBox(2) & ...
                            shiftedPos(:,1) <= (boundingBox(1)+boundingBox(3)) & ...
                            shiftedPos(:,2) <= (boundingBox(2)+boundingBox(4));
                        
                        % Cut the molecules
                        centroids = centroids(moleculesToKeep,:);
                        abs_position = abs_position(moleculesToKeep,:);
                        
                        % Calculate the brightness: First identify pixels
                        % to integrate
                        xBounds = repmat(round(centroids(:,1)), [1 2]) + repmat([-parameters.molNumPixelSum parameters.molNumPixelSum], ...
                            [size(centroids,1) 1]);
                        yBounds = repmat(round(centroids(:,2)), [1 2]) + repmat([-parameters.molNumPixelSum parameters.molNumPixelSum], ...
                            [size(centroids,1) 1]);
                        
                        brightness = zeros(1, size(xBounds,1));
                        for M=1:size(xBounds,1)
                            brightness(M) = sum(sum(frame(yBounds(M,1):yBounds(M,2), xBounds(M,1):xBounds(M,2))));
                        end
                        
                        % Create molecule list structure
                        newMolecules = struct('centroid', num2cell(centroids,2), 'abs_position', num2cell(abs_position,2), ...
                            'brightness', num2cell(brightness)', ...
                            'channel', find(strcmp({obj.dataOrganization.bitName},  dataChannelDescription{D,1})));
                            
                        % Append these molecules to the list of all
                        % molecules
                        mList = cat(1, mList, newMolecules);

                    end
                    
                end
                
                % Create display strings
                if obj.verbose
                    displayStrings{end+1} = ['...completed in ' num2str(toc(localTimer)) ' s'];
                    localTimer = tic;
                    % Flush buffer
                    if numlabs==1        
                        display(char(displayStrings));
                        displayStrings = {};
                    end
                end

                % Create display strings
                if obj.verbose
                    displayStrings = {};
                    displayStrings{end+1} = PageBreak('nodisplay');
                    displayStrings{end+1} = ['Parsing molecules in fov ' num2str(localFovID) ' at ' datestr(now)];
                    fovTimer = tic;
                    localTimer = tic;
                    % Flush buffer
                    if numlabs==1        
                        display(char(displayStrings));
                        displayStrings = {};
                    end
                end

                
                % Find all possible adjacent fov (to cut down on the
                % boundaries that need to be parsed)
                nnIDX = knnsearch(obj.fovPos, obj.fovPos(obj.fovIDs==localFovID,:), ...
                    'K', 9); % 8 neighbors + 1 duplicate of self
                fovIDsToSearch = obj.fovIDs(nnIDX);
                
                % Reduce the finalBoundaries
                goodFeatureInds = false(1, length(foundFeatures));
                for F=1:length(foundFeatures)
                    goodFeatureInds(F) = foundFeatures(F).InFov(fovIDsToSearch);
                end
                localFeatures = foundFeatures(goodFeatureInds);
                       
                
                % Handle the case that no molecules were found
                if isempty(mList)
                    % Create display strings
                    if obj.verbose
                        displayStrings{end+1} = ['...searching ' num2str(length(localFeatures)) ' boundaries']; 
                        displayStrings{end+1} = ['...complete in ' num2str(toc(localTimer)) ' s'];
                        displayStrings{end+1} = 'No molecules found. Skipping the parsing.';
                        displayStrings{end+1} = ['Completed analysis of fov ' obj.fov2str(localFovID) ' at ' datestr(now) ' in ' num2str(toc(fovTimer)) ' s'];
                        display(char(displayStrings)); % Final buffer flush
                    end
                    
                    % Continue
                    continue;
                end

                
                % Create display strings
                if obj.verbose
                    displayStrings{end+1} = ['...searching ' num2str(length(localFeatures)) ' boundaries']; 
                    displayStrings{end+1} = ['...completed in ' num2str(toc(localTimer)) ' s'];
                    displayStrings{end+1} = ['...parsing ' num2str(length(mList)) ' molecules'];
                    localTimer = tic;
                    % Flush buffer
                    if numlabs==1        
                        display(char(displayStrings));
                        displayStrings = {};
                    end
                end
                                                
                % Add additional fields to barcode list
                [mList(:).feature_id] = deal(int32(-1));            % The id of the feature to which the barcode was assigned 
                [mList(:).feature_dist] = deal(zeros(1, 'single')); % The distance to the nearest feature edge
                [mList(:).in_feature] = deal(zeros(1, 'uint8'));    % A boolean indicating whether or not (true/false) the RNA falls within the feature to which it was associated
                
                % Prepare a list of all barcode positions
                barcodePos = cat(1, mList.abs_position);

                %Discretize z positions for all barcodes
                zInds = discretize(barcodePos(:,3), [obj.zPos Inf]);
                
                %Loop over z indices
                for z=1:obj.numZPos
                    %Compile a composite list of boundaries in this z
                    %plane
                    combinedBoundaries = zeros(0,2);
                    isDilatedBoundary = zeros(0,1);
                    localFeatureIDs = zeros(0,1);
                    for F=1:length(localFeatures)
                        % Concatenate boundaries and dilated boundaries
                        combinedBoundaries = cat(1, combinedBoundaries, ...
                            localFeatures(F).abs_boundaries{z},...
                            localFeatures(F).DilateBoundary(z, obj.pixelSize/1000*obj.parameters.segmentation.dilationSize));
                        
                        % Concatenate indices for features
                        localFeatureIDs = cat(1, localFeatureIDs, ...
                            F*ones(length(localFeatures(F).abs_boundaries{z}),1), ...
                            F*ones(length(localFeatures(F).abs_boundaries{z}),1));
                        
                        % Concatenate flags for inside or outside feature
                        isDilatedBoundary = cat(1, isDilatedBoundary, ...
                            false(length(localFeatures(F).abs_boundaries{z}),1), ...
                            true(length(localFeatures(F).abs_boundaries{z}),1));
                            
                    end
                    
                    %Find the barcode indices in this zPlane
                    localBarcodeIndex = find(zInds == z);
                    
                    %Find the nearest neighbor point in the boundaries and
                    %the distance
                    [nnIDX, D] = knnsearch(combinedBoundaries, barcodePos(localBarcodeIndex,1:2));
                    
                    %Loop through barcodes to assign found values and
                    %determine if they are in the appropriate features
                    for b=1:length(localBarcodeIndex)
                        mList(localBarcodeIndex(b)).feature_id = int32(localFeatures(localFeatureIDs(nnIDX(b))).feature_id);
                        mList(localBarcodeIndex(b)).feature_dist = single(D(b));
                        mList(localBarcodeIndex(b)).in_feature = uint8(~isDilatedBoundary(nnIDX(b)));
                    end
                end
                                
                % Display progress
                if obj.verbose
                    displayStrings{end+1} = ['...completed in ' num2str(toc(localTimer)) ' s'];
                    displayStrings{end+1} = ['...writing molecule list'];

                    localTimer = tic;
                    % Flush buffer
                    if numlabs==1        
                        display(char(displayStrings));
                        displayStrings = {};
                    end
                end
                
                % Save the barcode list
                % Write binary file for all measured barcodes
                moleculeFile = [moleculeByFovPath 'fov_' obj.fov2str(localFovID) '_mlist.bin'];
                WriteBinaryFile(moleculeFile, mList);

                % Display progress
                if obj.verbose
                    displayStrings{end+1} = ['...wrote ' moleculeFile];
                    displayStrings{end+1} = ['...complete in ' num2str(toc(localTimer)) ' s'];
                    displayStrings{end+1} = ['Completed analysis of fov ' obj.fov2str(localFovID) ' at ' datestr(now) ' in ' num2str(toc(fovTimer)) ' s'];
                    display(char(displayStrings)); % Final buffer flush
                end
            end % End loop over fov
        end % End spmd loops
        
    end % End function
    
    % -------------------------------------------------------------------------
    % Combine features to create a final set of features
    % -------------------------------------------------------------------------
    function CombineFeatures(obj)
        % Combine found features determined in individual fov into a
        % single set of non-overlapping, real-world feature boundaries
        % 
        % CombineFeatures();
        
        % Display progress
        if obj.verbose
            PageBreak();
            disp('Combining found features from individual fov');
            totalTimer = tic;
        end
        
        % Create paths to intermediate save states
        allFoundFeaturesPath = [obj.normalizedDataPath obj.segmentationPath 'all_found_features.matb'];
        finalFoundFeaturesPath = [obj.normalizedDataPath obj.segmentationPath 'final_found_features.matb'];

        % Check to see if the combination has already been done and delete
        % if overwrite requested
        if obj.overwrite
            if exist(finalFoundFeaturesPath, 'file')
                if obj.verbose; disp(['...ovewriting final found features list']); end
                delete(finalFoundFeaturesPath);
            end
        end
        
        % Handle existing analysis with overwrite off
        if exist(finalFoundFeaturesPath, 'file')
            disp('...found existing final found features. Skipping analysis.');
            return;
        end
                
        % Compile partial boundaries by fov
        if ~exist(allFoundFeaturesPath, 'file')
            % Confirm that all fov have been properly segmented 
            foundFiles = BuildFileStructure([obj.normalizedDataPath obj.segmentationPath], ...
                'fileExt', 'matb', ...
                'regExp', 'found_features_fov_(?<fov>[0-9]+)', ...
                'fieldNames', {'fov'}, ...
                'fieldConv', {@str2num});

            % Display progress
            if obj.verbose
                disp(['...found ' num2str(length(foundFiles)) ' fov']);
            end

            % Check for missing fov
            missingFovIDs = setdiff(obj.fovIDs, [foundFiles.fov]);
            if ~isempty(missingFovIDs)
                disp(['...the following fov are missing boundaries'])
                disp(['... ' num2str(missingFovIDs)]);
                error('matlabFunctions:missingData', 'Combination cannot proceed until all fov have been segmented');
            end
        
            % Load and concatenate found features 
            if obj.verbose
                disp(['...loading found features']);
                localTimer = tic;
            end
            allFoundFeatures = [];
            for f=1:length(foundFiles)
                allFoundFeatures = cat(1,allFoundFeatures, ...
                    LoadByteStream(foundFiles([foundFiles.fov] == obj.fovIDs(f)).filePath, 'verbose', false));
            end
            
            % Display progress
            if obj.verbose
                disp(['...completed in ' num2str(toc(localTimer)) ' s']);
                saveTimer = tic;
            end
            
            % Save combined 'raw' features
            SaveSplitByteStream(allFoundFeaturesPath, allFoundFeatures, 5000, 'verbose', obj.verbose);
            
            % Delete tempory files for individual fov
            for f=1:length(foundFiles)
                delete(foundFiles(f).filePath);
            end
            
            % Display progress
            if obj.verbose
                disp(['...removed individual files and saved combined file in ' num2str(toc(saveTimer)) ' s']);
            end

        else % Load existing raw feature boundaries
            disp(['...found existing raw boundaries...loading']);
            loadTimer = tic;
            allFoundFeatures = LoadSplitByteStream(allFoundFeaturesPath, 'verbose', obj.verbose);
            disp(['...loaded ' num2str(length(allFoundFeatures)) ' features in ' num2str(toc(loadTimer)) ' s']);
        end
        
        % -------------------------------------------------------------------------
        % Find features that have broken boundaries 
        % -------------------------------------------------------------------------
        % Separate features that are closed from those that are
        % broken
        finalFeatures = allFoundFeatures(~[allFoundFeatures.is_broken]);
        featuresToFix = allFoundFeatures([allFoundFeatures.is_broken]);

        % Display progress
        if obj.verbose
            disp(['...found ' num2str(length(allFoundFeatures)) ' total features']);
            disp(['...found ' num2str(length(finalFeatures)) ' complete features']);
            disp(['...found ' num2str(length(featuresToFix)) ' features with one break']);
            joinTimer = tic;
            fovTimer = tic;
        end

        % Get the primary FOV ids for each cell: these will also serve
        % as a flag to mark features that have been already used to fix
        % other features
        primaryFovIDs = nan(1, length(featuresToFix));
        for F=1:length(featuresToFix)
            primaryFovIDs(F) = featuresToFix(F).ReturnPrimaryFovID();
        end

        for f=1:obj.numFov
            % Find local fov id
            localFovID = obj.fovIDs(f);

            % Find the features associated with this id
            localFeatures = featuresToFix(primaryFovIDs == localFovID);

            % Mark these features as used to prevent future assignment
            primaryFovIDs(primaryFovIDs==localFovID) = nan;

            % Find the surrounding fov ids
            [nnIDX, D] = knnsearch(obj.fovPos, obj.fovPos(f,:), 'K', 5); % 4 neighbors + 1 duplicate of self
            nnIDX = nnIDX(D>0); % Remove self reference

            % Find all possible pairs
            possibleFeaturesForJoin = featuresToFix(ismember(primaryFovIDs, obj.fovIDs(nnIDX)));

            % Loop over all features to fix
            for F=1:length(localFeatures)
                % Compute the join penalty for all possible pairs
                penaltyValues = nan(1, length(possibleFeaturesForJoin)+1);
                for J=1:length(possibleFeaturesForJoin)
                    penaltyValues(J) = localFeatures(F).CalculateJoinPenalty(possibleFeaturesForJoin(J));
                end

                % Compute the self penalty
                penaltyValues(end) = localFeatures(F).CalculateJoinPenalty();

                % Select the minimum
                [minPenalty, minInd] = min(penaltyValues);

                % Determine if the minPenalty is below a given threshold
                if minPenalty < obj.parameters.segmentation.maxEdgeDistance
                    % Determine whether this is a self join or not
                    if minInd <= length(possibleFeaturesForJoin)
                        % Join the feature
                        joinedFeature = localFeatures(F).JoinFeature(possibleFeaturesForJoin(minInd));

                        % Mark the joined feature as used
                        primaryFovIDs(strcmp({featuresToFix.uID}, possibleFeaturesForJoin(minInd).uID)) = nan;

                    else
                        % Join the feature with itself
                        joinedFeature = localFeatures(F).JoinFeature();
                    end

                    % Add the joined feature to the final features list
                    finalFeatures = cat(1, finalFeatures, joinedFeature);
                end
            end
            if obj.verbose
                disp(['...completed fov ' num2str(localFovID) ' in ' num2str(toc(fovTimer)) ' s']);
                fovTimer = tic;
            end
        end % End loop over fov
        
        % Display progress
        if obj.verbose
            % Calculate statistics on joining process
            selfJoined = sum([finalFeatures.num_joined_features]==1);
            pairJoined = sum([finalFeatures.num_joined_features]==2);
            unPaired = length(primaryFovIDs) - selfJoined - pairJoined;

            disp(['...completed all fov in ' num2str(toc(joinTimer)) ' s']);
            disp(['...self-joined ' num2str(selfJoined) ' features']);
            disp(['...pair-joined ' num2str(pairJoined) ' features']);
            disp(['...remaining unpaired features ' num2str(unPaired) ' features']);
            localTimer = tic;
            disp(['...saving final features']);
        end
                      
        % Assign a unique feature id (for fast indexing) to each final
        % feature
        for i=1:length(finalFeatures)
            finalFeatures(i).AssignFeatureID(i);
        end
        
        % Save final joined features
        SaveSplitByteStream(finalFoundFeaturesPath, finalFeatures, 5000, 'verbose', obj.verbose);
        
        % Export a flat file of unique feature ids
        fid = fopen([obj.normalizedDataPath obj.segmentationPath filesep 'featureNames.csv'], 'W');
        featureUIDs = {finalFeatures.uID};
        fprintf(fid, '%s\n', featureUIDs{:});
        fclose(fid);
        
        % Display progress
        if obj.verbose
            disp(['...completed save in ' num2str(toc(localTimer)) ' s']);
            disp(['...completed feature combination in ' num2str(toc(totalTimer)) ' s']);
        end
    end
            
    
    % -------------------------------------------------------------------------
    % Return found features
    % -------------------------------------------------------------------------
    function foundFeatures = GetFoundFeatures(obj)
        % Return found features if they have already been parsed
        % 
        % foundFeatures = GetFoundFeatures();
        
        % Determine if features have been found 
        if ~all(obj.CheckStatus([],'combine'))
            warning('matlabFunctions:incompleteAnalysis', 'Feature segmentation is not complete.');
            foundFeatures = FoundFeature();
        end
        
        % Load and return the found features
        foundFeatures = LoadSplitByteStream([obj.normalizedDataPath obj.segmentationPath 'final_found_features.matb'], ...
            'verbose', obj.verbose);
        
    end
    
    % -------------------------------------------------------------------------
    % Get image frame
    % -------------------------------------------------------------------------
    function [imageFrame, coordinates] = GetImage(obj, fovID, dataChannel, zStack)
        % Return the requested frame
        % 
        % imageFrame = obj.GetImage(fovID, dataChannel, zStack)
        
        % Check input
        if nargin < 3
            error('matlabFunctions:invalidArgument', 'A fov id, a data channel, and z stack must be provided');
        end
        if ~ismember(fovID, obj.fovIDs)
            error('matlabFunctions:invalidArgument', 'The provided fov id is invalid');
        end
        if ~ischar(dataChannel) || ~ismember(dataChannel, {obj.dataOrganization.bitName})
            error('matlabFunctions:invalidArgument', 'The provided data channel name is invalid');
        end
        if zStack < 1 || zStack > obj.numZPos
            error('matlabFunctions:invalidArgument', 'The provided zStack index is invalid');
        end
        
        % Calculate the position of the requested frame
        channelID = find(strcmp({obj.dataOrganization.bitName}, dataChannel));
        frameID = zStack + obj.numZPos*(channelID-1);
        
        % Load the frame
        imageFrame = imread([obj.normalizedDataPath obj.warpedDataPath 'fov_' obj.fov2str(fovID) '.tif'], frameID);
        
        % Define the coordinates
        lowerLeft = obj.Pixel2Abs([1 1], fovID);
        upperRight = obj.Pixel2Abs(obj.imageSize, fovID);
        coordinates.xLimits = [lowerLeft(1) upperRight(1)];
        coordinates.yLimits = [lowerLeft(2) upperRight(2)];

    end
    
    % -------------------------------------------------------------------------
    % Get processed image frame
    % -------------------------------------------------------------------------
    function [imageFrame, coordinates] = GetProcessedImage(obj, fovID, dataChannel, zStack)
        % Return the requested frame
        % 
        % imageFrame = obj.GetProcessedImage(fovID, dataChannel, zStack)
        
        % Check input
        if nargin < 3
            error('matlabFunctions:invalidArgument', 'A fov id, a data channel, and z stack must be provided');
        end
        if ~ismember(fovID, obj.fovIDs)
            error('matlabFunctions:invalidArgument', 'The provided fov id is invalid');
        end
        if ~ischar(dataChannel) || ~ismember(dataChannel, {obj.dataOrganization.bitName})
            error('matlabFunctions:invalidArgument', 'The provided data channel name is invalid');
        end
        if zStack < 1 || zStack > obj.numZPos
            error('matlabFunctions:invalidArgument', 'The provided zStack index is invalid');
        end
        
        % Calculate the position of the requested frame
        channelID = find(strcmp({obj.dataOrganization.bitName}, dataChannel));
        frameID = zStack + obj.numZPos*(channelID-1);
        
        % Load the frame
        imageFrame = imread([obj.normalizedDataPath obj.processedDataPath 'fov_' obj.fov2str(fovID) '.tif'], frameID);
        
        % Define the coordinates
        lowerLeft = obj.Pixel2Abs([1 1], fovID);
        upperRight = obj.Pixel2Abs(obj.imageSize, fovID);
        coordinates.xLimits = [lowerLeft(1) upperRight(1)];
        coordinates.yLimits = [lowerLeft(2) upperRight(2)];

    end

    % -------------------------------------------------------------------------
    % Get Decoded Image
    % -------------------------------------------------------------------------
    function [barcodeStack, magnitudeStack, coordinates] = GetDecodedImage(obj, fovID)
        % Return the decoded image stack
        % 
        % [barcodeStack, magnitudeStack, coordinates] = obj.GetDecodedImage(fovID)
        
        % Check input
        if nargin < 2
            error('matlabFunctions:invalidArgument', 'A fov id must be provided');
        end
        if ~ismember(fovID, obj.fovIDs)
            error('matlabFunctions:invalidArgument', 'The provided fov id is invalid');
        end

        % Check to see if the decoded image is available
        decodedImagePath = [obj.normalizedDataPath obj.barcodePath filesep 'decoded_images' filesep 'fov_' obj.fov2str(fovID) '.tif'];
        if ~exist(decodedImagePath)
            error('matlabFunctions:invalidData', 'The requested image does not exist.');
        end
                
        % Load the barcode stack
        for z=1:obj.numZPos
            % Dynamically find size and allocate memory
            if z==1
                temp = imread(decodedImagePath, z);
                barcodeStack = zeros([size(temp) obj.numZPos]);
                magnitudeStack = barcodeStack;
                barcodeStack(:,:,z) = temp;
                temp = [];
            else
                barcodeStack(:,:,z) = imread(decodedImagePath, z);
            end
        end
        
        % Load the magnitude stack
        for z=1:obj.numZPos
            magnitudeStack(:,:,z) = imread(decodedImagePath, z + obj.numZPos);
        end
        
        % Define the coordinates
        % Handle the crop (somewhat poorly coded)
        xC = 1:obj.imageSize(1);
        yC = 1:obj.imageSize(2);
        
        xC = xC((obj.parameters.decoding.crop+1):(end-obj.parameters.decoding.crop));
        yC = yC((obj.parameters.decoding.crop+1):(end-obj.parameters.decoding.crop));
        
        lowerLeft = obj.Pixel2Abs([xC(1) yC(1)], fovID);
        upperRight = obj.Pixel2Abs([xC(end) yC(end)], fovID);
        coordinates.xLimits = [lowerLeft(1) upperRight(1)];
        coordinates.yLimits = [lowerLeft(2) upperRight(2)];

    end

    
    % -------------------------------------------------------------------------
    % Return barcode lists
    % -------------------------------------------------------------------------
    function bList = GetBarcodeList(obj, fovID)
        % Read barcode lists as requested
        % 
        % bList = GetBarcodeList(fovID);
        
        % Confirm required input
        if ~exist('fovID', 'var') || ~ismember(fovID, obj.fovIDs)
            error('matlabFunctions:invalidArguments', 'A valid fovID must be provided.');
        end
        
        % Check to confirm that the requested fovID has been decoded
        if ~obj.CheckStatus(fovID, 'd')
            error('matlabFunctions:missingData', 'The requested fov has not yet been decoded.');
        end
        
        % Load and return barcode list
        barcodeByFovPath = [obj.normalizedDataPath obj.barcodePath filesep 'barcode_fov' filesep];
        barcodeFile = [barcodeByFovPath 'fov_' obj.fov2str(fovID) '_blist.bin'];

        bList = ReadBinaryFile(barcodeFile);

    end

    % -------------------------------------------------------------------------
    % Return molecule list
    % -------------------------------------------------------------------------
    function mList = GetMoleculeList(obj, fovID)
        % Read a list of molecules
        % 
        % mList = GetMoleculeList(fovID);
        
        % Confirm required input
        if ~exist('fovID', 'var') || ~ismember(fovID, obj.fovIDs)
            error('matlabFunctions:invalidArguments', 'A valid fovID must be provided.');
        end
                
        moleculeByFovPath = [obj.normalizedDataPath filesep 'smFISH' filesep];
        moleculeFile = [moleculeByFovPath 'fov_' obj.fov2str(fovID) '_mlist.bin'];

        if ~exist(moleculeFile)
            error('matlabFunctions:missingFile', 'The requested file does not exist');
        end
        
        mList = ReadBinaryFile(moleculeFile);

    end

    
    % -------------------------------------------------------------------------
    % Return barcode lists
    % -------------------------------------------------------------------------
    function bList = GetParsedBarcodeList(obj, fovID)
        % Read parsed barcode lists as requested
        % 
        % bList = GetParsedBarcodeList(fovID);
        
        % Confirm required input
        if ~exist('fovID', 'var') || ~ismember(fovID, obj.fovIDs)
            error('matlabFunctions:invalidArguments', 'A valid fovID must be provided.');
        end
        
        % Check to confirm that the requested fovID has been decoded
        if ~obj.CheckStatus(fovID, 'p')
            error('matlabFunctions:missingData', 'The requested fov has not yet been parsed.');
        end
        
        % Load and return barcode list
        barcodeByFovPath = [obj.normalizedDataPath obj.barcodePath filesep 'parsed_fov' filesep];
        barcodeFile = [barcodeByFovPath 'fov_' obj.fov2str(fovID) '_blist.bin'];

        bList = ReadBinaryFile(barcodeFile);

    end

    
    
    % -------------------------------------------------------------------------
    % Warp images
    % -------------------------------------------------------------------------
    function WarpFOV(obj, fovIDs)
        % Warp individual fov and produce a warped data stack
        % WarpFOV([]); % Warp all fov
        % WarpFOV(fovIDs); % Warp the fov that match the specified fovids
        
        % -------------------------------------------------------------------------
        % Determine properties of the requested fov ids
        % -------------------------------------------------------------------------
        if isempty(fovIDs)
            fovIDs = obj.fovIDs;
        elseif ~all(ismember(fovIDs, obj.fovIDs))
            error('matlabFunctions:invalidArguments', 'An invalid fov id has been requested');
        end
        
        % -------------------------------------------------------------------------
        % Make directories if they do not exist
        % -------------------------------------------------------------------------
        % Directory for molecule lists
        if ~exist([obj.normalizedDataPath obj.fiducialDataPath], 'dir')
            mkdir([obj.normalizedDataPath obj.fiducialDataPath]);
        end

        % Directory for warped data
        if ~exist([obj.normalizedDataPath obj.warpedDataPath], 'dir')
            mkdir([obj.normalizedDataPath obj.warpedDataPath]);
        end
       
        % -------------------------------------------------------------------------
        % Run processing on individual in parallel (if requested)
        % -------------------------------------------------------------------------
        spmd (obj.numPar) % Run in parallel if requested
            % Loop over requested fov
            for f=labindex:numlabs:length(fovIDs)
				% Determine local fov id
                localFovID = fovIDs(f);
								
                % Create display strings
                if obj.verbose
                    displayStrings = {};
                    displayStrings{end+1} = PageBreak('nodisplay');
                    displayStrings{end+1} = ['Started warping fov ' obj.fov2str(localFovID)];
                end

                % Check to see if the warped tiff file exists and if it is
                % complete
                warpedTiffFileName = ['fov_' obj.fov2str(localFovID) '.tif']; % Define name for warped tiff file
                tiffFileName = [obj.normalizedDataPath obj.warpedDataPath warpedTiffFileName];
				
                if exist(tiffFileName, 'file')
                    if obj.verbose
                        displayStrings{end+1} = ['Found existing warped data file for ' obj.fov2str(localFovID) '. ']; 
                    end
                    
                    if obj.overwrite % Repeat analysis
                        delete(tiffFileName);
                        if obj.verbose
                            displayStrings{end+1} = ['Overwriting'];
                        end
                    else
                        if obj.CheckTiffStack(tiffFileName, obj.numDataChannels*obj.numZPos); % Check for corrupt file
                            delete(tiffFileName);
                            if obj.verbose
                                displayStrings{end+1} = ['File appears to be corrupt. Deleting.'];
                            end
                        else % Skip analysis file exists
                            if obj.verbose
                                displayStrings{end+1} = ['File appears to be complete. Skipping analysis.'];
                                disp(char(displayStrings)); % Display the strings before exit
                            end
                            continue;
                        end
                    end
                end
                
                % Create display strings
                if obj.verbose
                    displayStrings{end+1} = ['Fitting fiducials'];
                    localTimer = tic;
                end

                % Switch based on the method for finding fiducials
                switch obj.parameters.warp.fiducialFitMethod
                    % Handle the daoSTORM case (the only supported option
                    % at this point)
                    case 'daoSTORM'                
                        % Fit Fiducials for each data channel
                        for c=1:obj.numDataChannels
                            % Skip fiducial fitting if no fiducial round is
                            % provided
                            if ~isfield(obj.dataOrganization(c), 'fiducialImageType')
                                continue;
                            end
                            
                            % Find fiducial information
                            if isfield(obj.dataOrganization(c), 'fiducialCameraID')
                                fileInd = find(...
                                    strcmp({obj.rawDataFiles.imageType}, obj.dataOrganization(c).fiducialImageType) & ...
                                    [obj.rawDataFiles.imagingRound] == obj.dataOrganization(c).fiducialImagingRound & ... 
                                    strcmp({obj.rawDataFiles.cameraID}, obj.dataOrganization(c).fiducialCameraID) & ...
                                    [obj.rawDataFiles.fov] == localFovID);

                            else
                                fileInd = find(...
                                    strcmp({obj.rawDataFiles.imageType}, obj.dataOrganization(c).fiducialImageType) & ...
                                    [obj.rawDataFiles.imagingRound] == obj.dataOrganization(c).fiducialImagingRound & ...
                                    [obj.rawDataFiles.fov] == localFovID);
                            end

                            % Check for consistency
                            if length(fileInd) ~= 1
                                disp(fileInd);
                                disp(c);
                                disp(obj.dataOrganization(c).fiducialImageType);
                                disp(obj.dataOrganization(c).fiducialImagingRound);
                                if isfield(obj.dataOrganization(c), 'fiducialCameraID')
                                    disp(obj.dataOrganization(c).fiducialCameraID);
                                end
                                disp(localFovID);
                                error('matlabFunctions:invalidFileInformation', ...
                                    'Either a file is missing or there are multiple files that match an expected pattern.');
                            end

                            % Determine file name
                            localFileName = obj.rawDataFiles(fileInd).name;
                            [~, baseName, ~] = fileparts(localFileName);
                            
                            % Create daostorm parameters file
                            WriteDaoSTORMParameters([obj.normalizedDataPath obj.fiducialDataPath baseName '_dao.xml'], ...  
                                'start_frame', obj.dataOrganization(c).fiducialFrame-1, ...
                                'max_frame', obj.dataOrganization(c).fiducialFrame, ...
                                'x_stop', obj.imageSize(1), 'y_stop', obj.imageSize(2), ...
                                'iterations', 1, ...
                                'pixel_size', obj.pixelSize, ...
                                'threshold', obj.parameters.warp.daoThreshold, ...
                                'sigma', obj.parameters.warp.sigmaInit, ...
                                'baseline', obj.parameters.warp.daoBaseline);

                            % Run daoSTORM
                            daoSTORM([obj.rawDataPath localFileName], ...
                                [obj.normalizedDataPath obj.fiducialDataPath baseName '_dao.xml'], ...
                                'overwrite', true, ...                      % Overwrite all files (overwrite protection is provided above)
                                'numParallel', 1, ...
                                'savePath', [obj.normalizedDataPath obj.fiducialDataPath], ...
                                'verbose', false, ...
                                'waitTime', 1, ...
                                'outputInMatlab', true, ...
                                'hideterminal', true);

                        end % End loop over number of data channels

                        % Create display strings
                        if obj.verbose
                            displayStrings{end+1} = ['... completed in ' num2str(toc(localTimer)) ' s'];
                            displayStrings{end+1} = PageBreak('nodisplay');
                            localTimer = tic; % Restart timer
                            displayStrings{end+1} = ['Constructing affine transforms'];
                        end

                        % Build affine transformation
                        for c=1:obj.numDataChannels
                            % Skip warping if no fiducial information is
                            % provided
                            if ~isfield(obj.dataOrganization(c), 'fiducialImageType')
                                localAffine(c) = affine2d();
                                localResiduals{c} = zeros(0,4);
                                continue;
                            end
                            
                            % Identify name of mList
                            if isfield(obj.dataOrganization(c), 'fiducialCameraID')
                                fileInd = ...
                                    strcmp({obj.rawDataFiles.imageType}, obj.dataOrganization(c).fiducialImageType) & ...
                                    [obj.rawDataFiles.imagingRound] == obj.dataOrganization(c).fiducialImagingRound & ...
                                    strcmp({obj.rawDataFiles.cameraID}, obj.dataOrganization(c).fiducialCameraID) & ...
                                    [obj.rawDataFiles.fov] == localFovID;
                            else
                                fileInd = ...
                                    strcmp({obj.rawDataFiles.imageType}, obj.dataOrganization(c).fiducialImageType) & ...
                                    [obj.rawDataFiles.imagingRound] == obj.dataOrganization(c).fiducialImagingRound & ...
                                    [obj.rawDataFiles.fov] == localFovID;
                            end

                            % Determine file name
                            localFileName = obj.rawDataFiles(fileInd).name;
                            [~, baseName, ~] = fileparts(localFileName);

                            % Fiducial file name
                            fiducialFileName = [obj.normalizedDataPath obj.fiducialDataPath baseName '_mList.bin'];

                            % Check to see if it exists
                            if ~exist(fiducialFileName, 'file')
                                disp(['Could not find ' fiducialFileName]);
                                error('matlabFunctions:missingFile', ...
                                    'Could not find a fiducial file.');
                            end

                            % Load molecule list
                            movList = ReadMasterMoleculeList(fiducialFileName, ...
                                'fieldsToLoad', {'x','y','xc','yc','frame'}, ...
                                'verbose', false);

                            % Define reference list as the list from the
                            % first data channel (i.e. bit 1)
                            if c==1
                                refList = movList;
                            end

                            % Generate transform
                            [localAffine(c), ~, localResiduals{c}] = MLists2Transform(refList, movList, ...
                                'ignoreFrames', true, ...
                                'controlPointMethod', 'kNNDistanceHistogram', ...
                                'histogramEdges', obj.parameters.warp.controlPointOffsetRange, ...
                                'numNN', obj.parameters.warp.numNN, ...
                                'pairDistTolerance', obj.parameters.warp.pairDistanceTolerance);
                            
                            % Generate additional input
                            if obj.verbose
                                displayStrings{end+1} = ['...completed affine transform for data channel ' num2str(c)];
                                displayStrings{end+1} = ['...out of ' num2str(length(refList.x)) ' reference molecules and ' num2str(length(movList.x)) ' moving molecules matched ' num2str(size(localResiduals{c},1))];
                            end
                        end % End construction of affine transformations

                        % Create display strings
                        if obj.verbose
                            displayStrings{end+1} = ['... completed in ' num2str(toc(localTimer)) ' s'];
                            disp(char(displayStrings)); % Display the strings
                            displayStrings = {};

                        end
                        
                        % Save data associated with transform for this fov
                        SaveAsByteStream([obj.normalizedDataPath obj.fiducialDataPath 'fov_' obj.fov2str(localFovID) '_affine.matb'], ...
                            localAffine, 'verbose', obj.verbose);
                        SaveAsByteStream([obj.normalizedDataPath obj.fiducialDataPath 'fov_' obj.fov2str(localFovID) '_residuals.matb'], ...
                            localResiduals, 'verbose', obj.verbose);

                    otherwise % Handle the case that an undefined fiducial/warp approach was defined
                        error('matlabFunctions:invalidArguments', ...
                            'The provided fiducial fit method is not supported');
                    
                end % End switch statement on fiducial fitting/affine transformation construction approach
               
                % Create display strings
                if obj.verbose
                    displayStrings{end+1} = PageBreak('nodisplay');
                    localTimer = tic; % Restart timer
                    displayStrings{end+1} = ['Writing warped tiff file'];
                end
                
                
                % Export a warped tiff stack of fiducials only: Run this
                % first as a complete warped tiff stack is the measure of a
                % complete warp process. 
                if obj.parameters.warp.exportWarpedBeads
                    
                    if obj.verbose
                        displayStrings{end+1} = ['Creating warped fidicual stack...'];
                    end
                    
                    % Define the folder for defaults if it does not exist
                    fiducialsPath = [obj.normalizedDataPath obj.warpedDataPath filesep 'warped_fiducials' filesep];
                    if ~exist(fiducialsPath, 'dir')
                        mkdir(fiducialsPath);
                    end
                    
                    % Define the bead tiff file name
                    fiducialsTiffFileName = [fiducialsPath 'fov_' obj.fov2str(localFovID) '.tif'];
                    
                    % Overwrite any previous analysis
                    if exist(fiducialsTiffFileName, 'dir')
                        if obj.verbose
                            displayStrings{end+1} = ['...found existing warped fidicual stack.... removing...'];
                        end
                        delete(fiducialsTiffFileName);
                    end

                    
                    % Create tiff file
                    fiducialsTiffFile = Tiff(fiducialsTiffFileName, 'w8');

                    % Create tiff tags
                    fiducialsTiffTagStruct.ImageLength = obj.imageSize(1);
                    fiducialsTiffTagStruct.ImageWidth = obj.imageSize(2);
                    fiducialsTiffTagStruct.Photometric = Tiff.Photometric.MinIsBlack;
                    fiducialsTiffTagStruct.BitsPerSample = 16;
                    fiducialsTiffTagStruct.SamplesPerPixel = 1;
                    fiducialsTiffTagStruct.RowsPerStrip = 16;
                    fiducialsTiffTagStruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
                    fiducialsTiffTagStruct.Software = 'MATLAB';
                    fiducialsTiffTagStruct.ImageDescription = sprintf(['ImageJ=1.47a\n' ...      % ImageJ label
                        'images=' num2str(obj.numDataChannels) '\n' ...     % The total number of images
                        'channels=1\n' ...                                              % The number of channels, 1 for now
                        'slices=1\n' ...                         % The number of z slices
                        'frames=' num2str(obj.numDataChannels) '\n' ...                     % The number of frames, i.e. data org entries
                        'hyperstack=true\n' ...
                        'loop=false\n' ...
                        ]);                
                    
                    % Loop over all data channels
                    try
                        for c=1:obj.numDataChannels
                            % Skip fiducial fitting if no fiducial round is
                            % provided
                            if ~isfield(obj.dataOrganization(c), 'fiducialImageType')
                                displayStrings{end+1} = ['...no fiducials for channel ' num2str(c)];
                                continue;
                            end
                            
                            % Find fiducial information
                            if isfield(obj.dataOrganization(c), 'fiducialCameraID')
                                fileInd = find(...
                                    strcmp({obj.rawDataFiles.imageType}, obj.dataOrganization(c).fiducialImageType) & ...
                                    [obj.rawDataFiles.imagingRound] == obj.dataOrganization(c).fiducialImagingRound & ... 
                                    strcmp({obj.rawDataFiles.cameraID}, obj.dataOrganization(c).fiducialCameraID) & ...
                                    [obj.rawDataFiles.fov] == localFovID);

                            else
                                fileInd = find(...
                                    strcmp({obj.rawDataFiles.imageType}, obj.dataOrganization(c).fiducialImageType) & ...
                                    [obj.rawDataFiles.imagingRound] == obj.dataOrganization(c).fiducialImagingRound & ...
                                    [obj.rawDataFiles.fov] == localFovID);
                            end

                            % Determine file name
                            localFileName = obj.rawDataFiles(fileInd).name;
                            
                            % Switch based on the imageExt to load fiducial
                            % image
                            localImage = [];
                            switch obj.imageExt
                                case 'dax'
                                    % Read single data frame
                                    localImage = ReadDax([obj.rawDataPath localFileName], ...
                                        'startFrame', obj.dataOrganization(c).fiducialFrame, ...
                                        'endFrame', obj.dataOrganization(c).fiducialFrame, ...
                                        'verbose', false);
                                case {'tiff','tif'}
                                    % Read tiff frame
                                    localImage = imread([obj.rawDataPath localFileName], ...
                                        obj.dataOrganization(c).fiducialFrame);
                                otherwise
                                    error('matlabFunctions:unsupportedFile', 'The specified image file ext is not yet supported.');
                            end
                            
                            % Check for image load problems
                            if isempty(localImage)
                                error('matlabFunctions:invalidData', 'The requested frame does not exist');
                            end

                            % Warp the image
                            ra = imref2d(size(localImage)); % Create a crop
                            localImage = imwarp(localImage, localAffine(c), 'OutputView', ra);

                            % Write tiff file
                            fiducialsTiffFile.setTag(fiducialsTiffTagStruct);
                            fiducialsTiffFile.write(localImage);
                            if c~= obj.numDataChannels
                                fiducialsTiffFile.writeDirectory(); % Write the directory for the next frame
                            end
                            
                            if obj.verbose
                                displayStrings{end+1} = ['...completed channel ' num2str(c) ' of ' num2str(obj.numDataChannels)];
                            end
                        end % end loop over channels
                    catch ME
                        % Close tiff file and delete existing file
                        fiducialsTiffFile.close();
                        delete(fiducialsTiffFileName);

                        % Update log
                        displayStrings{end+1} = ['Encountered error: ' ME.identifier];
                        disp(char(displayStrings));

                        % Rethrow err
                        rethrow(ME);
                    end
                    if obj.verbose
                        displayStrings{end+1} = ['...completed warped fiducial stack'];
                    end
                end % End if statement for creating a warped fiducial stack
                
                % Create tiff file
                tiffFile = Tiff(tiffFileName, 'w8');
                
                % Create tiff tags
                tiffTagStruct.ImageLength = obj.imageSize(1);
                tiffTagStruct.ImageWidth = obj.imageSize(2);
                tiffTagStruct.Photometric = Tiff.Photometric.MinIsBlack;
                tiffTagStruct.BitsPerSample = 16;
                tiffTagStruct.SamplesPerPixel = 1;
                tiffTagStruct.RowsPerStrip = 16;
                tiffTagStruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
                tiffTagStruct.Software = 'MATLAB';
                tiffTagStruct.ImageDescription = sprintf(['ImageJ=1.47a\n' ...      % ImageJ label
                    'images=' num2str(obj.numDataChannels*obj.numZPos) '\n' ...     % The total number of images
                    'channels=1\n' ...                                              % The number of channels, 1 for now
                    'slices=' num2str(obj.numZPos) '\n' ...                         % The number of z slices
                    'frames=' num2str(obj.numDataChannels) '\n' ...                     % The number of frames, i.e. data org entries
                    'hyperstack=true\n' ...
                    'loop=false\n' ...
                    ]);
                
                if obj.verbose
                    displayStrings{end+1} = ['Creating the warped tiff stack...'];
                end
                
                % Gracefully handle external kill commands
                try
                    % Load, warp, and write tiff files
                    for c=1:obj.numDataChannels
                        % Identify file properties based on file organization
                        imageType = obj.dataOrganization(c).imageType;
                        imagingRound = obj.dataOrganization(c).imagingRound;
                        color = obj.dataOrganization(c).color;
                        frames = obj.dataOrganization(c).frame;
                        localZPos = obj.dataOrganization(c).zPos;
                        if ~isfield(obj.dataOrganization(c), 'imagingCameraID')
                            localCameraID = '';
                        else
                            localCameraID = obj.dataOrganization(c).imagingCameraID;
                        end

                        % Identify dax file
                        localFile = obj.rawDataFiles([strcmp({obj.rawDataFiles.imageType}, imageType) & ...
                            [obj.rawDataFiles.fov] == localFovID & ...
                            strcmp({obj.rawDataFiles.cameraID}, localCameraID) & ...
                            [obj.rawDataFiles.imagingRound] == imagingRound]);

                        % Identify the affine transform to be used to adjust this color channel
                        localIDs = find(strcmp({obj.parameters.warp.colorTransforms.color}, ....
                            num2str(obj.dataOrganization(c).color))); % Find all transforms that match the specified color (they will be applied in order)
                        
                        % Extract this transform
                        if ~isempty(localIDs)
                            localTransforms = obj.parameters.warp.colorTransforms(localIDs);
                        else
                            localTransforms = [];
                        end
                        
                        % Loop over z and load dax
                        for z=1:length(obj.zPos)
                            % Find frame
                            frameInd = find(localZPos == obj.zPos(z));

                            % If the position doesn't exist, assume it is 1 (this
                            % will replicate frames for data sets that do not have
                            % the all z stacks
                            if isempty(frameInd)
                                frameInd = 1;
                            end

                            % Switch based on the imageExt
                            switch obj.imageExt
                                case 'dax'
                                    % Read single data frame
                                    localImage = ReadDax(localFile.filePath, 'startFrame', frames(frameInd), ...
                                        'endFrame', frames(frameInd), 'verbose', false);
                                case {'tiff','tif'}
                                    % Read tiff frame
                                    localImage = imread(localFile.filePath, frames(frameInd));
                                otherwise
                                    error('matlabFunctions:unsupportedFile', 'The specified image file ext is not yet supported.');
                            end
                            
                            % Check for image load problems
                            if isempty(localImage)
                                error('matlabFunctions:invalidData', 'The requested frame does not exist');
                            end
                            
                            % Apply a color based transform if needed
                            ra = imref2d(size(localImage)); % Create a crop
                            if ~isempty(localTransforms)
                                % Loop over chromatic transforms and apply
                                % in order
                                for T=1:length(localTransforms)
                                    switch localTransforms(T).type
                                        case {'similarity', 'Similarity', 'Euclidean', 'euclidean'}
                                            % Create the affine transform based
                                            % on the provided data
                                            trans = affine2d(localTransforms(T).transform);
                                            % Apply the transform
                                            localImage = imwarp(localImage, trans, 'OutputView', ra);
                                            
                                        case {'Invert', 'invert'}
                                            % Invert the axes of the image
                                            % as needed
                                            % Flip x axis
                                            if localTransforms(T).transform(1)
                                                localImage = localImage(:, end:-1:1);
                                            end
                                            % Flip y axis
                                            if localTransforms(T).transform(2)
                                                localImage = localImage(end:-1:1, :);
                                            end
                                            
                                        otherwise
                                            error('matlabFunctions:unsupportedFile', 'The transform type provided is not supported');
                                    end
                                end
                            end
                            
                            % Remove translation of the image due to stage
                            % alignment using the fiducial bead affine
                            % transform
                            localImage = imwarp(localImage, localAffine(c), 'OutputView', ra);
                            
                            % Reorient image to a fixed orientation
                            % Image width represents X,
                            % Image height represents Y, 
                            % X/Y increase with increasing pixel id
                            cameraOrientation = obj.parameters.warp.cameraOrientation;
                            if cameraOrientation(3)
                                localImage = transpose(localImage); % Exchange X/Y
                            end
                            if cameraOrientation(1)
                                localImage = flip(localImage, 2); % Invert the X axis (rows)
                            end
                            if cameraOrientation(2)
                                localImage = flip(localImage, 1); % Invert the Y axis (columns)
                            end
                            
                            % Write tiff file
                            tiffFile.setTag(tiffTagStruct);
                            tiffFile.write(localImage);
                            if z~=length(obj.zPos) || c~= obj.numDataChannels
                                tiffFile.writeDirectory(); % Write the directory for the next frame
                            end
                            
                            if obj.verbose
                                displayStrings{end+1} = ['...completed channel ' num2str(c) ' of ' num2str(obj.numDataChannels)];
                            end

                        end % End loop over z

                    end % End loop over channels
                catch ME
                    % Close tiff file and delete existing file
                    tiffFile.close();
                    delete(tiffFileName);
                    
                    % Update log
                    displayStrings{end+1} = ['Encountered error: ' ME.identifier];
                    display(char(displayStrings));
                    
                    % Rethrow err
                    rethrow(ME);
                end
                
                % Close tiff file
                tiffFile.close();
                
                % Create display strings
                if obj.verbose
                    displayStrings{end+1} = ['... completed in ' num2str(toc(localTimer)) ' s'];
                    displayStrings{end+1} = ['Completed warping of fov ' obj.fov2str(localFovID) ' at ' datestr(now)];
                    % Display all display commands
                    display(char(displayStrings));
                end

            end % End loop over fovIds
        end % end spmd loop
        
    end    
    
    % -------------------------------------------------------------------------
    % Preprocess Images
    % -------------------------------------------------------------------------
    function PreprocessFOV(obj, fovIDs)
        % Preprocess individual FOV and produce a processed tiff stack
        % PreprocessFOV([]); % Deconvolve all fov
        % PreprocessFOV(fovIDs); % Deconvolve the fov that match the specified fovids
        
        % -------------------------------------------------------------------------
        % Determine properties of the requested fov ids
        % -------------------------------------------------------------------------
        if isempty(fovIDs)
            fovIDs = obj.fovIDs;
        elseif ~all(ismember(fovIDs, obj.fovIDs))
            error('matlabFunctions:invalidArguments', 'An invalid fov id has been requested');
        end
        
        % -------------------------------------------------------------------------
        % Check to see if the pixel histogram field has already been
        % populated--no need to repeat analysis if it does
        % -------------------------------------------------------------------------
        if isempty(obj.pixelHistograms) || obj.overwrite
            % -------------------------------------------------------------------------
            % Make directories if they do not exist
            % -------------------------------------------------------------------------
            % Directory for molecule lists
            if ~exist([obj.normalizedDataPath obj.processedDataPath], 'dir')
                mkdir([obj.normalizedDataPath obj.processedDataPath]);
            end

            % Create Generic Tiff Tag
            tiffTagStruct.Photometric = Tiff.Photometric.MinIsBlack;
            tiffTagStruct.BitsPerSample = 16;
            tiffTagStruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
            tiffTagStruct.Software = 'MATLAB';
            tiffTagStruct.ImageLength = obj.imageSize(1);
            tiffTagStruct.ImageWidth = obj.imageSize(2);

            tiffTagStruct.ImageDescription = sprintf(['ImageJ=1.47a\n' ... % ImageJ label
                'images=' num2str(obj.numBits*obj.numZPos) '\n' ... % The total number of images
                'channels=1\n' ...             % The number of channels, 1 for now
                'slices=' num2str(obj.numZPos) '\n' ... % The number of z slices
                'frames=' num2str(obj.numBits) '\n' ... % The number of frames, i.e. data org entries
                'hyperstack=true\n' ...
                'loop=false\n' ...
                ]); % Convert \n to newline characters

            % Loop over fov in parallel
            spmd (obj.numPar)
                for f=labindex:numlabs:length(fovIDs)
                    % Determine local fovID
                    localFovID = fovIDs(f);

                    % Create display strings
                    if obj.verbose
                        displayStrings = {};
                        displayStrings{end+1} = PageBreak('nodisplay');
                        displayStrings{end+1} = ['Started preprocessing of fov ' obj.fov2str(localFovID)];
                        localTimer = tic;
                    end

                    % Create tiff to read file name
                    tiffName2Read = [obj.normalizedDataPath obj.warpedDataPath ...
                        'fov_' obj.fov2str(localFovID) '.tif'];
                    if ~exist(tiffName2Read, 'file')
                        error('matlabFunctions:missingFile', 'The requsted tiff stack is not present.');
                    end

                    % Create tiff to write file name
                    tiffName2Write = [obj.normalizedDataPath obj.processedDataPath ...
                        'fov_' obj.fov2str(localFovID) '.tif'];

                    % Create pixel histogram file and initialize
                    pixelHistogramFile = [obj.normalizedDataPath obj.processedDataPath 'pixel_histograms' filesep ...
                        'pixel_data_fov_' obj.fov2str(localFovID) '.matb'];

                    if ~exist([obj.normalizedDataPath obj.processedDataPath 'pixel_histograms' filesep], 'dir')
                        mkdir([obj.normalizedDataPath obj.processedDataPath 'pixel_histograms' filesep]);
                    end

                    pixelHistogram = zeros(obj.numBits, uint16(inf)); % Allocate memory and initialize to zero

                    if exist(tiffName2Write, 'file')
                        if obj.verbose
                            displayStrings{end+1} = ['Found existing warped data file for ' obj.fov2str(localFovID) '. ']; 
                        end

                        if obj.overwrite % Repeat analysis
                            if obj.verbose; displayStrings{end+1} = ['Overwriting']; end;
                            delete(tiffName2Write);
                        else
                            % Check for corrupt tiff file
                            if obj.CheckTiffStack(tiffName2Write, obj.numBits*obj.numZPos); 
                                if obj.verbose; displayStrings{end+1} = ['File appears to be corrupt. Deleting.']; end;
                                delete(tiffName2Write);
                            elseif ~exist(pixelHistogramFile, 'file') % Check for missing pixel histogram
                                if obj.verbose; displayStrings{end+1} = ['Pixel histogram is missing. Repeating analysis']; end;
                                delete(tiffName2Write);                            
                            else % Handle the case that the file is not corrupt and the pixel histogram file exists, i.e. the analysis is complete
                                if obj.verbose
                                    displayStrings{end+1} = ['File appears to be complete. Skipping analysis.'];
                                    display(char(displayStrings)); % Display the strings before exit
                                end
                                continue;
                            end
                        end
                    end

                    % Open tiff files
                    tiffToRead = Tiff(tiffName2Read, 'r');

                    % Create writing tiff
                    tiffToWrite = Tiff(tiffName2Write, 'w8');

                    % Gracefully handle interrupts
                    try
                        % Switch to run the specific algorithm
                        switch obj.parameters.preprocess.preprocessingMethod

                            case 'highPassDecon'
                                % Loop over imaging round
                                for b=1:obj.numBits
                                    % Loop over z
                                    for z=1:obj.numZPos
                                        % Set directory and load
                                        tiffToRead.setDirectory((b-1)*obj.numZPos + z);
                                        localFrame = uint16(tiffToRead.read());

                                        % High pass filter (and threshold on zero values b/c uint16)
                                        localFrame1 = uint16(localFrame) - ...
                                            uint16(imgaussfilt(localFrame, obj.parameters.preprocess.highPassKernelSize));

                                        % Deconvolve
                                        if obj.parameters.preprocess.numIterDecon > 0
                                            localFrame2 = deconvlucy(localFrame1, ...
                                                obj.parameters.preprocess.deconKernel, ...
                                                obj.parameters.preprocess.numIterDecon);
                                        else
                                            localFrame2 = localFrame1;
                                        end
                                        
                                        % Write frame
                                        tiffToWrite.setTag(tiffTagStruct);
                                        tiffToWrite.write(localFrame2);

                                        % Write directory for next frame (unless it is the last
                                        % one)
                                        if ((b-1)*obj.numZPos + z) ~= obj.numBits*obj.numZPos
                                            tiffToWrite.writeDirectory();
                                        end

                                        % Accumulate pixel histograms
                                        pixelHistogram(b,:) = pixelHistogram(b,:) + hist(localFrame2(:), 0:double(uint16(inf)-1));
                                    end 
                                end

                            case 'highPassErosion'
                                % Loop over imaging round
                                for b=1:obj.numBits
                                    % Loop over z
                                    for z=1:obj.numZPos
                                        % Set directory and load
                                        tiffToRead.setDirectory((b-1)*obj.numZPos + z);
                                        localFrame = uint16(tiffToRead.read());

                                        % High pass filter (and threshold on zero values b/c uint16)
                                        localFrame1 = uint16(localFrame) - ...
                                            uint16(imgaussfilt(localFrame, obj.parameters.preprocess.highPassKernelSize));

                                        % Erode
                                        localFrame2 = imerode(localFrame1, obj.parameters.preprocess.erosionElement);

                                        % Write frame
                                        tiffToWrite.setTag(tiffTagStruct);
                                        tiffToWrite.write(localFrame2);
                                        
                                        % Write directory for next frame (unless it is the last
                                        % one)
                                        if ((b-1)*obj.numZPos + z) ~= obj.numBits*obj.numZPos
                                            tiffToWrite.writeDirectory();
                                        end

                                        % Accumulate pixel histograms
                                        pixelHistogram(b,:) = pixelHistogram(b,:) + hist(localFrame2(:), 0:double(uint16(inf)-1));
                                    end 
                                end

                            case 'highPassDeconWB'
                                % Loop over imaging round
                                for b=1:obj.numBits
                                    % Loop over z
                                    for z=1:obj.numZPos
                                        % Set directory and load
                                        tiffToRead.setDirectory((b-1)*obj.numZPos + z);
                                        localFrame = uint16(tiffToRead.read());

                                        % High pass filter (and threshold on zero values b/c uint16)
                                        localFrame1 = uint16(localFrame) - ...
                                            uint16(imgaussfilt(localFrame, obj.parameters.preprocess.highPassKernelSize));

                                        % Deconvolve
                                        localFrame2 = WBDeconvolution(localFrame1, ...
                                            obj.parameters.preprocess.deconKernel, ...
                                            obj.parameters.preprocess.numIterDecon);

                                        % Write frame
                                        tiffToWrite.setTag(tiffTagStruct);
                                        tiffToWrite.write(localFrame2);
                                        
                                        % Write directory for next frame (unless it is the last
                                        % one)
                                        if ((b-1)*obj.numZPos + z) ~= obj.numBits*obj.numZPos
                                            tiffToWrite.writeDirectory();
                                        end

                                        % Accumulate pixel histograms
                                        pixelHistogram(b,:) = pixelHistogram(b,:) + hist(localFrame2(:), 0:double(uint16(inf)-1));
                                    end 
                                end

                            otherwise
                                error('matlabFunctions:invalidArguments', 'Invalid preprocessing method');
                        end

                    catch ME
                        % Close tiff file and delete partially constructed file
                        tiffToWrite.close();
                        delete(tiffName2Write);

                        % Update log
                        displayStrings{end+1} = ['Encountered error: ' ME.identifier];
                        display(char(displayStrings));

                        % Rethrow error
                        rethrow(ME);
                    end

                    % Close tiff files
                    tiffToWrite.close();
                    tiffToRead.close();

                    % Write the pixel histogram file
                    SaveAsByteStream(pixelHistogramFile, pixelHistogram, 'verbose', false);

                    % Create display strings
                    if obj.verbose
                        displayStrings{end+1} = ['... completed in ' num2str(toc(localTimer)) ' s'];
                        displayStrings{end+1} = ['Completed preprocessing of fov ' obj.fov2str(localFovID) ' at ' datestr(now)];
                        % Display all display commands
                        display(char(displayStrings));
                    end

                end
            end % spmd loop
        else
            warning('The pixel histograms field is not empty, indicating that preprocessing is complete on this data set.');
        end
    end
    
    % -------------------------------------------------------------------------
    % Check that various analysis steps have been completed on a fov 
    % -------------------------------------------------------------------------
    function complete = CheckStatus(obj, fovIDs, analysisType)
        % Check to see if a given analysis step has been completed on a set
        % of fovs by ID
        % complete = obj.CheckStatus(fovIDs, ...); % Check given fov ids
        % complete = obj.CheckStatus([], ...); % Check all fov ids
        % complete = obj.CheckStatus(..., 'warp'); % Is warping complete?
        % complete = obj.CheckStatus(..., 'process'); % Is preprocessing
        %    complete?
        % complete = obj.CheckStatus(..., 'optimize'); % Is optimization
        %    complete?
        % complete = obj.CheckStatus(..., 'decode'); % Is decoding
        %    complete?
        % complete = obj.CheckStatus(..., 'parse'); % Is parsing
        %    complete?
        % complete = obj.CheckStatus(..., 'segment'); % Is segmentation
        %    complete?
        % complete = obj.CheckStatus(..., 'combine'); % Is boundary
        %    combination complete?
        % complete = obj.CheckStatus(..., 'mosaic'); % Is the low
        %    resolution mosaic construction complete?
        
        % -------------------------------------------------------------------------
        % Determine properties of the requested fov ids
        % -------------------------------------------------------------------------
        if isempty(fovIDs)
            fovIDs = obj.fovIDs;
        elseif ~all(ismember(fovIDs, obj.fovIDs))
            error('matlabFunctions:invalidArguments', 'An invalid fov id has been requested');
        end
        
        % -------------------------------------------------------------------------
        % Loop over requested fov ids
        % -------------------------------------------------------------------------
        complete = false(1, length(fovIDs));
        for f=1:length(fovIDs)
            % Make local copy of fovID
            localFovID = fovIDs(f);
            
            % -------------------------------------------------------------------------
            % Switch based on analysis type
            % -------------------------------------------------------------------------
            switch analysisType
                case {'warp'}
                    % Create tiff to write file name
                    tiffName2Write = [obj.normalizedDataPath obj.warpedDataPath ...
                        'fov_' obj.fov2str(localFovID) '.tif'];
                                        
                    % Check for corrupt tiff file
                    complete(f) = ~obj.CheckTiffStack(tiffName2Write, obj.numDataChannels*obj.numZPos); % The pixel histogram has been created

                case {'process', 'preprocess', 'w'} % Handle preprocess                    
                    % Create tiff to write file name
                    tiffName2Write = [obj.normalizedDataPath obj.processedDataPath ...
                        'fov_' obj.fov2str(localFovID) '.tif'];
                    
                    % Check for compiled pixel histograms
                    histogramsCompiled = ~isempty(obj.pixelHistograms);
                    
                    % Name of pixel histogram file
                    pixelHistogramFile = [obj.normalizedDataPath obj.processedDataPath 'pixel_histograms' filesep ...
                        'pixel_data_fov_' obj.fov2str(localFovID) '.matb'];
                    
                    % Check for corrupt tiff file
                    complete(f) = ~obj.CheckTiffStack(tiffName2Write, obj.numBits*obj.numZPos) && ...% File exists and is not corrupt
                        (histogramsCompiled || exist(pixelHistogramFile, 'file')); % The pixel histogram has been created
                        
                % Check the decode state
                case {'decode', 'decoding', 'd'}
                    % Define barcode file path and check for uncorrupted existance
                    barcodeByFovPath = [obj.normalizedDataPath obj.barcodePath filesep 'barcode_fov' filesep];
                    barcodeFile = [barcodeByFovPath 'fov_' obj.fov2str(localFovID) '_blist.bin'];
                    
                    % Check if exists and is not corrupt
                    if exist(barcodeFile, 'file')
                        isCorrupt = false;
                        try     
                            fileHeader = ReadBinaryFileHeader(barcodeFile);
                            isCorrupt = fileHeader.isCorrupt;
                        catch
                            isCorrupt = true;
                        end
                        complete(f) = ~isCorrupt;
                    else
                        complete(f) = false;
                    end
                    
                % Check optimization
                case {'optimize', 'o', 'opt'}
                    % Update the scale factors field (another instance of
                    % this decoder could have written this value)
                    obj.LoadField('scaleFactors');

                    % Determine if the scale factors have been set
                    complete(f) = ~isempty(obj.scaleFactors);
                    
                % Check segmentation
                case {'segment', 'segmentation', 's'}
                    % Create file path and check for existance
                    foundFeaturesPath = [obj.normalizedDataPath obj.segmentationPath ...
                        'found_features_fov_'  obj.fov2str(localFovID) '.matb'];
                    complete(f) = exist(foundFeaturesPath, 'file') || ...   % Either the file
                        exist([obj.normalizedDataPath obj.segmentationPath 'final_found_features.matb'], 'file'); % or the combined information

                % Check combining boundaries
                case {'combine', 'c', 'comb'}
                    complete(f) = exist([obj.normalizedDataPath obj.segmentationPath 'final_found_features.matb'], ...
                        'file');
                    
                % Check the parse state
                case {'parse', 'parsing', 'p'}
                    % Define barcode file path and check for uncorrupted existance
                    parsedBarcodePath = [obj.normalizedDataPath obj.barcodePath filesep 'parsed_fov' filesep];
                    barcodeFile = [parsedBarcodePath 'fov_' obj.fov2str(localFovID) '_blist.bin'];
                    
                    % Check if exists and is not corrupt
                    if exist(barcodeFile, 'file')
                        isCorrupt = false;
                        try     
                            fileHeader = ReadBinaryFileHeader(barcodeFile);
                            isCorrupt = fileHeader.isCorrupt;
                        catch
                            isCorrupt = true;
                        end
                        complete(f) = ~isCorrupt;
                    else
                        complete(f) = false;
                    end
                    
                % Check the raw summation state
                case {'sum', 'r'}
                    totalSignalFilePath = [obj.normalizedDataPath obj.summationPath ...
                            'total_signal_fov_'  obj.fov2str(localFovID) '.matb'];
                
                    numberOfPixelsFilePath = [obj.normalizedDataPath obj.summationPath ...
                            'total_pixels_fov_'  obj.fov2str(localFovID) '.matb'];
                  
                    combinedSignalFilePath = [obj.normalizedDataPath obj.summationPath ...
                            'total_signal.csv'];
                        
                    combinedPixelsFilePath = [obj.normalizedDataPath obj.summationPath ...
                            'total_pixels.csv'];
                        
                    if exist(totalSignalFilePath, 'file') && exist(numberOfPixelsFilePath, 'file') 
                        complete(f) = true;
                    elseif exist(combinedSignalFilePath, 'file') && exist(combinedPixelsFilePath, 'file')
                        complete(f) = true;
                    else
                        complete(f) = false;
                    end
                    
                case {'combineSum', 'combine_sum', 'i'}
                    combinedSignalFilePath = [obj.normalizedDataPath obj.summationPath ...
                            'total_signal.csv'];
                        
                    combinedPixelsFilePath = [obj.normalizedDataPath obj.summationPath ...
                            'total_pixels.csv'];
                        
                    if exist(combinedSignalFilePath, 'file') && exist(combinedPixelsFilePath, 'file')
                        complete(f) = true;
                    else
                        complete(f) = false;
                    end
                    
                case {'numbers', 'n'}
                    % Determine the expected file names
                    reportsPath = [obj.normalizedDataPath 'reports' filesep];
                    fileNames = {'countsPerCellExactIn.csv', ...
                        'countsPerCellCorrectedIn.csv', ...
                        'countsPerCellExactOut.csv', ...
                        'countsPerCellCorrectedOut.csv'};
                    
                    % Determine which files exist
                    doesExist = false(1, length(fileNames));
                    for d=1:length(fileNames)
                        doesExist(d) = exist([reportsPath fileNames{d}]);
                    end
                    
                    % They must all exist to be marked as complete
                    complete(f) = all(doesExist);
                    
                case {'mosaic', 'l'}
                    % Determine the number of slices
                    numSlices = length(obj.sliceIDs);
                    isSlicePresent = false(1, numSlices);
					
                    % Check for each tiff stack
                    for s=1:numSlices
						% If the ids are present, then has the tif been
                        % made?
                        isSlicePresent(s) = exist([obj.normalizedDataPath obj.mosaicPath 'slice_' num2str(s) '.tif'], 'file');

						% If the slice is not present check if the fovs are actually present (if not, then the slice was meant to be skipped)
						if ~isSlicePresent(s) & s>1
							if ~isempty(setdiff(obj.sliceIDs{s}, obj.fovIDs))
								isSlicePresent(s) = true;
							end
						end
                        
                    end
                    
                    % Return the status
                    complete(f) = all(isSlicePresent);
					
                otherwise
                    error('matlabFunctions:invalidArguments', 'The requested analysis type does not exist.');
            end
            
        end
    end
    
    % -------------------------------------------------------------------------
    % Combine and generate report for affine transforms
    % -------------------------------------------------------------------------
    function report = GenerateWarpReport(obj)
        % Generate a report on the fiducial warping
        
        % -------------------------------------------------------------------------
        % Map and load transforms and residuals if needed
        % -------------------------------------------------------------------------
        if isempty(obj.affineTransforms) && isempty(obj.residuals)
            
            % Display progress
            if obj.verbose
                PageBreak();
                display(['Searching for warping files in ' obj.normalizedDataPath obj.fiducialDataPath]);
                localTimer = tic;
            end

            % Map the location of affine transform and residuals
            affineFiles = BuildFileStructure(...
                [obj.normalizedDataPath obj.fiducialDataPath], ...
                'regExp', 'fov_(?<fov>[0-9]+)_affine', ...
                'fileExt', 'matb', ...
                'fieldNames', {'fov'}, ...
                'fieldConv', {@str2num});

            residualFiles = BuildFileStructure(...
                [obj.normalizedDataPath obj.fiducialDataPath], ...
                'regExp', 'fov_(?<fov>[0-9]+)_residual', ...
                'fileExt', 'matb', ...
                'fieldNames', {'fov'}, ...
                'fieldConv', {@str2num});

            if obj.verbose
                display(['Loading transformation files for ' num2str(length(affineFiles)) ' fov']);
            end
            
            % Check for existing files
            if isempty(affineFiles) || isempty(residualFiles)
                error('matlabFunctions:missingFiles', 'Could not find residual or affine files.');
            end

            % Confirm that an affine transform and a residual file exists for
            % each fovID
            if ~isempty(setdiff(obj.fovIDs, [affineFiles.fov]))
                error('matlabFunctions:missingFiles', 'Some affine files are missing.');
            end
            if ~isempty(setdiff(obj.fovIDs, [residualFiles.fov]))
                error('matlabFunctions:missingFiles', 'Some affine files are missing.');
            end

            % Load and combine transforms
            
            % Sort files in order of increasing fov id
            [~, sind] = sort([affineFiles.fov], 'ascend');
            affineFiles = affineFiles(sind);
            [~, sind] = sort([residualFiles.fov], 'ascend');
            residualFiles = residualFiles(sind);

            % Allocate memory
            obj.affineTransforms = repmat(affine2d(), [obj.numDataChannels obj.numFov]);
            obj.residuals = cell([obj.numDataChannels obj.numFov]);

            % Load and store affine files and residuals
            for i=1:length(affineFiles)
                obj.affineTransforms(:,i) = LoadByteStream(affineFiles(i).filePath, 'verbose', false);
                obj.residuals(:,i) = LoadByteStream(residualFiles(i).filePath, 'verbose', false);
            end
            
            % Display progress
            if obj.verbose
                display(['...completed in ' num2str(toc(localTimer)) ' s']);
            end
            
            % Save these new fields in the merfish decoder
            obj.UpdateField('affineTransforms', 'residuals');

            % Delete these files
            if obj.verbose
                display(['Removing all fiducial files...']);
                localTimer = tic;
            end
            
            rmdir([obj.normalizedDataPath obj.fiducialDataPath], 's'); % Delete folder
            
            if obj.verbose
                display(['...completed in ' num2str(toc(localTimer)) ' s']);
            end
            
        end
        
        if isempty(obj.geoTransformReport)
            % Update progress
            if obj.verbose
                display(['Preparing warp report']);
                localTimer = tic;
            end

            % -------------------------------------------------------------------------
            % Create warp report
            % -------------------------------------------------------------------------
            % Check to see if directory exists
            if ~exist([obj.normalizedDataPath obj.reportPath], 'dir');
                mkdir([obj.normalizedDataPath obj.reportPath]);
            end

            % Generate Transform Report and Save
            [report, figHandles] = GenerateGeoTransformReport(obj.affineTransforms,...
                obj.residuals, ...
                'edges', obj.parameters.warp.geoTransformEdges);

            % Update progress
            if obj.verbose
                display(['...completed in ' num2str(toc(localTimer)) ' s']);
            end

            % Save figure handles
            for h=1:length(figHandles)
                SaveFigure(figHandles(h), 'parameters', obj.parameters.display, ...
                    'savePath', [obj.normalizedDataPath obj.reportPath]);
                close(figHandles(h));
            end

            % Save report
            obj.geoTransformReport = report;
            obj.UpdateField('geoTransformReport');
        else
            display(['A warp report already exists. Clear geoTransformReport to allow recreation']);
        end

    end

    % -------------------------------------------------------------------------
    % Reset the scale factors used to normalize different data channels
    % -------------------------------------------------------------------------
    function ResetScaleFactors(obj)
        % Reset scale factors
        % obj.ResetScaleFactors();      % Reset scale factors
        
        % Reset scale factors
        obj.scaleFactors = [];
        if obj.verbose
            PageBreak();
            disp('Cleared existing scale factors');
        end
    end
    
    % -------------------------------------------------------------------------
    % Combine and generate report for affine transforms
    % -------------------------------------------------------------------------
    function OptimizeScaleFactors(obj, numFOV, varargin)
        % Optimize scale factors
        % OptimizeScaleFactors([]);      % Optimize on all FOV
        % OptimizeScaleFactors(numFOV);  % Optimize on numFOV randomly
        % selected FOV ids
        % OptimizeScaleFactors(numFOV, 'overwrite', true); % Overwrite
        % existing optimization

        % -------------------------------------------------------------------------
        % Handle varargin
        % -------------------------------------------------------------------------
        defaults = cell(0,3);
        defaults(end+1,:) = {'overwrite', 'boolean', false};    % Overwrite existing scale factors
        defaults(end+1,:) = {'useBlanks', 'boolean', true};     % Use blank barcodes in the optimization process or not
        defaults(end+1,:) = {'blankFunc', 'func', @(x)isempty(x.id)}; % Function to identify blanks based on codebook entry
        
        parameters = ParseVariableArguments(varargin, defaults, mfilename);
        
        % -------------------------------------------------------------------------
        % Clear existing scale factors if desired
        % -------------------------------------------------------------------------
        if parameters.overwrite
            obj.ResetScaleFactors();
        end
        
        % -------------------------------------------------------------------------
        % Check to see if optimization has already been completed
        % -------------------------------------------------------------------------
        if ~isempty(obj.scaleFactors)
            error('matlabFunctions:overwrite', 'Scale factors have been initialized. Clear before rerunning optimization');
        end
        
        % -------------------------------------------------------------------------
        % Check to see if the initial scale factors have been defined
        % -------------------------------------------------------------------------
        if isempty(obj.initScaleFactors)
            error('matlabFunctions:missingValues', 'Scale factors have not been initialized');
        end
        
        % -------------------------------------------------------------------------
        % Select fov IDs for optimization
        % -------------------------------------------------------------------------
        if isempty(numFOV)
            obj.optFovIDs = obj.fovIDs;
        elseif numFOV > obj.numFov
            warning('More fov ids were requested than are available. Using all fov ids');
            obj.optFovIDs = obj.fovIDs;
        else
            obj.optFovIDs = obj.fovIDs(randperm(obj.numFov, numFOV));
        end
        
        % Display fov IDs
        if obj.verbose
            display(['Using the following fov IDs for optimization:']);
            display(obj.optFovIDs);
        end

        % -------------------------------------------------------------------------
        % Generate decoding matrices
        % -------------------------------------------------------------------------
        % Generate the exact barcodes
        exactBarcodes = obj.GenerateDecodingMatrices();
        
        % Cut blanks (if requested)
        if ~parameters.useBlanks
            % Identify blanks
            isBlank = arrayfun(parameters.blankFunc, obj.codebook);
            
            % Remove blanks
            exactBarcodes = exactBarcodes(~isBlank, :);        
            
            if obj.verbose
                display(['Not using ' num2str(sum(isBlank)) ' blank barcodes for optimization']);
            end
        else
            if obj.verbose
                display(['Using all barcodes, including blanks, for optimization.']);
            end
        end
        
        % -------------------------------------------------------------------------
        % Initialize scale factors and other quantities to report
        % -------------------------------------------------------------------------
        numIter = obj.parameters.optimization.numIterOpt;
        localScaleFactors = zeros(numIter, obj.numBits);
        localScaleFactors(1,:) = obj.initScaleFactors;                        
        onBitIntensity = zeros(numIter-1, obj.numBits);
        allCounts = zeros(numIter-1, size(exactBarcodes,1));

        % -------------------------------------------------------------------------
        % Iterate
        % -------------------------------------------------------------------------
        for i=1:numIter
            % Display progress
            if obj.verbose
                PageBreak();
                display(['Starting ' num2str(i) ' of ' num2str(numIter) ' iterations']);
                display(['Utilizing ' num2str(obj.numPar) ' parallel workers']);
                iterTimer = tic;
            end
            
            % Loop over all files for analysis: perform in parallel to decrease time
            spmd (obj.numPar)
                % Create memory for accumulated pixel traces
                accumPixelTraces = zeros(size(exactBarcodes,1),  obj.numBits); % The accumulated pixel traces for all barcodes
                localCounts = zeros(1,size(exactBarcodes,1)); % The accumulated number of barcodes

                % Loop over files to analyze
                for f=labindex:numlabs:length(obj.optFovIDs)  
                    % Determine the local fovID
                    localFovID = obj.optFovIDs(f);
                    
                    % Display progress
                    if obj.verbose
                        displayStrings = {};
                        displayStrings{end+1} = ['Optimizing ' num2str(localFovID)];
                        displayStrings{end+1} = ['Loading stack and decoding'];
                        localTimer = tic;
                    end
                    
                    % Create tiff to read file name
                    tiffName2Read = [obj.normalizedDataPath obj.processedDataPath ...
                        'fov_' obj.fov2str(localFovID) '.tif'];
                    if ~exist(tiffName2Read, 'file')
                        error('matlabFunctions:missingFile', 'The requested tiff stack is not present.');
                    end

                    % Read (and low pass filter and crop) the tiff stack
                    localData = ReadAndFilterTiffStack(obj, tiffName2Read);
                    
                    % Decode data
                    [decodedImage, localMagnitude, pixelVectors] = obj.DecodePixels(localData, ...
                        localScaleFactors(i,:), exactBarcodes, obj.parameters.decoding.distanceThreshold);
                  
                    % Save memory by clearing local data
                    localData = [];
                    
                    % Set low intensity barcodes to zero
                    decodedImage(localMagnitude < obj.parameters.decoding.minBrightness) = 0;
                    
                    % Display progress
                    if obj.verbose
                        displayStrings{end+1} = ['...completed in ' num2str(toc(localTimer)) ' s'];
                        displayStrings{end+1} = ['Compiling pixel traces'];
                        localTimer = tic;
                    end
                    
                    % Loop over codebook entries
                    for b=1:size(exactBarcodes,1) 
                        % Define connected regions
                        conn = bwconncomp(decodedImage == b, obj.parameters.decoding.connectivity);

                        % Identify connected regions
                        properties = regionprops(conn, 'Area', 'PixelIdxList');

                        % Remove regions smaller than a minimum area
                        properties = properties([properties.Area] >= obj.parameters.decoding.minArea);

                        % Place additional cuts on properties for
                        % optimization 
                        properties = properties(...
                            [properties.Area] >= obj.parameters.optimization.areaThresh);
                                
                        % Accumulate the number of each barcode
                        localCounts(b) = localCounts(b) + length(properties);

                        % Accumulate the pixel traces
                        for l=1:length(properties)
                            localMag = localMagnitude(properties(l).PixelIdxList);
                            localTraces = pixelVectors(properties(l).PixelIdxList,:).* ...
                                repmat(localMag, [1 obj.numBits]);
                            localPixelTrace = mean(localTraces,1);
                            localPixelTrace = localPixelTrace/sqrt(sum(localPixelTrace.*localPixelTrace)); % Renormalize
                            accumPixelTraces(b,:) = accumPixelTraces(b,:) + localPixelTrace; % Accumulate
                        end                            
                    end % End loop over barcodes
                    
                    % Display progress
                    if obj.verbose
                        displayStrings{end+1} = ['...completed in ' num2str(toc(localTimer)) ' s'];
                        display(char(displayStrings)); % Display
                    end

                end % End loop over fov
            
            end % End spmd loop
            
            % ------------------------------------------------------------------------
            % Combine composite objects
            %--------------------------------------------------------------------------
            temp = localCounts{1};
            for j=2:length(localCounts)
                temp = temp + localCounts{j};
            end
            allCounts(i,:) = temp;

            temp = accumPixelTraces{1};
            for j=2:length(accumPixelTraces)
                temp = temp + accumPixelTraces{j};
            end
            accumPixelTraces = temp;
            
            % ------------------------------------------------------------------------
            % Compute new scale factors
            %--------------------------------------------------------------------------
            % Normalize and zero pixel traces
            normPixelTraces = accumPixelTraces./repmat(allCounts(i,:)', [1 obj.numBits]);
            normPixelTraces(exactBarcodes(:) == 0) = nan;

            % Compute the average intensity of the onBitIntensity
            onBitIntensity(i,:) = nanmean(normPixelTraces,1);
            refactors = onBitIntensity(i,:)./mean(onBitIntensity(i,:));

            % Record new scale factors
            if i<numIter
                localScaleFactors(i+1,:) = localScaleFactors(i,:).*refactors;
            end
            if obj.verbose
                display(['Completed iteration ' num2str(i) ' in ' num2str(toc(iterTimer)) ' s']);
            end

        end % End iteration loop
        
        % Save scale factors
        obj.scaleFactors = localScaleFactors(i,:);
        
        % Create and display optimization report
        obj.GenerateOptimizationReport(localScaleFactors, onBitIntensity, allCounts);
    end
    
    % -------------------------------------------------------------------------
    % Combine and generate report for affine transforms
    % -------------------------------------------------------------------------
    function SetScaleFactors(obj, scaleFactors)
        % Set scale factors manually
        % SetScaleFactors(scaleFactors);      % Optimize on all FOV

        % Check dimensions of provided scaleFactors
        dim = size(scaleFactors);
        
        if dim(1) ~= 1 || dim(2) ~= obj.numBits  || length(dim) > 2
            error('matlabFunctions:invalidInput', 'Provided scale factors are not the correct size');
        end
        
        % Set the scale factors
        obj.scaleFactors = scaleFactors;
        
    end
    
    % -------------------------------------------------------------------------
    % DecodeFOV
    % -------------------------------------------------------------------------
    function DecodeFOV(obj, fovIDs)
        % Decode individual images
        % DecodeFOV([]);                 % Decode all FOV
        % DecodeFOV(fovIDs);             % Decode specified fovIDs
                
        % -------------------------------------------------------------------------
        % Determine properties of the requested fov ids
        % -------------------------------------------------------------------------
        if isempty(fovIDs)
            fovIDs = obj.fovIDs;
        elseif ~all(ismember(fovIDs, obj.fovIDs))
            error('matlabFunctions:invalidArguments', 'An invalid fov id has been requested');
        end
        
        % -------------------------------------------------------------------------
        % Make directories if they do not exist
        % -------------------------------------------------------------------------
        % Base directory for barcodes
        if ~exist([obj.normalizedDataPath obj.barcodePath], 'dir')
            mkdir([obj.normalizedDataPath obj.barcodePath]);
        end
        % Directory for barcodes by fov
        barcodeByFovPath = [obj.normalizedDataPath obj.barcodePath filesep 'barcode_fov' filesep];
        if ~exist(barcodeByFovPath, 'dir')
            mkdir(barcodeByFovPath);
        end

        % -------------------------------------------------------------------------
        % Generate decoding matrices
        % -------------------------------------------------------------------------
        [exactBarcodes, singleBitErrorBarcodes] = obj.GenerateDecodingMatrices();
        
        % -------------------------------------------------------------------------
        % Initialize and check scale factors
        % -------------------------------------------------------------------------
        localScaleFactors = obj.scaleFactors; % Use the scale factors created by compiling the results of all optimized FOV
        if isempty(localScaleFactors)
            error('matlabFunctions:nonexistingVariable', 'Decoding cannot be run until scale factors have been initialized');
        end
        
        % -------------------------------------------------------------------------
        % Decode individual FOV
        % -------------------------------------------------------------------------
        spmd (obj.numPar)
            
            % Loop over individual fov
            for f=labindex:numlabs:length(fovIDs)
                % Determine local fovID
                localFovID = fovIDs(f);
                
                % Create display strings
                if obj.verbose
                    displayStrings = {};
                    displayStrings{end+1} = PageBreak('nodisplay');
                    displayStrings{end+1} = ['Started decoding of fov ' obj.fov2str(localFovID)];
                end
                
                % Define barcode file path and check for uncorrupted existance
                barcodeFile = [barcodeByFovPath 'fov_' obj.fov2str(localFovID) '_blist.bin'];
                
                % Erase if overwrite
                if obj.overwrite
                    if exist(barcodeFile, 'file')
                        delete(barcodeFile);
                        if obj.verbose; displayStrings{end+1} = ['Overwriting...']; end 
                    end
                end
                
                % Check if corrupt and erase if it is
                if exist(barcodeFile, 'file')
                    if obj.verbose;  displayStrings{end+1} = ['Found existing barcode file'];end 
                    isCorrupt = false;
                    try
                        fileHeader = ReadBinaryFileHeader(barcodeFile);
                        isCorrupt = fileHeader.isCorrupt;
                    catch
                        isCorrupt = true;
                    end
                    if isCorrupt
                        if obj.verbose; displayStrings{end+1} = ['File is corrupt. Overwriting...']; end; 
                        delete(barcodeFile);
                    end
                end
                
                % Skip if the file exists
                if exist(barcodeFile, 'file')
                    if obj.verbose; displayStrings{end+1} = ['File is complete. Skipping analysis']; end;
                    continue;
                end            

                % Display progress
                if obj.verbose
                    displayStrings{end+1} = ['Loading preprocessed stack'];
                    localTimer = tic;
                end

                % Create tiff to read file name
                tiffName2Read = [obj.normalizedDataPath obj.processedDataPath ...
                    'fov_' obj.fov2str(localFovID) '.tif'];
                if ~exist(tiffName2Read, 'file')
                    error('matlabFunctions:missingFile', 'The requested tiff stack is not present.');
                end
                
                % Read (and low pass filter) the tiff stack
                localData = ReadAndFilterTiffStack(obj, tiffName2Read);
                
                % Create display strings
                if obj.verbose
                    displayStrings{end+1} = ['... completed in ' num2str(toc(localTimer)) ' s'];
                    displayStrings{end+1} = ['Starting decoding'];
                    displayStrings{end+1} = ['... starting barcode assignment'];
                    localTimer = tic;
                    assignmentTimer = tic;
                end
                
                % Decode data
                [decodedImage, localMagnitude, pixelVectors, D] = obj.DecodePixels(localData, ...
                    localScaleFactors, exactBarcodes, obj.parameters.decoding.distanceThreshold);
                
                % Clear local data to open up memory
                localData = [];
                
                % Set low intensity barcodes to zero
                decodedImage(localMagnitude < obj.parameters.decoding.minBrightness) = 0;
                
                % Create display strings
                if obj.verbose
                    displayStrings{end+1} = ['... ... completed assignment in ' num2str(toc(assignmentTimer)) ' s'];
                    displayStrings{end+1} = ['... saving decoded image'];
                    saveTimer = tic;
                end
                
                % Save the decoded image and the magnitude image
                obj.SaveDecodedImageAndMagnitudeImage(decodedImage, reshape(localMagnitude, size(decodedImage)), localFovID);
                
                % Create display strings
                if obj.verbose
                    displayStrings{end+1} = ['... ... completed save in ' num2str(toc(saveTimer)) ' s'];
                    displayStrings{end+1} = ['... starting metadata assembly'];
                    metadataTimer = tic;
                end

                
                % Clear measured barcodes
                measuredBarcodes = [];
                for b=1:obj.numBarcodes % Loop over codebook entries
                    % Define connected regions
                    conn = bwconncomp(decodedImage == b, obj.parameters.decoding.connectivity);

                    % Identify connected regions
                    properties = regionprops(conn, 'Area', ...
                        'Centroid', 'PixelIdxList');

                    % Remove regions smaller than a minimum area
                    properties = properties([properties.Area] >= obj.parameters.decoding.minArea);
                    
                    % Compile properties of measured barcodes
                    measuredBarcodes = cat(2, measuredBarcodes, ...
                        obj.GenerateBarcodes(properties, localMagnitude, singleBitErrorBarcodes, ...
                        pixelVectors, D, b, localFovID));                        
                end % End loop over barcodes
                
                % Create display strings
                if obj.verbose
                    displayStrings{end+1} = ['... ... completed metadata assembly in ' num2str(toc(metadataTimer)) ' s'];
                    displayStrings{end+1} = ['... completed decoding in ' num2str(toc(localTimer)) ' s'];
                    displayStrings{end+1} = ['... saving ' num2str(length(measuredBarcodes)) ' barcodes'];
                    localTimer = tic;
                end
                
                % Write binary file for all measured barcodes
                barcodeFile = [barcodeByFovPath 'fov_' obj.fov2str(localFovID) '_blist.bin'];
                WriteBinaryFile(barcodeFile, measuredBarcodes);
                
                % Finish and display the progress strings
                if obj.verbose
                    displayStrings{end+1} = ['... completed in ' num2str(toc(localTimer)) ' s'];
                    display(char(displayStrings));
                end
                
            end
        end
    end
    
    % -------------------------------------------------------------------------
    % SaveDecodedImageAndMagnitudeImage
    % -------------------------------------------------------------------------
    function SaveDecodedImageAndMagnitudeImage(obj, decodedImage, localMagnitude, fovID)
        % Helper function: save these images as tiff stacks
        
        % Create path if necessary
        imagePath = [obj.normalizedDataPath obj.barcodePath filesep 'decoded_images' filesep];
        if ~exist(imagePath, 'dir')
            mkdir(imagePath);
        end
        
        % Create tiff Tag structure
        tiffTagStruct.ImageLength = size(decodedImage,1);
        tiffTagStruct.ImageWidth = size(decodedImage,2);
        tiffTagStruct.Photometric = Tiff.Photometric.MinIsBlack;
        tiffTagStruct.BitsPerSample = 32;
        tiffTagStruct.SamplesPerPixel = 1;
        tiffTagStruct.RowsPerStrip = 16;
        tiffTagStruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
        tiffTagStruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tiffTagStruct.Software = 'MATLAB';
        tiffTagStruct.ImageDescription = sprintf(['ImageJ=1.47a\n' ...      % ImageJ label
            'images=' num2str(obj.numZPos*2) '\n' ...                       % The total number of images
            'channels=1\n' ...                                              % The number of channels, 1 for now
            'slices=' num2str(obj.numZPos) '\n' ...                         % The number of z slices
            'frames=' num2str(2) '\n' ...                                   % The number of frames, the first is the decoded image the second is the magnitude image
            'hyperstack=true\n' ...
            'loop=false\n' ...
            ]);
        
        % Write tiff for the decodedImage
        tiffImage = Tiff([imagePath 'fov_' obj.fov2str(fovID) '.tif'], 'w8');
        
        % Write the decoded image
        for z=1:obj.numZPos
            tiffImage.setTag(tiffTagStruct);
            tiffImage.write(single(decodedImage(:,:,z)));
            tiffImage.writeDirectory();
        end
        % Write the magnitude image
        for z=1:obj.numZPos
            tiffImage.setTag(tiffTagStruct);
            tiffImage.write(single(localMagnitude(:,:,z)));
            if z~=obj.numZPos
                tiffImage.writeDirectory();
            end
        end
        % Close tiff files
        tiffImage.close();

    end
    
    % -------------------------------------------------------------------------
    % Convert pixel coordinates into absolute, real-world coordinates
    % -------------------------------------------------------------------------
    function absPos = Pixel2Abs(obj, pixelPos, fovID)
        % Convert a set of pixel coordinates in a given fov to the absolute coordinate system
        %
        % absPos = obj.Pixel2Abs(pixelPos, fovID)
        
        % Convert x and y
        absPos = obj.fovPos(obj.fovIDs == fovID,:) + ...
            obj.pixelSize/1000*obj.parameters.decoding.stageOrientation.* ...
            (double(pixelPos(:,1:2)) - obj.imageSize/2);
        
        % Convert z via linear interpolation
        if size(pixelPos,2) == 3
            absPos(:,3) = interp1(1:obj.numZPos, obj.zPos, double(pixelPos(:,3)));
        end
        
    end
    
    % -------------------------------------------------------------------------
    % GenerateBarcodes
    % -------------------------------------------------------------------------
    function measuredBarcodes = GenerateBarcodes(obj, properties, localMagnitude, singleBitErrorBarcodes, pixelTraces, D, b, fovID)
        
        % Define the basic barcode structure.... just for reference
        %         measuredBarcodes = repmat(struct(...
        %             'barcode', uint64(0), ...
        %             'barcode_id', uint16(0), ...   
        %             'fov_id', uint16(0), ...  
        %             'total_magnitude', single(0), ...
        %             'pixel_centroid', single(zeros(1, 2 + floor(obj.numZPos-1))), ...
        %             'weighted_pixel_centroid', single(zeros(1, 2 + floor(obj.numZPos-1))), ...
        %             'abs_position', single(zeros(1,3)), ...
        %             'area', uint16(0), ...
        %             'pixel_trace_mean', single(zeros(1, obj.numBits)), ...
        %             'pixel_trace_std', single(zeros(1, obj.numBits)), ...
        %             'is_exact', uint8(0), ...
        %             'error_bit', uint8(0), ...
        %             'error_dir', uint8(0), ...
        %             'av_distance', single(0)), ...
        %             [1 0]);

        measuredBarcodes = [];        
        % Loop over all barcodes
        for p=1:length(properties)
            % Transfer barcode properties
            measuredBarcodes(end+1).barcode = uint64(obj.codebook(b).barcode);                          % The barcode
            measuredBarcodes(end).barcode_id = uint16(b);                                           % The order of the entry in the codebook

            % Transfer fov id
            measuredBarcodes(end).fov_id = uint16(fovID);

            % Compute weighted centroid
            [y,x,z] = ind2sub([obj.imageSize obj.numZPos] - ...
                [2*obj.parameters.decoding.crop, 2*obj.parameters.decoding.crop, 0], ...
                properties(p).PixelIdxList);
            magnitude = localMagnitude(properties(p).PixelIdxList);
            totalMagnitude = sum(magnitude);
            weightedX = sum(x.*localMagnitude(properties(p).PixelIdxList))/totalMagnitude;
            weightedY = sum(y.*localMagnitude(properties(p).PixelIdxList))/totalMagnitude;
            weightedZ = sum(z.*localMagnitude(properties(p).PixelIdxList))/totalMagnitude;

            % Compute total magnitude
            measuredBarcodes(end).total_magnitude = single(totalMagnitude);                         % The total brightness of all pixel traces

            % Transfer position properties
            if length(properties(p).Centroid) == 2 % Handle single z slice case
                measuredBarcodes(end).pixel_centroid = uint16(properties(p).Centroid ...
                    + (obj.parameters.decoding.crop)*[1 1]);                                                     % The absolute pixel location in the image
                measuredBarcodes(end).weighted_pixel_centroid = single(...
                    [weightedX weightedY] ...
                    + (obj.parameters.decoding.crop)*[1 1]);                                                     % The X/Y/Z position in camera coordinates weighted by magnitude
            elseif length(properties(p).Centroid) == 3
                measuredBarcodes(end).pixel_centroid = uint16(properties(p).Centroid ...
                    + (obj.parameters.decoding.crop)*[1 1 0]);                                                     % The absolute pixel location in the image
                measuredBarcodes(end).weighted_pixel_centroid = single(...
                    [weightedX weightedY weightedZ] ...
                    + (obj.parameters.decoding.crop)*[1 1 0]);                                                     % The X/Y/Z position in camera coordinates weighted by magnitude
            end

            % Calculate absolute position
            absXYPos = obj.fovPos(obj.fovIDs == fovID,:) + ...
                obj.pixelSize/1000*obj.parameters.decoding.stageOrientation.* ...
                (measuredBarcodes(end).weighted_pixel_centroid(1:2) - obj.imageSize/2);

            weightedZ = round(weightedZ, 2); % Handle round-off error
            zLow = obj.zPos(floor(weightedZ)); % Handle weighted Z points inbetween z planes
            zHigh = obj.zPos(ceil(weightedZ));
            absZPos = zLow + (zHigh - zLow)*(weightedZ-floor(weightedZ));

            measuredBarcodes(end).abs_position = single([absXYPos absZPos]);                        % The absolute position in the stage coordinates

            % Transfer area
            measuredBarcodes(end).area = uint16(properties(p).Area);                                % The number of pixels

            % Compute average pixel traces
            localPixelTrace = mean(pixelTraces(properties(p).PixelIdxList,:),1);
            measuredBarcodes(end).pixel_trace_mean = single(localPixelTrace);
            measuredBarcodes(end).pixel_trace_std = single(std(pixelTraces(properties(p).PixelIdxList,:),1,1));

            % Compute location within the Hamming sphere
            nnID = knnsearch(squeeze(singleBitErrorBarcodes(b,:,:)), localPixelTrace,'K',1);
            measuredBarcodes(end).is_exact = uint8(nnID == 1);                                       % Whether or not this barcode best matches the exact barcode
            measuredBarcodes(end).error_bit = uint8(nnID-1);                                         % The bit at which an error occurred (0 if no error occurred)
            if measuredBarcodes(end).error_bit > 0
                measuredBarcodes(end).error_dir = ...
                    uint8(singleBitErrorBarcodes(b,1,measuredBarcodes(end).error_bit)>0);            % The identity of the original bit (1 = 1->0 error; 0 = 0->1 error)
            else
                measuredBarcodes(end).error_dir = uint8(0);
            end

            % Record average distance
            measuredBarcodes(end).av_distance = single(mean(D(properties(p).PixelIdxList)));         % The average distance from the pixel traces to the matched barcode

        end
    end
    
    % -------------------------------------------------------------------------
    % Combine pixel histograms
    % -------------------------------------------------------------------------
    function localData = ReadAndFilterTiffStack(obj, tiffName2Read)
        % Open and read a tiff stack
        
        % Open tiff files
        tiffToRead = Tiff(tiffName2Read, 'r');

        % Allocate memory for image
        localData = zeros(obj.imageSize(2), obj.imageSize(1), obj.numZPos, obj.numBits, 'double');

        % Loop over frames, loading, accumulating pixel values, and
        % filtering (if requested)
        for b=1:obj.numBits
            % Loop over z-stacks
            for z=1:obj.numZPos
                % Determine frame location in stack
                frame = (b-1)*obj.numZPos + z;
                % Set directory and load
                tiffToRead.setDirectory(frame);
                
                %localTimer = tic; %% TEMPORARY CODE
                localFrame = tiffToRead.read(); 
                %display(['Loaded frame ' num2str(frame) ' in ' num2str(toc(localTimer))]); %% TEMPORARY CODE
                %localTimer = tic; %% TEMPORARY CODE
                
                % Low pass filter
                if obj.parameters.decoding.lowPassKernelSize > 0
                    localData(:,:,z,b) = imgaussfilt(localFrame, obj.parameters.decoding.lowPassKernelSize);
                else
                    localData(:,:,z,b) = localFrame;
                end
                %display(['Filtered frame ' num2str(frame) ' in ' num2str(toc(localTimer))]); %% TEMPORARY CODE

            end
        end

        % Crop edges
        localData = localData((obj.parameters.decoding.crop+1):(end-obj.parameters.decoding.crop), ...
            (obj.parameters.decoding.crop+1):(end-obj.parameters.decoding.crop), :, :);

        % Close tiff
        tiffToRead.close();
    end    
    % -------------------------------------------------------------------------
    % Initialize scale factors
    % -------------------------------------------------------------------------
    function InitializeScaleFactors(obj)
        % Initialize scale factors and combine pixel histograms
        
        % Combine pixel histograms
        if isempty(obj.pixelHistograms)
            % Display progress
            if obj.verbose
                PageBreak();
                display('Combining pixel histograms');
                combineTimer = tic;
            end

            % Identify the saved pixel histogram files
            foundFiles = BuildFileStructure([obj.normalizedDataPath obj.processedDataPath 'pixel_histograms' filesep], ...
                'regExp', 'pixel_data_fov_(?<fov>[0-9]+)', ...
                'fileExt', 'matb', ...
                'fieldNames', {'fov'}, ...
                'fieldConv', {@str2num});
            fileType = 'matb';
            
            % Allow for csv files (and different naming convention)
            if isempty(foundFiles)
                foundFiles = BuildFileStructure([obj.normalizedDataPath obj.processedDataPath 'pixel_histograms' filesep], ...
                    'regExp', 'fov_(?<fov>[0-9]+)', ...
                    'fileExt', 'csv', ...
                    'fieldNames', {'fov'}, ...
                    'fieldConv', {@str2num});
                fileType = 'csv';
            end

            % Display progress
            if obj.verbose
                display(['...found ' num2str(length(foundFiles)) ' pixel histogram files']);
            end

            % Check that all fov are present
            if ~isempty(setdiff([obj.fovIDs], [foundFiles.fov]))
                display(['The following fov ids do not have pixel_histograms:']);
                display(setdiff([obj.fovIDs], [foundFiles.fov]));
                error('matlabFunctions:missingFiles', 'A pixel histogram for all fov ids was not found.');
            end

            % Initialize and allocate memory
            pixelHistograms =  zeros(obj.numBits, uint16(inf)); % Allocate memory and initialize to zero

            % Loop over all fov ids
            for f=1:obj.numFov
                localFovID = obj.fovIDs(f);
                switch fileType
                    case 'matb'
                        localData = LoadByteStream(foundFiles([foundFiles.fov] == localFovID).filePath, 'verbose', false);
                    case 'csv'
                        localData = csvread(foundFiles([foundFiles.fov] == localFovID).filePath);
                    otherwise
                        error('matlabFunctions:invalidArguments', 'Unrecognized file extension requested for pixel histograms');
                end
                pixelHistograms = pixelHistograms + localData; % Combine pixel histogram
            end
            
            % Save pixel histograms
            obj.UpdateField('pixelHistograms');

            % Update progress
            if obj.verbose
                display(['...completed in ' num2str(toc(combineTimer)) ' s']);
            end
            
            % Remove pixel histogram files
            if obj.verbose
                display(['Deleting intermediate files...']);
                localTimer = tic;
            end
            rmdir([obj.normalizedDataPath obj.processedDataPath 'pixel_histograms' filesep], 's');
            if obj.verbose
                display(['...completed in ' num2str(toc(localTimer)) ' s']);
            end
            
            % Update pixel histograms in decoder
            obj.pixelHistograms = pixelHistograms;
        end
        
        % Create report path for saving the initial scale factor selection
        if ~exist([obj.normalizedDataPath obj.reportPath], 'dir')
            mkdir([obj.normalizedDataPath obj.reportPath]);
        end
        
        % Create figure to display selections
        figHandle = figure('Name', ['Intensity normalization'],...
            'Color', 'w', ...
            'visible', obj.parameters.display.visibleOption);

        markers = {'r', 'g', 'b', 'y', 'm', 'c',...
            'r--', 'g--', 'b--', 'y--', 'm--', 'c--', ...
            'r.-', 'g.-', 'b.-', 'y.-', 'm.-', 'c.-'};

        % Loop over bits
        minValue = inf;
        for i=1:obj.numBits
            cumSum = cumsum(obj.pixelHistograms(i,1:end));
            cumSum = cumSum/cumSum(end); % Normalize to 1

            [~, localScale] = min(abs(cumSum - obj.parameters.optimization.quantileTarget));
            initScaleFactors(i) = localScale + 1; % Add 1 b/c 0 is the first bin in the histogram

            semilogx(cumSum, markers{mod(i-1, length(markers))+1}); hold on;

            % Record minimum value for improved display
            minValue = min([minValue cumSum(1)]);
            quantilePos(i) = cumSum(initScaleFactors(i));
        end
        
        % Plot selection points
        for i=1:obj.numBits
            semilogx([1 1]*initScaleFactors(i) , [0 1]*quantilePos(i), 'k--');
        end

        % Add plot annotations
        xlabel('Intensity');
        ylabel('Cumulative probability');
        legend(cellstr(num2str([1:obj.numBits]')), 'Location', 'SouthEast');
        ylim([minValue 1]);
        
        % Save and close report figure
        SaveFigure(figHandle, 'parameters', obj.parameters.display, ...
            'savePath', [obj.normalizedDataPath obj.reportPath]);
        close(figHandle);
        
        % Normalize the scale factors to 1 to facilitate real brightness
        % measures
        if obj.parameters.optimization.normalizeToOne
            if obj.verbose
                disp('Normalizing scale factors to 1');
            end
            initScaleFactors = initScaleFactors/mean(initScaleFactors);
        end
        
        % Record the initial scale factors
        obj.initScaleFactors = initScaleFactors;
        obj.UpdateField('initScaleFactors');
        
        % Display progress
        if obj.verbose
            display(['Using the following scale factors: ']);
            display([num2str(obj.initScaleFactors)]);
        end

    end

    % -------------------------------------------------------------------------
    % GenerateDecodingMatrices
    % -------------------------------------------------------------------------
    function [weightedBarcodes, singleBitErrorBarcodes] = GenerateDecodingMatrices(obj)
        % Generate matrics used in the decoding process from the loaded
        % codebook
        
        % Display progress
        if obj.verbose
            PageBreak();
            display('Generating decoding matrices');
            localTimer = tic;
        end
        
        % Extract the binary barcodes
        binaryBarcodes = de2bi([obj.codebook.barcode], obj.numBits);
        
        % Calculate magnitude and normalize to unit length
        magnitudeTemp = sqrt(sum(binaryBarcodes.*binaryBarcodes,2));
        weightedBarcodes = binaryBarcodes./repmat(magnitudeTemp, [1 size(binaryBarcodes,2)]);

        % Create useful function to generate bit flips
        bitFlip = @(x,n)bitset(x, n, ~bitget(x,n)); 

        % Compute the single bit error matrices for computing error rates
        for b=1:obj.numBarcodes
            % Create an array of the correct and all single bit flip barcodes
            singleBitFlipIntegers = [obj.codebook(b).barcode ...
                arrayfun(@(n)bitFlip(obj.codebook(b).barcode,n), 1:obj.numBits)];
            singleBitFlipBinary = de2bi(singleBitFlipIntegers,obj.numBits);

            % Normalize
            localMag = sqrt(sum(singleBitFlipBinary.*singleBitFlipBinary,2));
            singleBitErrorBarcodes(b,:,:) = singleBitFlipBinary./repmat(localMag, [1 size(singleBitFlipBinary,2)]);
        end      
        
        if obj.verbose
            PageBreak();
            display(['...completed in ' num2str(toc(localTimer)) ' s']);
        end
    end
                
    % -------------------------------------------------------------------------
    % SetParallel
    % -------------------------------------------------------------------------
    function GenerateOptimizationReport(obj, localScaleFactors, onBitIntensity, allCounts)
        % Check to see if directory exists
        if ~exist([obj.normalizedDataPath obj.reportPath], 'dir')
            mkdir([obj.normalizedDataPath obj.reportPath]);
        end
                
        % Display report
        figHandle = figure('Name', ['Optimization report'], ...
            'Color', 'w', 'visible', obj.parameters.display.visibleOption);

        % Compile data to plot
        data = {localScaleFactors, onBitIntensity, log10(allCounts)};
        titles = {'Scale factor', 'On-bit intensity', 'Counts (log_{10})'};
        
        % Create subplots
        for p=1:length(data)
            subplot(1, length(data), p);
            pcolor(data{p}); hold on;
            shading flat;
            ylabel('Iteration');
            title(titles{p});
            xlim([1 size(data{p},2)]);
            ylim([1 size(data{p},1)]);
            colorbar;
        end

        % Save and close report figure
        SaveFigure(figHandle, 'parameters', obj.parameters.display, ...
            'savePath', [obj.normalizedDataPath obj.reportPath]);
        close(figHandle);
    end
    
    % -------------------------------------------------------------------------
    % SetParallel
    % -------------------------------------------------------------------------
    function SetParallel(obj, p)
        % Set or update the parallel.Pool object
        % obj.setParallel(p) 
        % obj.setParallel([]) removes the existing pool
        
        % -------------------------------------------------------------------------
        % Check validity
        % -------------------------------------------------------------------------
        if ~isempty(p) && ~isa(p, 'parallel.Pool')
            error('matlabFunctions:invalidArgument', 'Must provide a valid parallel.Pool object');
        end
        
        obj.parallel = p;
        if isempty(obj.parallel)
            obj.numPar = 0;
        else
            obj.numPar = obj.parallel.NumWorkers;
        end
    end
    
    % -------------------------------------------------------------------------
    % Save Function
    % -------------------------------------------------------------------------
    function Save(obj, varargin)
        % Save the MERFISHDesigner object in a directory specified by dirPath
        % obj.Save(dirPath)
        
        % -------------------------------------------------------------------------
        % Check directory validity
        % -------------------------------------------------------------------------
        % Set default save path
        if isempty(varargin)
            dirPath = [obj.normalizedDataPath obj.mDecoderPath]; % Make path an absolute path
        else
            dirPath = varargin{1};
        end
        
        % Coerce provided path
        if dirPath(end) ~= filesep
            dirPath(end+1) = filesep;
        end
        
        % Make the directory if needed
        status= mkdir(dirPath);
        if ~status
            error('matlabFunctions:invalidArguments', 'Invalid directory path');
        end
        
        % Update MERFISH decoder path
        % obj.mDecoderPath = dirPath;  % This will make the path an absolute rather than a relative path...
        
        % -------------------------------------------------------------------------
        % Define fields to save
        % -------------------------------------------------------------------------
        fieldsToSave = properties(obj);
        fieldsToSave = setdiff(fieldsToSave, {'parallel', 'numPar'});
        
        % -------------------------------------------------------------------------
        % Save fields
        % -------------------------------------------------------------------------
        for i=1:length(fieldsToSave)
            switch fieldsToSave{i}
                case '' % Reserved for future use
                    
                otherwise
                    SaveAsByteStream([dirPath fieldsToSave{i} '.matb'], ...
                        obj.(fieldsToSave{i}), 'verbose', obj.verbose);
            end
        end
    end
    
    % -------------------------------------------------------------------------
    % SetParameter
    % -------------------------------------------------------------------------
    function SetParameter(obj, varargin)
        % Set fields in the parameters structure
        
        % Define default parameter sets
        defaultParameters = cell(0,2);
        defaultParameters(end+1,:) = {'warp', MERFISHDecoder.DefaultWarpParameters()};
        defaultParameters(end+1,:) = {'preprocess', MERFISHDecoder.DefaultPreprocessingParameters()};
        defaultParameters(end+1,:) = {'decoding', MERFISHDecoder.DefaultDecodingParameters()};
        defaultParameters(end+1,:) = {'optimization', MERFISHDecoder.DefaultOptimizationParameters()};
        defaultParameters(end+1,:) = {'display', MERFISHDecoder.DefaultDisplayParameters()};
        defaultParameters(end+1,:) = {'segmentation', MERFISHDecoder.DefaultSegmentationParameters()};
        defaultParameters(end+1,:) = {'quantification', MERFISHDecoder.DefaultQuantificationParameters()};
        defaultParameters(end+1,:) = {'summation', MERFISHDecoder.DefaultSummationParameters()};
        defaultParameters(end+1,:) = {'molecules', MERFISHDecoder.DefaultMoleculeParameters()};
        
        % Create flag to catch parameters not updated
        requestedFields = varargin(1:2:end);
        
        % Loop over parameter sets
        for p=1:size(defaultParameters,1)
            % Find matches to defaults
            localDefaults = defaultParameters{p,2};
            matchInd = find(ismember(varargin(1:2:end), localDefaults(:,1)));
            matchInd = sort([(2*matchInd -1) 2*matchInd]);
            if ~isempty(matchInd)
                parameters = ParseVariableArguments(varargin(matchInd), localDefaults, mfilename);
                % Transfer fields
                fieldsToUpdate = varargin(matchInd(1:2:end));
                for f=1:length(fieldsToUpdate)
                    % Store parameters
                    obj.parameters.(defaultParameters{p,1}).(fieldsToUpdate{f}) = parameters.(fieldsToUpdate{f});
                end
                % Mark fields as updated
                requestedFields = setdiff(requestedFields, fieldsToUpdate);
            end
        end
        
        % Raise warning
        if ~isempty(requestedFields)
            warning('matlabFunctions:invalidParameters', 'One or more of the requested parameters are not recognized');
            for r=1:length(requestedFields)
                disp(requestedFields{r});
            end
        end
        
    end

    
    % -------------------------------------------------------------------------
    % Update normalized data path
    % -------------------------------------------------------------------------
    function UpdateNormalizedDataPath(obj,newPath)
        % Update the normalized data path. Useful when the directory has
        % been copied to a new location
        % obj.UpdateNormalizedDataPath(newPath);
        
        % -------------------------------------------------------------------------
        % Check validity of arguments
        % -------------------------------------------------------------------------
        if nargin < 1 | ~exist(newPath, 'dir')
            error('matlabFunctions:invalidPath', 'A valid path must be provided');
        end
        
        % -------------------------------------------------------------------------
        % Update path
        % -------------------------------------------------------------------------
        obj.normalizedDataPath = newPath;
        
        % -------------------------------------------------------------------------
        % Handle filesep conversion for all relative paths
        % -------------------------------------------------------------------------
        % Define all relative paths
        relativePaths = {'mDecoderPath', 'fiducialDataPath', 'warpedDataPath', ...
            'processedDataPath', 'barcodePath', 'reportPath', ...
            'segmentationPath', 'summationPath', 'mosaicPath'};
        
        % Loop over relative paths
        for p=1:length(relativePaths)
            tempPath = obj.(relativePaths{p});
            tempPath(end) = filesep;
            obj.(relativePaths{p}) = tempPath;
        end
        
    end

    % -------------------------------------------------------------------------
    % Update Field
    % -------------------------------------------------------------------------
    function UpdateField(obj, varargin)
        % Save the MERFISHDesigner object in a directory specified by dirPath
        % obj.UpdateField(fields...)
                        
        % -------------------------------------------------------------------------
        % UpdateFields
        % -------------------------------------------------------------------------
        validFields = properties(obj);
        
        % -------------------------------------------------------------------------
        % Save fields
        % -------------------------------------------------------------------------
        for i=1:length(varargin)
            if ismember(varargin{i}, validFields)
                switch varargin{i}
                    case '' % Reserved for future use
                    otherwise
                        if obj.verbose
                            display(['Updating ' varargin{i}]);
                        end
                        SaveAsByteStream([obj.normalizedDataPath obj.mDecoderPath varargin{i} '.matb'], ...
                            obj.(varargin{i}), 'verbose', obj.verbose);
                end
            else
                warning('%s is not a valid field', varargin{i});
            end
        end
    end

    % -------------------------------------------------------------------------
    % LoadField
    % -------------------------------------------------------------------------
    function LoadField(obj, varargin)
        % Update the current MERFISHDecoder with revised values saved to
        % disk
        % obj.LoadField(fields...)
                        
        % -------------------------------------------------------------------------
        % UpdateFields
        % -------------------------------------------------------------------------
        validFields = properties(obj);
        
        % -------------------------------------------------------------------------
        % Save fields
        % -------------------------------------------------------------------------
        for i=1:length(varargin)
            if ismember(varargin{i}, validFields)
                switch varargin{i}
                    case '' % Reserved for future use
                    otherwise
                        if obj.verbose
                            display(['Updating ' varargin{i}]);
                        end
                        obj.(varargin{i}) = LoadByteStream([obj.normalizedDataPath obj.mDecoderPath varargin{i} '.matb'], ...
                            'verbose', obj.verbose);
                end
            else
                warning('%s is not a valid field', varargin{i});
            end
        end
    end

    % -------------------------------------------------------------------------
    % Downsample dataset: This hidden function allows datasets to be
    % downsampled for the purpose of illustration
    % -------------------------------------------------------------------------
    function Downsample(obj, bitNamesToRemove, fovIDsToRemove)
        % Identify the bits to remove
        
        shouldKeepBit = ~ismember({obj.dataOrganization.bitName}, bitNamesToRemove);
        shouldKeepFOV = ~ismember(obj.fovIDs, fovIDsToRemove);
        
        % Cut the data organization file
        obj.dataOrganization = obj.dataOrganization(shouldKeepBit);
        
        % Cut the nubmer of data channels
        obj.numDataChannels = length(obj.dataOrganization);
        
        % Save the original largest fovID for the purposes of matching the
        % pad in fov strings
        obj.originalMaxFovID = max(obj.fovIDs);
        
        % Cut the fov
        obj.fovIDs = obj.fovIDs(shouldKeepFOV);
        obj.fovPos = obj.fovPos(shouldKeepFOV,:);
        obj.numFov = length(obj.fovIDs);
        
        % Clear values associated with optimization
        obj.scaleFactors = [];
        obj.optFovIDs = [];
        obj.rawDataFiles = [];
        obj.affineTransforms = obj.affineTransforms(shouldKeepBit, shouldKeepFOV);
        obj.residuals = obj.residuals(shouldKeepBit, shouldKeepFOV);
        obj.geoTransformReport = [];
        obj.pixelHistograms = [];
        obj.initScaleFactors = [];
        
    end
   
end
% -------------------------------------------------------------------------
% Static methods
% -------------------------------------------------------------------------
methods (Static)
    % -------------------------------------------------------------------------
    % Build a TRDesigner object from a saved version
    % -------------------------------------------------------------------------
    function obj = Load(dirPath, varargin)
        % obj = MERFISHDecoder.Load(dirPath)

        % -------------------------------------------------------------------------
        % Check provided path
        % -------------------------------------------------------------------------
        if dirPath(end) ~= filesep
            dirPath(end+1) = filesep;
        end
        if ~isfolder(dirPath)
            error('matlabFunctions:invalidArguments', 'The provided path is not valid');
        end
        
        % -------------------------------------------------------------------------
        % Handle varargin
        % -------------------------------------------------------------------------
        defaults = cell(0,3);
        defaults(end+1,:) = {'verbose', 'boolean', true};   % Display load progress
        parameters = ParseVariableArguments(varargin, defaults, mfilename);
        
        % -------------------------------------------------------------------------
        % Create empty object (to define fields to load)
        % -------------------------------------------------------------------------
        obj = MERFISHDecoder();
        
        % -------------------------------------------------------------------------
        % Define fields to load
        % -------------------------------------------------------------------------
        fieldsToLoad = setdiff(properties(obj), {'parallel', 'numPar', 'mDecoderPath', ...
            'fov2str'}); %% Update me
        
        % -------------------------------------------------------------------------
        % Display progress
        % -------------------------------------------------------------------------
        if parameters.verbose
            PageBreak();
            display(['Loading MERFISH Decoder from ' dirPath]);
            loadTimer = tic;
        end
        
        % -------------------------------------------------------------------------
        % Check to see if valid -- all previous versions have verbose field
        % -------------------------------------------------------------------------
        if ~exist([dirPath obj.mDecoderPath 'verbose.matb'])
            error('The mDecoder appears to be corrupt!');
        end
        % -------------------------------------------------------------------------
        % Load the version number
        % -------------------------------------------------------------------------
        try
            version = LoadByteStream([dirPath obj.mDecoderPath 'version.matb'], 'verbose', false);
        catch % Handle the case that the saved version is not in matb format or is non-existent
            version = '0.1';
        end
        
        % -------------------------------------------------------------------------
        % Load properties/data
        % -------------------------------------------------------------------------
        switch version
            case {'0.1', '0.2', '0.3','0.4', '0.5', '0.6'}
                for i=1:length(fieldsToLoad)
                    switch fieldsToLoad{i}
                        case '' % Reserved for future use
                        otherwise
                            if exist([dirPath obj.mDecoderPath fieldsToLoad{i} '.matb'], 'file')
                                obj.(fieldsToLoad{i}) = LoadByteStream([dirPath obj.mDecoderPath fieldsToLoad{i} '.matb'], ...
                                'verbose', false);
                            else
                                warning('Did not find a default field: %s', fieldsToLoad{i});
                            end
                    end
                end
            otherwise
                error('matlabFunctions:unsupportedVersion', 'The version is not supported');
        end
    
        if parameters.verbose
            display(['Version: ' obj.version]);
            display(['FOV: ' num2str(obj.numFov)]);
            display(['Bits: ' num2str(obj.numBits)]);
            display(['Data channels: ' num2str(obj.numDataChannels)]);
            display(['Number of cameras: ' num2str(obj.numCameraIDs)]);
            display(['Loaded in ' num2str(toc(loadTimer)) ' s']);
        end
        
        % -------------------------------------------------------------------------
        % Handle updated parameters fields (version up-conversion)
        % -------------------------------------------------------------------------
        % Define default parameter sets
        [~, defaultParameters] = obj.InitializeParameters();
        
        % Compare to defaults to loaded parameters
        for c=1:size(defaultParameters,1) % Loop over parameters sub sets
            % Get field type
            parametersType = defaultParameters{c,1};
            % Convert to a structure
            parametersStruct = ParseVariableArguments({}, defaultParameters{c,2});
            % Check to see if the loaded parameters have this sub-set
            if ~isfield(obj.parameters, parametersType) 
                obj.parameters.(parametersType) = parametersStruct;
                display(['Set missing ' parametersType ' parameters to the default values.']);
            else % Investigate fields one by one with this sub-set (if it exists)
                % Get parameter names
                localParameterFields = fieldnames(parametersStruct);
                % Loop over parameter names
                for f=1:length(localParameterFields)
                    % If the parameter is missing, add it and give it the
                    % default value
                    if ~isfield(obj.parameters.(parametersType), localParameterFields{f})
                        obj.parameters.(parametersType).(localParameterFields{f}) = ...
                            parametersStruct.(localParameterFields{f});
                        display(['Set missing value ' parametersType '.' localParameterFields{f} ' to the default value.']);
                    end
                end
            end     
        end
        % -------------------------------------------------------------------------
        % Handle up-conversion of rawDataFiles
        % -------------------------------------------------------------------------
        if ~isempty(obj.rawDataFiles) && ~isfield(obj.rawDataFiles(1), 'cameraID')
            disp(['Adding a cameraID field to rawDataFiles to maintain version compatibility']);
            for F=1:length(obj.rawDataFiles)
                obj.rawDataFiles(F).cameraID = '';
            end
        end
        
        % -------------------------------------------------------------------------
        % Handle updated parameters fields (important for some version up-conversions)
        % -------------------------------------------------------------------------
        padNum2str = @(x,y)num2str(x,['%0',num2str(ceil(log10(y+1))),'d']);
        if ~isempty(obj.originalMaxFovID)
            obj.fov2str = @(x)padNum2str(x, obj.originalMaxFovID);
        else
            obj.fov2str = @(x)padNum2str(x, max(obj.fovIDs));
        end
        % -------------------------------------------------------------------------
        % Handle the location of the mDecoder (the normalizedDataPath is
        % the load path)
        % -------------------------------------------------------------------------
        obj.UpdateNormalizedDataPath(dirPath);
    end
    
    % -------------------------------------------------------------------------
    % Generate Warp defaults
    % -------------------------------------------------------------------------
    function defaults = DefaultWarpParameters()
        % Generate the default parameters information for the warping method
        % defaultCell = MERFISHDecoder.WarpDefaultParameters();
    
        defaults = cell(0,3);
        defaults(end+1,:) = {'warpingDataPath', ...             % Path to saved data for fiducial information
            'filePath', []};
        defaults(end+1,:) = {'fiducialFitMethod', ...           % Method for fiducial fitting
            {'daoSTORM'}, ...
            'daoSTORM'};
        defaults(end+1,:) = {'controlPointOffsetRange', ...
            'array', [-60:.5:60]};                              % The histogram properties for determining crude offset
        defaults(end+1,:) = {'numNN', 'positive', 10};          % The number of nearest neighbors to include in the search
        defaults(end+1,:) = {'pairDistanceTolerance', 'positive', 3}; % The multiple of the histogram distance used to judge paired beads
        defaults(end+1,:) = {'pixelSize', 'positive', 109};     % Pixel size in nm/pixel
        defaults(end+1,:) = {'sigmaInit', 'positive', 1.6};     % Initial guess for PSF size for fiducial images
        defaults(end+1,:) = {'daoThreshold', 'positive', 500};  % The minimum brightness for fitting fiducials
        defaults(end+1,:) = {'daoBaseline', 'positive', 100};   % The assumed baseline of the camera
        
        defaults(end+1,:) = {'exportWarpedBeads', 'boolean', true}; % Export a warped set of bead images
        
        defaults(end+1,:) = {'cameraOrientation', 'array', ...  % Set the orientation of the camera to the stage:
            [0 0 0]};                                           % The first/second element invert X/Y if '1'; the third element transposes X/Y
        
        defaults(end+1,:) = {'geoTransformEdges', 'cell', ...   % The histogram bins for position calculating position dependent bias in warping
            {0:25:2048, 0:25:2048}}; 
        
        defaults(end+1,:) = {'colorTransforms', 'struct', ...   % A structure array that includes affine transforms for all color channels. 
            struct('color', '', 'transform', [], 'type', '')};  % The color entries must match those provided in the data organization file
        
    end
    
    % -------------------------------------------------------------------------
    % Generate Preprocessing Default Parameters
    % -------------------------------------------------------------------------
    function defaults = DefaultPreprocessingParameters()
        % Generate the default parameters information for the preprocessing method
        % defaultCell = MERFISHDecoder.PreprocessingDefaultParameters();
    
        defaults = cell(0,3);
        % Select the preprocessing method
        defaults(end+1,:) = {'preprocessingMethod', ...         % The method for preprocessing
            {'highPassDecon', 'highPassErosion', ...
            'highPassDeconWB'}, ...
            'highPassDecon'};
        
        % Parameters for high pass
        defaults(end+1,:) = {'highPassKernelSize', ...          % The size of the Gaussian filter used for high pass filtering
            'positive', 3};       
        
        % Parameters for deconvolution
        defaults(end+1,:) = {'deconKernel', ...                 % The kernel to use for Lucy Richardson Deconvolution (no decon applied if empty)
            'array', fspecial('gaussian', 10, 2)};
        defaults(end+1,:) = {'numIterDecon', 'nonnegative', 20};
        
        % Parameters for GPU deconvolution: DEPRECATED FOR NOW                      
%         defaults(end+1,:) = {'deconGPUsigma', 'nonnegative', 2};  % The sigma for the Lucy Richardson kernel for the GPU decon
%         defaults(end+1,:) = {'deconGPUkernel', 'nonnegative', 10}; % The kernel size for the Lucy Richardson kernel for the GPU decon
        
        % Parameters for erosion
        defaults(end+1,:) = {'erosionElement', ...              % The morphological structuring element used for image erosion 
            'array', strel('disk', 1)};
        
    end
    
    % -------------------------------------------------------------------------
    % Generate Decoding Default Parameters
    % -------------------------------------------------------------------------
    function defaults = DefaultDecodingParameters()
        % Generate the default parameters for the decoding process
   
        % Create empty cell
        defaults = cell(0,3);
        
        % Parameters for additional preprocessing
        defaults(end+1,:) = {'lowPassKernelSize', ...           % The size of kernel for an intial low pass average (0 indicates on averaging)
            'nonnegative', 1};
        defaults(end+1,:) = {'crop', 'nonnegative', 40};        % The number of pixels to crop from each edge of the image

        % Parameters for decoding
        defaults(end+1,:) = {'decodingMethod', ...              % The method to employ for decoding
            {'distanceHS1'}, 'distanceHS1'};
        defaults(end+1,:) = {'distanceThreshold', ...           % The distance defining the Hamming Sphere around each barcode to assign to that barcode
            'positive', 0.5176};
        
        % Parameters for pre-saving cuts
        defaults(end+1,:) = {'minBrightness', ...               % The minimum brightness to call a pixel as a barcode
            'nonnegative', 10^0};
        defaults(end+1,:) = {'minArea', ...                     % The minimum area to save a barcode
            'nonnegative', 1};
        defaults(end+1,:) = {'connectivity', ...                % The connectivity matrix used to connect objects: 
            'array', cat(3, ...                                 % The default is to NOT connect objects between z-planes
            zeros(3,3), conndef(2, 'maximal'), zeros(3,3))};

        % Parameters for absolute coordinates
        defaults(end+1,:) = {'stageOrientation', ...            % The orientation of the different stage axis wrt to the camera
            'array', [1 1]};
        defaults(end+1,:) = {'pixelSize', 'positive', 109};     % Pixel size in nm/pixel
        
        
    end
    
    % -------------------------------------------------------------------------
    % Generate Segmentation Default Parameters
    % -------------------------------------------------------------------------
    function defaults = DefaultSegmentationParameters()
        % Generate the default parameters for the feature segmentation process
   
        % Create empty cell
        defaults = cell(0,3);
        
        % Parameters for determining segmentation method
        defaults(end+1,:) = {'segmentationMethod', ...          % The size of kernel for an intial low pass average (0 indicates on averaging)
            {'seededWatershed'}, ...                                % A method in which a seed frame is used to create required features for a watershed approach
            'seededWatershed'};
        
        % Parameters defining the location of image data
        defaults(end+1,:) = {'watershedSeedChannel', ...      % The frame (or frames) used to define the seed associated with each watershed
            'string', 'DAPI'};
        defaults(end+1,:) = {'watershedChannel', ...          % The frame (or frames) used to define the watershed. 
            'string', 'polyT'};
        
        % Parameters for segmentation
        defaults(end+1,:) = {'seedFrameFilterSize', ...         % Size of a guassian kernal to filter frame
            'nonnegative', 5};
        defaults(end+1,:) = {'seedFrameErosionKernel', ...      % The kernal for erosion of the seed frame
            'freeType', [strel('sphere', 1) strel('disk', 20)]};
        defaults(end+1,:) = {'seedThreshold', ...               % The threshold for the seed frames
            'freeType', 'adaptive'};
        defaults(end+1,:) = {'seedConnectionKernel', ...        % The kernel for connecting nearby seeds
            'freeType', [strel('sphere', 1) strel('disk', 20)]};
        defaults(end+1,:) = {'seedDilationKernel', ...          % The kernel for dilating seed centroids prior to watershed
            'freeType', [strel('sphere', 1) strel('disk', 5)]};
        defaults(end+1,:) = {'minCellSize', ...                 % The minimum number of voxels in a cell
            'nonnegative', 100};
        defaults(end+1,:) = {'watershedFrameFilterSize', ...    % Size of guassian kernel to filter frame
            'nonnegative', 5};
        defaults(end+1,:) = {'watershedFrameThreshold', ...     % The brightness threshold for being in a cell
            'freeType', 'adaptive'};
        defaults(end+1,:) = {'ignoreZ', ...                     % Ignore z in the segmentation process
            'boolean', false};
        
        % Parameters for converting to real world coordinates
        defaults(end+1,:) = {'boundingBox', ...                 % The bounding box in microns centered on the middle of the fov
            'array', [-100 -100 200 200]};
        
        % Parameters defining stitching together features from different
        % FOV
        defaults(end+1,:) = {'maxEdgeDistance', ...             % The maximum distance between the end of an edge in one frame and that of another to be connected
            'nonnegative', 4};
        defaults(end+1,:) = {'maxFeatureCentroidDistance', ...
            'nonnegative', 5};
        
        % Parameters for parsing of barcodes into features
        defaults(end+1,:) = {'dilationSize', ...                % The fraction of a pixel size by which all boundaries will be expanded outwards to facilitate rapid parsing of barcodes
            'nonnegative', 0.1};
        
        % Parameters for display of results/archival of results
        defaults(end+1,:) = {'saveSegmentationReports', ...     % Should segmentation reports be generated and saved?
            'boolean', false};
                
    end
    
    % -------------------------------------------------------------------------
    % Default summation parameters
    % -------------------------------------------------------------------------
    function defaults = DefaultSummationParameters()
        % Generate the default parameters for the summation of raw data
   
        % Create empty cell
        defaults = cell(0,3);
        
        % Parameters for summation
        defaults(end+1,:) = {'areaBounds', ...                  % The lower and upper bounds on the area of individual features for calculation. The bounds are not inclusive.
            'nonnegative', [0 500]};
        defaults(end+1,:) = {'dcIndsForSummation', ...          % The indices for the data channels to use for summation
            'array', 17:40};
        defaults(end+1,:) = {'zIndForSummation', ...            % The indices associated with the z-channels for summation: REMOVE ME! NO LONGER USED!
            'array', []};
        
    end

    % -------------------------------------------------------------------------
    % Default smFISH molecule parameters
    % -------------------------------------------------------------------------
    function defaults = DefaultMoleculeParameters()
        % Generate the default parameters for the fast identification of
        % individual molecules
   
        % Create empty cell
        defaults = cell(0,3);
        
        % Parameters for summation
        defaults(end+1,:) = {'molLowPassfilterSize', ...   % The size of the low pass guassian filter in pixels
            'nonnegative', 5};
        defaults(end+1,:) = {'molIntensityThreshold', ...  % Intensity threshold
            'nonnegative', 1000};
        defaults(end+1,:) = {'molNumPixelSum', ...         % The number of pixels to sum in each direction for the brightness of the spot
            'nonnegative', 1};
        defaults(end+1,:) = {'molDataChannels', ...        % Information on the data channels to analyze
            'cell', {'RS0763', 4; 'RS1199', 4; 'RS1040', 4}};

    end

    
    % -------------------------------------------------------------------------
    % Generate Optimization Default Parameters
    % -------------------------------------------------------------------------
    function defaults = DefaultOptimizationParameters()
        % Generate the default parameters for the optimization process
        
        defaults = cell(0,3);
        defaults(end+1,:) = {'weightingOptimizationMethod', ... % Method for optimizing image weights
            {'equalOnBits'}, 'equalOnBits'};
        defaults(end+1,:) = {'quantileTarget', ...              % The quantile to set to 1 for initial weighting of histograms
            'nonnegative', 0.9};
        defaults(end+1,:) = {'areaThresh', ...                  % The area threshold for barcodes to be used in optimization
            'nonnegative', 4};
        defaults(end+1,:) = {'optNumFov', ...                   % The number of fov to use in the optimization process
            'positive', 50};
        defaults(end+1,:) = {'numIterOpt', 'nonnegative', 10};     % The number of iterations to perform in the optimization of weighting
        defaults(end+1,:) = {'normalizeToOne', 'boolean', false};   % Normalize the scale factors to 1?
    end    
    
    % -------------------------------------------------------------------------
    % Generate Default Display Parameters
    % -------------------------------------------------------------------------
    function defaults = DefaultDisplayParameters()
        % Generate the default parameters for display
        
        defaults = cell(0,3);
        defaults(end+1,:) = {'visibleOption', ...               % Display figures as they are created             
            {'on', 'off'}, 'on'};
        defaults(end+1,:) = {'overwrite', 'boolean', true};     % Overwrite existing figures
        defaults(end+1,:) = {'formats', ...                     % The figure formats to save
            'cell', {'fig', 'png'}};
        defaults(end+1,:) = {'useExportFig', 'boolean', false}; % Use the export_fig package (not always available)
        
        % Parameters for creating low resolution mosaics
        defaults(end+1,:) = {'downSample', 'nonnegative', 10};  % The value for downsampling images when creating low resolution mosaics
        defaults(end+1,:) = {'mosaicZInd', 'nonnegative', 4};  % The z index to select for generating mosaic indices
    end
    
    % -------------------------------------------------------------------------
    % Generate Default Quantification Parameters
    % -------------------------------------------------------------------------
    function defaults = DefaultQuantificationParameters()
        % Generate the default parameters for quantification of data
        
        defaults = cell(0,3);
        defaults(end+1,:) = {'minimumBarcodeArea', ...
            'nonnegative', 4};
        defaults(end+1,:) = {'minimumBarcodeBrightness', ...
            'nonnegative', 10^0.75};
        defaults(end+1,:) = {'minimumDistanceToFeature', ...
            'nonnegative', inf};
        defaults(end+1,:) = {'zSliceRange', ...
            'array', []};
        
    end
    
    % -------------------------------------------------------------------------
    % InitializeParameters
    % -------------------------------------------------------------------------
    function [parameters, defaultParameters] = InitializeParameters(varargin)
        % Create a default parameters structure and set fields in this
        % structure if specified
        %
        % parameters = obj.InitializeParameters(); % Return all defaults
        % parameters = obj.InitializeParameters('name', value, ...);
        % [~, defaultParameters] = obj.InitializeParameters(); % Return a
        %    cell array containing each parameters descriptor
        
        % Define default parameter sets
        defaultParameters = cell(0,2);
        defaultParameters(end+1,:) = {'warp', MERFISHDecoder.DefaultWarpParameters()};
        defaultParameters(end+1,:) = {'preprocess', MERFISHDecoder.DefaultPreprocessingParameters()};
        defaultParameters(end+1,:) = {'decoding', MERFISHDecoder.DefaultDecodingParameters()};
        defaultParameters(end+1,:) = {'optimization', MERFISHDecoder.DefaultOptimizationParameters()};
        defaultParameters(end+1,:) = {'display', MERFISHDecoder.DefaultDisplayParameters()};
        defaultParameters(end+1,:) = {'segmentation', MERFISHDecoder.DefaultSegmentationParameters()};
        defaultParameters(end+1,:) = {'quantification', MERFISHDecoder.DefaultQuantificationParameters()};
        defaultParameters(end+1,:) = {'summation', MERFISHDecoder.DefaultSummationParameters()};
        defaultParameters(end+1,:) = {'molecules', MERFISHDecoder.DefaultMoleculeParameters()};
        
        % Loop over parameter sets
        parameters = [];
        for p=1:size(defaultParameters,1)
            % Find matches to defaults
            localDefaults = defaultParameters{p,2};
            matchInd = find(ismember(varargin(1:2:end), localDefaults(:,1)));
            matchInd = sort([(2*matchInd -1) 2*matchInd]);
            subParameters = ParseVariableArguments(varargin(matchInd), localDefaults, mfilename);

            % Store parameters
            parameters.(defaultParameters{p,1}) = subParameters;
        end
    end
    
    % -------------------------------------------------------------------------
    % Combine and generate report for affine transforms
    % -------------------------------------------------------------------------
    function isCorrupt = CheckTiffStack(tiffFileName, expectedNumFrames)
        % Check the status of a tiff stack -- is it complete?

        isCorrupt = false; % Assume it is not corrupt
        try % See if one can open information on the tiff stack
            
            % Check for validity
            tiffInfo = imfinfo(tiffFileName);
            
            % Check to see if the number of frames are correct
            if length(tiffInfo) < expectedNumFrames
                isCorrupt = true; % If it can be opened, but there are not enough frames, it is corrupt
            end
            
        catch
            isCorrupt = true; % If the info file cannot be opened, it is corrupt
        end
    end
    
    % -------------------------------------------------------------------------
    % Decode an image
    % -------------------------------------------------------------------------
    function [decodedImage, localMagnitude, pixelTraces, D] = DecodePixels(imageStack, scaleFactors, decodingVectors, distanceThreshold)
        % ------------------------------------------------------------------------
        % [decodedImage, localMagnitude, imageStack] = obj.DecodePixels(imagePath, scaleFactors, decodingVectors, distanceThreshold);
        % This function takes an image stack (width, height, z, numFrames) and 
        % decodes it by comparing the individual normalized pixel vectors 
        % to the normalized vectors provided in the decodingVectors matrix. Pixels
        % outside of a given distance threshold to the nearest barcode are
        % discarded.
        %

        % Reshape imageStack to create pixel traces
        if ndims(imageStack) == 3 % Handle case of no z position
            imageHeight = size(imageStack,1);
            imageWidth = size(imageStack,2);
            stackLength = size(imageStack,3);
            numZPos = [];
            pixelTraces = reshape(imageStack, [imageHeight*imageWidth stackLength]);
        else
            imageHeight = size(imageStack,1);
            imageWidth = size(imageStack,2);
            numZPos = size(imageStack, 3);
            stackLength = size(imageStack,4);
            pixelTraces = reshape(imageStack, [imageHeight*imageWidth*numZPos stackLength]);
        end

        % Determine number of pixels
        numPixels = size(pixelTraces,1);

        % Type cast
        pixelTraces = double(pixelTraces);

        % Equalize bit brightness distributions via provided scale factors
        pixelTraces = pixelTraces./repmat(scaleFactors, [numPixels 1]);

        % Calculate the magnitude to normalize
        localMagnitude = sqrt(sum(pixelTraces.*pixelTraces, 2));

        % Create the pixel traces
        goodInds = localMagnitude > 0; % Handle the divide by zero case
        pixelTraces(goodInds,:) = pixelTraces(goodInds,:)./repmat(localMagnitude(goodInds), [1 size(pixelTraces,2)]);

        % Find nearest neighbor barcodes (on the N-dimensional unit sphere)
        [barcodeID, D] = knnsearch(decodingVectors, pixelTraces, 'K', 1);

        % Associate pixels to barcodes for which they are within a N-1 sphere
        % defined by the distanceThreshold
        exactInds = D <= distanceThreshold;

        % Decode image
        decodedImage = zeros([imageHeight imageWidth numZPos]);
        decodedImage(exactInds) = barcodeID(exactInds);
    end
    
    

end

end % end classdef
    
