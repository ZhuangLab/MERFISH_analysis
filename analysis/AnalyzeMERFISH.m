function [words, totalImageData, totalFiducialData, parameters] = AnalyzeMERFISH(dataPath, varargin)
% ------------------------------------------------------------------------
% [words, parameters] = AnalyzeMERFISH(dataPath, varargin)
% This function analyzes a series of raw conventional images in the 
%   specified directory and creates a words structure which represents all
%   of the identified words in the data.
%--------------------------------------------------------------------------
% Necessary Inputs
%--------------------------------------------------------------------------
% Outputs
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 3, 2014
%--------------------------------------------------------------------------
% Based on various scripts by Alistair Boettiger
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);

% Parameters for parsing file names
defaults(end+1,:) = {'imageTag', 'string', 'STORM'};        % Base tag for all images
defaults(end+1,:) = {'imageMListType', 'string', 'alist'};  % Flag for image mlist
defaults(end+1,:) = {'fiducialMListType', 'string','list'}; % Flag for fiducial mlists

% Parameters for parsing file names
defaults(end+1,:) = {'fileExt', 'string', 'bin'};           % Delimiters for bin files
defaults(end+1,:) = {'fieldNames', 'cell', ...              % Labels for fields in image name
    {'movieType', 'hybNum', 'cellNum', 'isFiducial', 'binType'}}; 
defaults(end+1,:) = {'fieldConv', 'cell', ...               % Conversion functions for fields in image name
    {@char, @str2num, @str2num, @(x)strcmp(x,'c2'), @char}};
defaults(end+1,:) = {'appendExtraFields', 'bool', true};    % How to handle names that don't match this pattern

% Parameters for the cells to analyze
defaults(end+1,:) = {'cellsToAnalyze', 'array', []};        % List of cell/FOV ids to analyze

% Parameters on the number of hybridizations
defaults(end+1,:) = {'numHybs', 'nonnegative', 16};         % Number of hybridizations
defaults(end+1,:) = {'bitOrder', 'boolean', 1:16};          % Order of bits

% Parameters for fiducial tracking
defaults(end+1,:) = {'maxD', 'nonnegative', 8};             % Maximum distance for fiducial tracking
defaults(end+1,:) = {'fiducialFrame', 'nonnegative', 1};    % Reference frame for fiducial markers
defaults(end+1,:) = {'fiducialWarp2Hyb1', 'boolean', false};

% Parameters for constructing words from spots
defaults(end+1,:) = {'maxDtoCentroid', 'nonnegative', 1};   % Distance between spots in different rounds

% Parameters for decoding words
defaults(end+1,:) = {'codebookPath', 'filePath', ''};       % Path to codebook
defaults(end+1,:) = {'codebook', 'struct', []};             % Codebook structure
defaults(end+1,:) = {'exactMap', 'map', []};                % containers.Map for decoding exact matches
defaults(end+1,:) = {'correctableMap', 'map', []};          % containers.Map for decoding correctable matches
defaults(end+1,:) = {'errCorrFunc', 'function', ...         % Error correction function
    @SECDEDCorrectableWords};
defaults(end+1,:) = {'keyType', {'int', 'binStr'}, ...      % Display type for binary word, e.g. binary or decimal
    'binStr'};

% Parameters for progress reports and intermediate figures
defaults(end+1,:) = {'savePath', 'fileDir', ''};            % Path to save incidental figures
defaults(end+1,:) = {'reportsToGenerate', 'cell', ...       % List of flags for generating different reports
    cell(0,2)}; 
defaults(end+1,:) = {'verbose', 'boolean', false};          % Display progress?
defaults(end+1,:) = {'printedUpdates', 'boolean', true};    % Display additional forms of progress?

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1 || ~(exist(dataPath) == 7) % 7=folder
    error('matlabFunctions:invalidArguments', 'A valid data path is required.');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Import Java utility for generating random IDs
% -------------------------------------------------------------------------
import java.util.UUID;
    
% -------------------------------------------------------------------------
% Provide overview of analysis
% -------------------------------------------------------------------------
parameterFieldsToDisplay = {'codebookPath', 'imageTag', 'imageMListType', ...
    'fiducialMListType', 'fileExt', 'numHybs', 'bitOrder', 'maxD', ...
    'wordConstMethod', 'savePath', ...
    };

if parameters.printedUpdates
    display('--------------------------------------------------------------');
    display('Analyzing Multiplexed FISH Data')
    display('--------------------------------------------------------------');
    display('Analysis parameters');
    for i=1:length(parameterFieldsToDisplay)
        display(['    ' parameterFieldsToDisplay{i} ': ' ...
            num2str(parameters.(parameterFieldsToDisplay{i}))]);
    end
end

% -------------------------------------------------------------------------
% Find data for analysis
% -------------------------------------------------------------------------
if parameters.printedUpdates
    display('--------------------------------------------------------------');
    display(['Finding cells in ' dataPath]);
end

foundFiles = BuildImageDataStructures(dataPath, 'parameters', parameters);
numCells = max([foundFiles.cellNum]);

if parameters.printedUpdates
    display(['    Found ' num2str(numCells) ' cells']);
end

if numCells == 0
    error('No valid cells found');
end

% -------------------------------------------------------------------------
% Load codebook and generate maps
% -------------------------------------------------------------------------
if isempty(parameters.codebook) && ~isempty(parameters.codebookPath)
    parameters.codebook = fastaread(parameters.codebookPath);
end

if ~isempty(parameters.codebook)
    parameters.exactMap = CodebookToMap(parameters.codebook, ...
        'keyType', parameters.keyType);
    if ~isempty(parameters.errCorrFunc)
        parameters.correctableMap = CodebookToMap(parameters.codebook, ...
            'keyType', parameters.keyType, ...
            'errCorrFunc', parameters.errCorrFunc, ...
            'mapContents', 'correctable');
    end
end

if parameters.printedUpdates
    display('--------------------------------------------------------------');
    if isempty(parameters.exactMap)
        display('No codebook provided. Found words will not be decoded.');
    elseif isempty(parameters.correctableMap)
        display('No error correction will be applied.');
    end
end

% -------------------------------------------------------------------------
% Prepare loop variables
% -------------------------------------------------------------------------
words = [];
totalImageData = [];
totalFiducialData = [];

% -------------------------------------------------------------------------
% Determine cells to analyze
% -------------------------------------------------------------------------
if isempty(parameters.cellsToAnalyze)
    cellIDs = 1:numCells;
else
    cellIDs = parameters.cellsToAnalyze(parameters.cellsToAnalyze >= 1 & ...
        parameters.cellsToAnalyze <= numCells);
end

% -------------------------------------------------------------------------
% Loop over all cells
% -------------------------------------------------------------------------
for i=cellIDs
    % ---------------------------------------------------------------------
    % Display cell number
    % ---------------------------------------------------------------------
    if parameters.printedUpdates
        display('--------------------------------------------------------------');
        display(['Analyzing data for cell ' num2str(i) ' of ' num2str(numCells)]);
    end
    
    % ---------------------------------------------------------------------
    % Identify all files for this cell
    % ---------------------------------------------------------------------
    imageData = foundFiles(strcmp({foundFiles.movieType}, parameters.imageTag) & ...
        strcmp({foundFiles.binType}, parameters.imageMListType) & ...
        [foundFiles.cellNum] == cellIDs(i) & ...
        ~[foundFiles.isFiducial]);
    
    fiducialData = foundFiles([strcmp({foundFiles.movieType}, parameters.imageTag)] & ...
        strcmp({foundFiles.binType}, parameters.fiducialMListType) & ...
        [foundFiles.cellNum] == cellIDs(i) & ...
        [foundFiles.isFiducial]);

    if parameters.printedUpdates & parameters.verbose
        display(['    Found ' num2str(length(imageData)) ' image files']);
        for j=1:length(imageData)
            display(['       ' imageData(j).filePath]);
        end
        display(['    Found ' num2str(length(fiducialData)) ' fiducial files']);
        for j=1:length(imageData)
            display(['       ' fiducialData(j).filePath]);
        end
    end
    
    % ---------------------------------------------------------------------
    % Cross checks on found files
    % ---------------------------------------------------------------------
    if length(imageData) ~= length(fiducialData)
        warning('matlabFunctions:AnalyzeMultiFISH', ...
            ['Cell ' num2str(i) ' does not have equal numbers of data and fiducial images']);
        continue;
    end
    if length(imageData) < parameters.numHybs
        warning('matlabFunctions:AnalyzeMultiFISH', ...
            ['Cell ' num2str(i) ' has fewer data images than hybs']);
        continue;
    end
    if length(fiducialData) < parameters.numHybs
        warning('matlabFunctions:AnalyzeMultiFISH', ...
            ['Cell ' num2str(i) ' has fewer fiducial images than hybs']);
        continue;
    end
    
    % ---------------------------------------------------------------------
    % Sort files based on hyb number: Almost certainly not necessary
    % ---------------------------------------------------------------------
    [~, sind] = sort([imageData.hybNum]);
    imageData = imageData(sind);
    [~, sind] = sort([fiducialData.hybNum]);
    fiducialData = fiducialData(sind);
    
    % ---------------------------------------------------------------------
    % Generate and append unique IDs for each image and fiducial data set
    % ---------------------------------------------------------------------
    for j=1:length(imageData)
        imageData(j).uID = char(UUID.randomUUID());
    end
    for j=1:length(fiducialData)
        fiducialData(j).uID = char(UUID.randomUUID());
    end
    
    % ---------------------------------------------------------------------
    % Load and transfer information on the corresponding dax
    % ---------------------------------------------------------------------
    imageData = TransferInfoFileFields(imageData, 'parameters', parameters);
%    fiducialData = TransferInfoFileFields(fiducialData, 'parameters', parameters);
    
    % ---------------------------------------------------------------------
    % Generate a measure of focus lock quality for all images
    % ---------------------------------------------------------------------
    %imageData = GenerateFocusLockQuality(imageData, 'parameters', parameters);

    % ---------------------------------------------------------------------
    % Load Molecule Lists and Fiducial Positions
    % ---------------------------------------------------------------------
    if parameters.printedUpdates & parameters.verbose
        display('--------------------------------------------------------------');
        display('Loading molecules lists');
    end
    for j=1:parameters.numHybs
        imageData(j).mList = ReadMasterMoleculeList(imageData(j).filePath, ...
            'compact', true, 'transpose', true, 'verbose', false);
        fiducialData(j).mList = ReadMasterMoleculeList(fiducialData(j).filePath, ...
            'compact', true, 'transpose', true, 'verbose', false);
        
        if parameters.printedUpdates & parameters.verbose
            display(['    ' imageData(j).name ': ' num2str(length(imageData(j).mList.x)) ' molecules']);
            display(['    ' fiducialData(j).name ': ' num2str(length(fiducialData(j).mList.x)) ' beads']);
        end
    end
     
    % ---------------------------------------------------------------------
    % Create Tiled Image (if desired)
    % ---------------------------------------------------------------------
    %GenerateTiledImage(imageData, 'parameters', parameters);
    
    % ---------------------------------------------------------------------
    % Add Transforms to Fiducial Data
    % ---------------------------------------------------------------------
    [fiducialData, parameters] = AlignFiducials(fiducialData, ...
        'parameters', parameters);
    
    % ---------------------------------------------------------------------
    % Transform Image Data and Transfer Fiducial Data
    % ---------------------------------------------------------------------
    [imageData, parameters] = TransformImageData(imageData,fiducialData,'parameters',parameters);
    
    % -------------------------------------------------------------------------
    % Create Words from Spots
    % -------------------------------------------------------------------------
    [wordsByCell, parameters] = CreateWords(imageData, 'parameters', parameters);

    % -------------------------------------------------------------------------
    % Decode words
    % -------------------------------------------------------------------------
    if ~isempty(parameters.codebook)
       [wordsByCell, parameters] = DecodeWords(wordsByCell, parameters.exactMap, ...
           parameters.correctableMap, 'parameters', parameters);
           if parameters.printedUpdates
                display(['    Found ' num2str(sum([wordsByCell.isExactMatch])) ' exact matches']);
                display(['    Found ' num2str(sum([wordsByCell.isCorrectedMatch])) ' corrected matches']);
           end
    end
    
    % ---------------------------------------------------------------------
    % Create Composite image with words
    % ---------------------------------------------------------------------
    GenerateCompositeImage(wordsByCell, imageData, 'parameters', parameters);
    
    % ---------------------------------------------------------------------
    % Create Cell By Cell On Bit Histogram
    % ---------------------------------------------------------------------
    GenerateOnBitHistograms(wordsByCell, 'parameters', parameters);
    
    % -------------------------------------------------------------------------
    % Append Words and imageData
    % -------------------------------------------------------------------------
    words = [words wordsByCell];
    totalImageData = [totalImageData imageData];
    totalFiducialData = [totalFiducialData fiducialData];
end

if parameters.printedUpdates
    display('--------------------------------------------------------------');
    display('Completed Multiplexed FISH Analysis');
end
