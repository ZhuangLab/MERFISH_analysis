function parameters = MERFISHPerformanceMetrics(normalizedDataPath, varargin)
% ------------------------------------------------------------------------
% MERFISHPerformanceMetrics(normalizedDataPath, varargin)
% This function generates a series of important quantifications of the
% MERFISH data as well as several useful performance metrics. 

%--------------------------------------------------------------------------
% Necessary Inputs: 
%   normalizedDataPath -- A valid path to a normalized MERFISH data path
%--------------------------------------------------------------------------
% Outputs: 
%   --None
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
% 
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% June 5, 2016
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2018.


% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);

% Parameters for displaying progress
defaults(end+1,:) = {'verbose', 'boolean', false};      % Display progress?
defaults(end+1,:) = {'logProgress', 'boolean', true};   % Create log file?
defaults(end+1,:) = {'archive', 'boolean', true};       % Create copies of the utilized functions

% Parameters for location of barcodes
defaults(end+1,:) = {'barcodePath', 'fileDir', []};         % The path to the barcodes

% Parameters for location of saved reports
defaults(end+1,:) = {'outputPath', 'fileDir', []};          % Path to save generated reports and metrics

% Parameters for analyzing barcodes
defaults(end+1,:) = {'codebookPath', ...                % Path to the codebook
    'filePath', []};
defaults(end+1,:) = {'abundDataPath', ...               % Path to FPKM data
    'filePath', []};
defaults(end+1,:) = {'cellBoundariesPath', ...          % Path to a file containing cell boundaries
    'filePath', []};

% Parameters for parallel processing
defaults(end+1,:) = {'parallel', 'parallel', []};       % A parallel.pool object

% Parameters for handling cells vs fov
defaults(end+1,:) = {'cellIDMethod', ...                % The method for assigning RNAs to cells
    {'fov', 'cellID'}, 'fov'};
defaults(end+1,:) = {'blockSize', ...                   % The number of barcodes to load at a time
    'nonnegative', 1e5};

% Parameters for selecting barcodes
defaults(end+1,:) = {'brightnessThreshold', ...         % The minimum brightness to save a barcode
    'nonnegative', 10^0.75};
defaults(end+1,:) = {'areaThreshold', ...               % The minimum area to save a barcode
    'nonnegative', 2};
defaults(end+1,:) = {'stageOrientation', ...            % The orientation of the different stage axis wrt to the camera
    'array', [1 -1]};

% Parameters for generating reports
defaults(end+1,:) = {'visibleOption', ...               % Display figures or just generate, save, and close
    {'on', 'off'}, 'on'};
defaults(end+1,:) = {'brightnessBins', ...              % The histogram bins for log10 brightness
    'array', 0:.025:4};
defaults(end+1,:) = {'areaBins', ...                    % The histogram bins for area
    'array', 1:1:15};

% Parameters for finding blanks
defaults(end+1,:) = {'blankFnc', 'function', ...        % Function to identify blank controls in codebook
    @(x)~isempty(regexp(x, 'Blank-', 'once'))};

% Parameters for barcode density report                 % The size of 2D histogram bins in pixels
defaults(end+1,:) = {'barcodeDensityBinSize', ...
    'nonnegative', 20};

% Parameters for SaveFigure
defaults(end+1,:) = {'overwrite', 'boolean', true};     % Options for SaveFigure
defaults(end+1,:) = {'formats', 'cell', {'png', 'fig'}}; 
defaults(end+1,:) = {'useExportFig', 'boolean', false};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1 || ~exist(normalizedDataPath)
    error('matlabFunctions:invalidArguments', ...
        'A valid normalized data path must be provided');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Define paths to various items
% -------------------------------------------------------------------------
% Define default path and (alternate default) for codebook
if isempty(parameters.codebookPath)
    % Look for the a codebook file
    foundFiles = BuildFileStructure(normalizedDataPath, ...
        'regExp', '(?<codebookName>\w+)_codebook', ...
        'fileExt', 'csv', ...
        'fieldNames', {'codebookName'}, 'fieldConv', {@char});
    if isempty(foundFiles)
        error('matlabFunctions:missingItem', 'Could not find a valid codebook');
    elseif length(foundFiles) > 2
        error('matlabFunctions:missingItem', 'Found two many files that look like codebooks');
    end
    PageBreak();
    disp(['Utilizing a found codebook for ' foundFiles(1).codebookName]);
    disp(['... ' foundFiles(1).filePath]);
    parameters.codebookPath = foundFiles(1).filePath;
end

% Define path to barcodes
if isempty(parameters.barcodePath)
    parameters.barcodePath = [normalizedDataPath filesep 'barcodes' filesep];
end

% Define default output path
if isempty(parameters.outputPath)
    parameters.outputPath = [parameters.barcodePath ...
        filesep 'performance' filesep];
end

% Create this path if it does not exist
if ~exist(parameters.outputPath)
    mkdir(parameters.outputPath);
end

% -------------------------------------------------------------------------
% Start log (if requested)
% -------------------------------------------------------------------------
if parameters.logProgress
    logFilePath = [parameters.outputPath 'performance.log'];
    if exist(logFilePath) % Delete existing file
        diary off;
        delete(logFilePath);
    end
    diary(logFilePath); % Set diary file
    diary on;

    % Display information
    PageBreak();
    disp(['Started log: ' datestr(now)]);
end

% ------------------------------------------------------------------------
% Display progress
% -------------------------------------------------------------------------
PageBreak();
disp(['Calculating performance metrics for ' normalizedDataPath]);
disp(['Started: ' datestr(now)]);

% -------------------------------------------------------------------------
% Determine properties of parallel.pool if provided
% -------------------------------------------------------------------------
if ~isempty(parameters.parallel)
    numPar = parameters.parallel.NumWorkers;
    disp(['Utilizing ' num2str(numPar) ' parallel workers']);
else
    numPar = 0;
end

% ------------------------------------------------------------------------
% Load codebook and create decoding maps
%-------------------------------------------------------------------------
% Load codebook
[codebook, codebookHeader] = LoadCodebook(parameters.codebookPath, 'verbose', true, ...
    'barcodeConvFunc', @(x)bi2de(x=='1'));
bits = codebookHeader.bit_names;

% Archive numbers
numBarcodes = length(codebook);
numBits = length(bits);

% Prepare variables for abundance correlation correlation
if exist(parameters.abundDataPath) 
    PageBreak();
    disp(['Loading abundance data for correlation from: ' parameters.abundDataPath]);
    abundData = LoadByteStream(parameters.abundDataPath);
    geneNames = {codebook.name};

    commonNames = intersect(geneNames, {abundData.geneName});

    disp(['...found ' num2str(length(commonNames)) ' values that overlap with codebook entries']);
else
    PageBreak();
    disp(['No abundance data provided.']);
end

% ------------------------------------------------------------------------
% Load the mDecoder
%-------------------------------------------------------------------------
mDecoder = MERFISHDecoder.Load(normalizedDataPath);

% ------------------------------------------------------------------------
% Load and compile data either via different methods for handling different
% 'cells'
%-------------------------------------------------------------------------
% Archive progress
PageBreak();
disp(['Compiling properties by cell using the ' parameters.cellIDMethod ' method']);
totalTimer = tic;

% Define the path based on the method, define how these will be indexed as
% well
switch parameters.cellIDMethod
    case 'cellID' % Handle the case that cells were parsed and barcodes contain cellIDs
        % Define path to parsed barcodes
        bListPath = [parameters.barcodePath filesep 'parsed' ...
            filesep 'assigned_blist.bin'];
                
        % Read flat header of assigned barcode list to determine number of entries
        flatHeader = ReadBinaryFileHeader(bListPath);

        % Break indices into blocks of desired size
        inds = [0:parameters.blockSize:(flatHeader.numEntries-1) flatHeader.numEntries];
        numObjects = length(inds)-1;
        
        % Define cell boundaries path
        if isempty(parameters.cellBoundariesPath)
            parameters.cellBoundariesPath = [normalizedDataPath ...
                filesep 'segmentation' filesep 'rawCellBoundaries.matb'];
        end
        
        % Load cell boundaries
        cellBoundaries = LoadByteStream(parameters.cellBoundariesPath, ...
            'verbose', parameters.verbose);
        
        % Define the number of cells
        numCells = length(cellBoundaries);
        
        % Define the y-labels for plot below
        FPKMCorrYLabel = 'Counts/Cell';
        
    case 'fov'
        % Define path to barcodes for each fov
        bListPath = [parameters.barcodePath filesep 'barcode_fov' filesep];
        
        % Build the file structure for these barcodes
        foundFiles = BuildFileStructure(bListPath, ...
            'fileExt', 'bin', ...
            'regExp', 'fov_(?<fov>[0-9]+)', ...
            'fieldNames', {'fov'}, ...
            'fieldConv', {@str2num});
        
        % Display progress
        disp(['Found ' num2str(length(foundFiles)) ' barcode files']);
        
        % Define properties for iteration
        numObjects = length(foundFiles);
        numCells = length(foundFiles);
        
        % Define the y-labels for plot below
        FPKMCorrYLabel = 'Counts/FOV';

end    


% Load and process blocks of barcodes or barcodes by fov in parallel
spmd (numPar)

    % Define local variables for accumulation per worker
    countsPerCellExact = zeros(numBarcodes, numCells);
    countsPerCellCorrected = zeros(numBarcodes, numCells);

    numZero2One = zeros(numBarcodes, numBits+1); % 0 is also included to signify no error
    numOne2Zero = zeros(numBarcodes, numBits+1);

    brightnessAreaHist = zeros(length(parameters.brightnessBins), ...
        length(parameters.areaBins));

    distToNucleusOutNucleus = zeros(numBarcodes, numCells); % Average distance per RNA to nucleus if the RNA is outside of the nucleus
    distToNucleusInNucleus = zeros(numBarcodes, numCells); % Average distance per RNA to nucleus if the RNA is inside the nucleus
    fractionInNucleus = zeros(numBarcodes, numCells); % Fraction of RNAs inside the nucleus

    barcodeDensity = zeros(length(1:parameters.barcodeDensityBinSize:mDecoder.imageSize(1)), ... % Histogram of barcode density per fov
        length(1:parameters.barcodeDensityBinSize:mDecoder.imageSize(2)), mDecoder.numZPos);  % The histogram will be calculated for all z positions individually
    
    % Loop over a set of blocks per worker
    for b=labindex:numlabs:numObjects
        % Display progress
        displayStrings = {};

        % Load barcodes 
        loadTimer = tic;
        switch parameters.cellIDMethod
            case 'cellID'
                displayStrings{end+1} = (['Loading block ' num2str(b) ' of ' num2str(length(inds)-1)]);
                aList = ReadBinaryFile(bListPath, 'first', inds(b)+1, ...
                            'last', inds(b+1));
            case 'fov'
                displayStrings{end+1} = (['Loading fov ' num2str(foundFiles(b).fov) ' of ' num2str(numObjects)]);
                aList = ReadBinaryFile(foundFiles(b).filePath);
        end

        % Check for empty barcodes
        if isempty(aList)
            continue;
        end

        displayStrings{end+1} = (['Loaded ' num2str(length(aList)) ' barcodes in ' num2str(toc(loadTimer)) ' s']); 

        % Compute area/brightness histogram prior to cuts
        histogramTimer = tic;

        brightnessAreaHist = brightnessAreaHist + hist3(cat(2, log10([aList.total_magnitude]./double([aList.area]))', ...
            double([aList.area])'), 'Edges', {parameters.brightnessBins, parameters.areaBins});

        displayStrings{end+1} = (['Computed brightness and area histogram in ' num2str(toc(histogramTimer)) ' s']); 

        % Cut barcodes
        cutTimer = tic;
        aList = aList([aList.area] >= parameters.areaThreshold & ...
            [aList.total_magnitude]./double([aList.area]) >= parameters.brightnessThreshold);

        displayStrings{end+1} = (['Cut to ' num2str(length(aList)) ' barcodes in ' num2str(toc(cutTimer)) ' s']); 

        % Check for empty barcodes
        if isempty(aList)
            continue;
        end

        % Extract cellID (or fov_id)
        switch parameters.cellIDMethod
            case 'cellID'
                cellID = [aList.cellID];
            case 'fov'
                cellID = [aList.fov_id];
        end
        % Compute histograms for counts per cell: exact and corrected
        histogramTimer = tic;
        data = cat(2,[aList([aList.is_exact]==1).barcode_id]', ...
            cellID([aList.is_exact]==1)');
        if ~isempty(data)
            countsPerCellExact = countsPerCellExact + hist3(data, 'Edges', {1:numBarcodes, 1:numCells});
        end
        data = cat(2,[aList([aList.is_exact]==0).barcode_id]', ...
            cellID([aList.is_exact]==0)');
        if ~isempty(data)
            countsPerCellCorrected = countsPerCellCorrected + hist3(data, 'Edges', {1:numBarcodes, 1:numCells});
        end
        displayStrings{end+1} = (['Computed counts per cell histograms in ' num2str(toc(histogramTimer)) ' s']); 

        % Accumulate errors at each bit for each barcode
        histogramTimer = tic;
        data = cat(2, [aList([aList.error_dir]==0).barcode_id]', ...
            [aList([aList.error_dir]==0).error_bit]');
        if ~isempty(data)
            numZero2One = numZero2One + hist3(data, 'Edges', {1:numBarcodes, 0:numBits});
        end

        data = cat(2, [aList([aList.error_dir]==1).barcode_id]', ...
            [aList([aList.error_dir]==1).error_bit]');
        if ~isempty(data)
            numOne2Zero = numOne2Zero + hist3(data, 'Edges', {1:numBarcodes, 0:numBits});
        end

        displayStrings{end+1} = (['Computed errors in ' num2str(toc(histogramTimer)) ' s']); 

        % Accumulate the barcode density as a function of location
        densityTimer = tic;
    
        barcodeCenter = cat(1, aList.weighted_pixel_centroid); % Extract the centers of all barcodes
        % Build a list of z indices for all barcodes
        if size(barcodeCenter,2) == 3
            zInds = discretize(barcodeCenter(:,3), [1:mDecoder.numZPos inf]);
        else
            zInds = ones(1, size(barcodeCenter,1));
        end
        for z=1:mDecoder.numZPos
            barcodeDensity(:,:,z) = barcodeDensity(:,:,z) + hist3(barcodeCenter(zInds == z,1:2), ...
                'Edges', {1:parameters.barcodeDensityBinSize:mDecoder.imageSize(1), ...
                1:parameters.barcodeDensityBinSize:mDecoder.imageSize(2)});
        end
        
        displayStrings{end+1} = (['Computed barcode density in ' num2str(toc(densityTimer)) ' s']); 
        
        % Display complete progress
        disp(char(displayStrings));

    end
end
% Compile composite objects
PageBreak();
disp(['Combining composite objects']);
timer = tic;

temp = countsPerCellExact{1};
for i=2:length(countsPerCellExact)
    temp = temp + countsPerCellExact{i};
end
countsPerCellExact = temp;

temp = countsPerCellCorrected{1};
for i=2:length(countsPerCellCorrected)
    temp = temp + countsPerCellCorrected{i};
end
countsPerCellCorrected = temp;

temp = numZero2One{1};
for i=2:length(numZero2One)
    temp = temp + numZero2One{i};
end
numZero2One = temp;

temp = numOne2Zero{1};
for i=2:length(numOne2Zero)
    temp = temp + numOne2Zero{i};
end
numOne2Zero = temp;

temp = brightnessAreaHist{1};
for i=2:length(brightnessAreaHist)
    temp = temp + brightnessAreaHist{i};
end
brightnessAreaHist = temp;

temp = distToNucleusOutNucleus{1};
for i=2:length(distToNucleusOutNucleus)
    temp = temp + distToNucleusOutNucleus{i};
end
distToNucleusOutNucleus = temp;

temp = distToNucleusInNucleus{1};
for i=2:length(distToNucleusInNucleus)
    temp = temp + distToNucleusInNucleus{i};
end
distToNucleusInNucleus = temp;

temp = fractionInNucleus{1};
for i=2:length(fractionInNucleus)
    temp = temp + fractionInNucleus{i};
end
fractionInNucleus = temp;

temp = barcodeDensity{1};
for i=2:length(barcodeDensity)
    temp = temp + barcodeDensity{i};
end
barcodeDensity = temp;
disp(['...completed in ' num2str(toc(timer)) ' s']);

disp(['Completed compiling performance statistics at ' datestr(now)]);
disp(['...in ' num2str(toc(totalTimer)) ' s']);

% Normalize the distances
% Add all counts
totalCounts = countsPerCellExact + countsPerCellCorrected;

% Compute number out nucleus
numOutNucleus = totalCounts - fractionInNucleus;

% Normalize distances
distToNucleusInNucleus = distToNucleusInNucleus./fractionInNucleus; % fractionInNucleus contains the number of counts in nucleus at this point
distToNucleusOutNucleus = distToNucleusOutNucleus./numOutNucleus;

% Normalize fraction
fractionInNucleus = fractionInNucleus./totalCounts;

% Handle division by zero
distToNucleusInNucleus(isnan(distToNucleusInNucleus)) = 0;
distToNucleusOutNucleus(isnan(distToNucleusOutNucleus)) = 0;
fractionInNucleus(isnan(fractionInNucleus)) = 0;

% ------------------------------------------------------------------------
% Save the calculated data
%-------------------------------------------------------------------------
if ~exist(parameters.outputPath)
    mkdir(parameters.outputPath);
end

PageBreak();
disp(['Writing data']);
csvwrite([parameters.outputPath 'countsPerCellExact.csv'], countsPerCellExact);
disp(['...wrote ' parameters.outputPath 'countsPerCellExact.csv']);

csvwrite([parameters.outputPath 'countsPerCellCorrected.csv'], countsPerCellCorrected);
disp(['...wrote ' parameters.outputPath 'countsPerCellCorrected.csv']);

csvwrite([parameters.outputPath 'numZero2One.csv'], numZero2One);
disp(['...wrote ' parameters.outputPath 'numZero2One.csv']);

csvwrite([parameters.outputPath 'numOne2Zero.csv'], numOne2Zero);
disp(['...wrote ' parameters.outputPath 'numOne2Zero.csv']);

csvwrite([parameters.outputPath 'brightnessAreaHist.csv'], brightnessAreaHist);
disp(['...wrote ' parameters.outputPath 'brightnessAreaHist.csv']);

for z=1:size(barcodeDensity,3) % Write all of the z position densities separately
    csvwrite([parameters.outputPath 'barcodeDensity-' num2str(z) '.csv'], barcodeDensity);
    disp(['...wrote ' parameters.outputPath 'barcodeDensity-' num2str(z) '.csv']);
end

if strcmp(parameters.cellIDMethod, 'cellID')
    csvwrite([parameters.outputPath 'distToNucleusInNucleus.csv'], distToNucleusInNucleus);
    disp(['...wrote ' parameters.outputPath 'distToNucleusInNucleus.csv']);

    csvwrite([parameters.outputPath 'distToNucleusOutNucleus.csv'], distToNucleusOutNucleus);
    disp(['...wrote ' parameters.outputPath 'distToNucleusOutNucleus.csv']);

    csvwrite([parameters.outputPath 'fractionInNucleus.csv'], fractionInNucleus);
    disp(['...wrote ' parameters.outputPath 'fractionInNucleus.csv']);
end

% -------------------------------------------------------------------------
% Create area-brightness report
%--------------------------------------------------------------------------
% Create figure handle
figHandle = figure('Name', 'Area and brightness histograms', 'Color', 'w', ...
    'visible', parameters.visibleOption, 'Position', [1 1 1618 420]);

% Area distribution
subplot(1,3,1);
bar(parameters.areaBins, sum(brightnessAreaHist,1), 1, 'EdgeColor', 'none'); hold on;
xlabel('Area (pixels)');
ylabel('Counts');
plot([1 1]*parameters.areaThreshold, ylim, '--', 'Color', [0.75 0.75 0.75]);
xlim([parameters.areaBins(1)-1 parameters.areaBins(end)+1]);

% Brightness distribution (log10)
subplot(1,3,2);
bar(parameters.brightnessBins, sum(brightnessAreaHist,2), 1, 'EdgeColor', 'none'); hold on;
xlabel('Brightness (log_{10})');
ylabel('Counts');
plot([1 1]*log10(parameters.brightnessThreshold), ylim, '--', 'Color', [0.75 0.75 0.75]);
xlim([parameters.brightnessBins(1) parameters.brightnessBins(end)] + ...
    [1 -1]*mean(diff(parameters.brightnessBins)));

% Area/brightness distributions
subplot(1,3,3);
for i=1:length(parameters.areaBins)
    p = brightnessAreaHist(:,i)/max(brightnessAreaHist(:,i)); % Normalize to 1
    plot([-p' fliplr(p') -p(1)]/2.5 + parameters.areaBins(i), ...
        [parameters.brightnessBins fliplr(parameters.brightnessBins) parameters.brightnessBins(1)], ...
        'b'); hold on;
end
% Add thresholds
xlabel('Area (pixels)');
ylabel('Brightness (log_{10})');
xlim([parameters.areaBins(1)-1 parameters.areaBins(end)+1]);
plot(xlim,  log10(parameters.brightnessThreshold)*[1 1],'--', 'Color', [0.75 0.75 0.75]);
plot((parameters.areaThreshold-0.5)*[1 1], ylim, '--',  'Color', [0.75 0.75 0.75]);

% Save figure
SaveFigure(figHandle, 'parameters', parameters, 'savePath', parameters.outputPath);
close(figHandle);


% -------------------------------------------------------------------------
% Create barcode density report
%--------------------------------------------------------------------------
% Create separate reports for all z planes
for z=1:size(barcodeDensity,3)
    figHandle = figure('Name', ['Barcode density z ' num2str(z)], 'Color', 'w', ...
        'visible', parameters.visibleOption, 'Position', [1 1 1000 1000]);

    % Plot the density
    imagesc('XData', 1:parameters.barcodeDensityBinSize:mDecoder.imageSize(2), ...
        'YData', 1:parameters.barcodeDensityBinSize:mDecoder.imageSize(1), ...
        'CData', barcodeDensity(:,:,z))
    xlim([1 mDecoder.imageSize(2)]);
    ylim([1 mDecoder.imageSize(1)]);
    xlabel('X Position (pixels)');
    ylabel('Y Position (pixels)');
    cb = colorbar;
    cb.Label.String = 'Number';

    % Save figure
    SaveFigure(figHandle, 'parameters', parameters, 'savePath', parameters.outputPath);
    close(figHandle);
end

% Create a combined report 
figHandle = figure('Name', ['Barcode density'], 'Color', 'w', ...
    'visible', parameters.visibleOption, 'Position', [1 1 1000 1000]);

% Plot the density
imagesc('XData', 1:parameters.barcodeDensityBinSize:mDecoder.imageSize(2), ...
    'YData', 1:parameters.barcodeDensityBinSize:mDecoder.imageSize(1), ...
    'CData', sum(barcodeDensity,3))
xlim([1 mDecoder.imageSize(2)]);
ylim([1 mDecoder.imageSize(1)]);
xlabel('X Position (pixels)');
ylabel('Y Position (pixels)');
cb = colorbar;
cb.Label.String = 'Number';

% Save figure
SaveFigure(figHandle, 'parameters', parameters, 'savePath', parameters.outputPath);
close(figHandle);

% -------------------------------------------------------------------------
% Create FPKM correlation plot
%--------------------------------------------------------------------------
if ~isempty(parameters.abundDataPath) % Check for FPKM data
    % Create figure handle
    figHandle = figure('Name', 'FPKM correlation', 'Color', 'w', ...
        'visible', parameters.visibleOption, 'Position', [1 1 1000 667]);

    % Sort codebook and abundance data to match names
    [~, sortIndBarcodes, sindB] = intersect(geneNames, {abundData.geneName});
    sortedFPKM = [abundData(sindB).FPKM];
    
    % Exact FPKM Correlation
    subplot(2,3,1);
    nA = mean(countsPerCellExact,2);
    nA1 = nA(sortIndBarcodes); % Sort to FPKM data
    goodInds = nA1' > 0 & sortedFPKM > 0;
    if sum(goodInds) > 2
        loglog(sortedFPKM(goodInds), nA1(goodInds), '.'); hold on;
        [R,P] = corrcoef(log10(sortedFPKM(goodInds)), log10(nA1(goodInds)));
        title(['Exact: \rho_{10}: ' num2str(R(1,2),2) '; P: ' num2str(P(1,2),2)]);
        ylim([min(nA)*0.8 max(nA)*1.2]);
        xlim([min(sortedFPKM)*0.8 max(sortedFPKM)*1.2]);
    end
    xlabel('FPKM');
    ylabel(FPKMCorrYLabel);
    axis square;

    % Enrichment/deenrichment
    subplot(2,3,4);
    if sum(goodInds) > 2
        ratio = log2(nA1(goodInds)'./sortedFPKM(goodInds));
        [~, sortInd] = sort(nA1(goodInds), 'descend');
        plot(ratio(sortInd) - mean(ratio), '.');
        title([num2str(std(ratio), 2)]);
    end
    xlabel('Barcode');
    ylabel('Ratio (log_{2})');
    axis square;

    % Corrected FPKM Correlation
    subplot(2,3,2);
    nA = mean(countsPerCellCorrected,2);
    nA1 = nA(sortIndBarcodes); % Sort to FPKM data
    goodInds = nA1' > 0 & sortedFPKM > 0;
    if sum(goodInds) > 2
        loglog(sortedFPKM(goodInds), nA1(goodInds), '.'); hold on;
        [R,P] = corrcoef(log10(sortedFPKM(goodInds)), log10(nA1(goodInds)));
        title(['Corrected: \rho_{10}: ' num2str(R(1,2),2) '; P: ' num2str(P(1,2),2)]);
        ylim([min(nA)*0.8 max(nA)*1.2]);
        xlim([min(sortedFPKM)*0.8 max(sortedFPKM)*1.2]);
    end
    xlabel('FPKM');
    ylabel(FPKMCorrYLabel);
    axis square;

    % Enrichment/deenrichment
    subplot(2,3,5);
    if sum(goodInds)>2
        ratio = log2(nA1(goodInds)'./sortedFPKM(goodInds));
        [~, sortInd] = sort(nA1(goodInds), 'descend');
        plot(ratio(sortInd) - mean(ratio), '.');
        title([num2str(std(ratio), 2)]);
    end
    xlabel('Barcode');
    ylabel('Ratio (log_{2})');
    axis square;

    % Total FPKM Correlation
    subplot(2,3,3);
    nA = mean(countsPerCellCorrected+countsPerCellExact,2);
    nA1 = nA(sortIndBarcodes); % Sort to FPKM data
    goodInds = nA1' > 0 & sortedFPKM > 0;
    if sum(goodInds) > 2
        loglog(sortedFPKM(goodInds), nA1(goodInds), '.'); hold on;
        [R,P] = corrcoef(log10(sortedFPKM(goodInds)), log10(nA1(goodInds)));
        title(['All: \rho_{10}: ' num2str(R(1,2),2) '; P: ' num2str(P(1,2),2)]);
        ylim([min(nA)*0.8 max(nA)*1.2]);
        xlim([min(sortedFPKM)*0.8 max(sortedFPKM)*1.2]);
    end
    xlabel('FPKM');
    ylabel(FPKMCorrYLabel);
    axis square;

    % Enrichment/deenrichment
    subplot(2,3,6);
    if sum(goodInds) > 2
        ratio = log2(nA1(goodInds)'./sortedFPKM(goodInds));
        [~, sortInd] = sort(nA1(goodInds), 'descend');
        plot(ratio(sortInd) - mean(ratio), '.');
        title([num2str(std(ratio), 2)]);
    end
    xlabel('Barcode');
    ylabel('Ratio (log_{2})');
    axis square;

    %  Save figure
    SaveFigure(figHandle, 'parameters', parameters, 'savePath', parameters.outputPath);
    close(figHandle);
end

% -------------------------------------------------------------------------
% Per-bit error report
%--------------------------------------------------------------------------
% Create figure handle
figHandle = figure('Name', 'Error report', 'Color', 'w', 'visible', parameters.visibleOption, ...
    'Position', [1 1 1135 420]);

% Plot barcode numbers
subplot(2,4,1);
bar([sum(countsPerCellCorrected(:)) sum(countsPerCellExact(:))]);
set(gca, 'XTick', 1:2, 'XTickLabel', {'Corrected', 'Exact'}, 'XTickLabelRotation', 45);
ylabel('Counts');
title([num2str(sum(countsPerCellCorrected(:))/(sum(countsPerCellCorrected(:)+countsPerCellExact(:)))*100,2) '%']);

% Plot num above blanks for different data sets
dataSets = {countsPerCellCorrected, countsPerCellExact, countsPerCellCorrected + countsPerCellExact};
labels = {'Corrected', 'Exact', 'All'};

% Determine blank/non-blank inds
blankInds = find(cellfun(parameters.blankFnc, {codebook.name}));
nonBlankInds = find(~cellfun(parameters.blankFnc, {codebook.name}));

for i=1:length(dataSets) % Loop over the different data combinations
    subplot(2,4,i+1);
    nA = sum(dataSets{i},2);

    [sortednA, sind] = sort(nA, 'descend');
    x = 1:length(codebook);
    x = x(sind);
    localBlankInds = find(ismember(x, blankInds));
    bar(1:length(codebook), sortednA, 1, 'b', 'EdgeColor', 'none'); hold on;
    bar(localBlankInds, sortednA(localBlankInds), 1, 'r', 'EdgeColor', 'none');
    set(gca,'YScale', 'log');
    xlabel('Barcode ID');
    ylabel(['Counts ' labels{i}]);
    xlim([0 length(codebook)+1]);
    
    maxBlank = max(nA(blankInds));
    title(['Number above: ' num2str(sum(nA(nonBlankInds) > maxBlank))]);
end

% Plot confidence ratio
confidenceRatio = sum(countsPerCellExact,2)./(sum(countsPerCellCorrected,2)+sum(countsPerCellExact,2));

subplot(2,4,5);
[sortedCR, sind] = sort(confidenceRatio, 'descend');
x = 1:length(codebook);
x = x(sind);
localBlankInds = find(ismember(x, blankInds));
bar(1:length(codebook), sortedCR, 1, 'b', 'EdgeColor', 'none'); hold on;
bar(localBlankInds, sortedCR(localBlankInds), 1, 'r', 'EdgeColor', 'none');
set(gca,'YScale', 'log');
xlabel('Barcode ID');
ylabel(['Confidence ratio']);
xlim([0 length(codebook)+1]);

maxBlank = max(confidenceRatio(blankInds));
title(['Number above: ' num2str(sum(confidenceRatio(nonBlankInds) > maxBlank))]);

% Plot per bit error rates
normalizedOne2Zero = numOne2Zero(:,2:end)./repmat(sum(countsPerCellExact+countsPerCellCorrected,2),[1 numBits]);
normalizedZero2One = numZero2One(:,2:end)./repmat(sum(countsPerCellExact+countsPerCellCorrected,2),[1 numBits]);

% 1->0
subplot(2,4,6);
bar(nanmean(normalizedOne2Zero,1), 'b');
ylabel('Error rate');
xlabel('Bit');
title(['1 \rightarrow 0: ' num2str(nanmean(nanmean(normalizedOne2Zero,1)),2)]);
xlim([0 length(bits)+1]);

% 0->1
subplot(2,4,7);
bar(nanmean(normalizedZero2One,1), 'b');
ylabel('Error rate');
xlabel('Bit');
title(['0 \rightarrow 1: ' num2str(nanmean(nanmean(normalizedZero2One,1)),2)]);
xlim([0 length(bits)+1]);

SaveFigure(figHandle, 'parameters', parameters, 'savePath', parameters.outputPath);
close(figHandle);

% ------------------------------------------------------------------------
% Archival
% -------------------------------------------------------------------------
% Record completion
PageBreak();
disp(['Completed decoding of ' normalizedDataPath]);
disp(['...at ' datestr(now)]);

% Archive code
if parameters.archive
    PageBreak();
    functionPath = [mfilename('fullpath') '.m'];
    disp(['Archiving ' functionPath ' to ' parameters.outputPath]);
    copyfile(functionPath, [parameters.outputPath mfilename '.m']);
    SaveAsByteStream([parameters.outputPath 'performance_parameters.matb'], parameters, 'verbose', parameters.verbose);
end

% Stop logging
if parameters.logProgress
   diary off;
end

