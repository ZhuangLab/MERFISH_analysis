% Example MERFISH Analysis Script
% Jeffrey Moffitt
% January 30, 2016
% lmoffitt@mcb.harvard.edu
% -------------------------------------------------------------------------
% Purpose: To illustrate the analysis of MERFISH data. 
% -------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.

%% Setup Paths
analysisBasePath = '....'; % Insert path to folder for saving analysis
exampleDataPath = '....\MERFISH_Examples\'; % Insert path to example data 
% Example data can be downloaded from http://zhuang.harvard.edu/merfish/MERFISHData/MERFISH_Examples.zip

%% Setup parameters
% Setup parameters for parsing the name of image files
parameters.imageTag = 'STORM';          % Initial name--typically describes the scope
parameters.imageMListType = 'alist';    % Tag for molecule list file for found RNA spots
parameters.fiducialMListType = 'list';  % Tag for molecule list file for found fiducial spots

% Setup parameters of the run
parameters.numHybs = 16;                % The number of hybridization/imaging rounds
parameters.bitOrder = fliplr(1:16);     % The order in which bits are imaged. The example was run in reverse

% Setup parameters for constructing words
parameters.wordConstMethod = 'perLocalization'; % A flag to specify the word construction method. This is the preferred method.

% Setup parameters for decoding data
parameters.codebookPath = [exampleDataPath 'codebook\codebook.fasta'];           % Insert path to the codebook
parameters.exactMap = CodebookToMap(parameters.codebookPath, ...
        'keyType', 'binStr');
parameters.errCorrFunc = @SECDEDCorrectableWords;   % The function used to generate a map for correctable words
parameters.FPKMData = LoadByteStream(...
    [exampleDataPath '\FPKM_data\FPKMData.matb']);  % Insert path to the FPKMdata file

% Setup FOV/cells to analyze
parameters.cellsToAnalyze = [];         % Analyze all cells if empty

% Setup parameters for saving results
parameters.savePath = SetFigureSavePath(analysisBasePath, 'makeDir', true);

% Configure parameters for generating different reports
parameters.reportsToGenerate = cell(0,2);
parameters.reportsToGenerate(end+1,:) = {'fiducialReport2', 'off'}; % {report name, 'off'/'on' do not/do display figure}
parameters.reportsToGenerate(end+1,:) = {'numOnBitsHistByCell', 'off'}; 
parameters.reportsToGenerate(end+1,:) = {'focusLockReportByCell', 'off'};
parameters.reportsToGenerate(end+1,:) = {'totalFPKMReport', 'off'};
parameters.reportsToGenerate(end+1,:) = {'cellByCellFPKMReport', 'off'};
parameters.reportsToGenerate(end+1,:) = {'cellWithWordsImage', 'off'};
parameters.reportsToGenerate(end+1,:) = {'molStats', 'off'}; 
parameters.reportsToGenerate(end+1,:) = {'molDistStats', 'off'};
parameters.reportsToGenerate(end+1,:) = {'compositeHybImage', 'off'};
parameters.reportsToGenerate(end+1,:) = {'hamming1DReportAllGenes', 'off'};
parameters.reportsToGenerate(end+1,:) = {'bitFlipProbabilitiesAverage', 'off'};
parameters.reportsToGenerate(end+1,:) = {'bitFlipProbabilitiesAllGenes', 'off'};
parameters.reportsToGenerate(end+1,:) = {'confidenceRatioReport', 'off'};

parameters.overwrite = true;                % Overwrite existing files
parameters.figFormats = {'fig', 'png'};     % Output formats
parameters.useSubFolderForCellReport = true; 
parameters.saveAndClose = true; % Save figure once created, then close it

%% Run analysis!
[words, imageData, fiducialData, parameters] = AnalyzeMERFISH([exampleDataPath 'example_data'], ...
    'parameters', parameters);

%% Generate different summary reports
GenerateFPKMReport(words, parameters.FPKMData, 'parameters', parameters, ...
    'reportsToGenerate', {'totalFPKMReport', 'off'}, ...
    'FPKMReportExactMatchOnly', false, 'showNames', false); 
GenerateOnBitHistograms(words, 'parameters', parameters, 'reportsToGenerate',...
    {'numOnBitsHistAllCells', 'off'});
moleculeStats = GenerateMoleculeStatsReport(words, 'parameters', parameters);
bitFlipReport = GenerateBitFlipReport(words, parameters.exactMap, 'parameters', parameters);
GenerateHammingSphereReport(words, parameters.exactMap, 'parameters', parameters);
   
%% Save output
SaveAsByteStream([parameters.savePath 'imageData.matb'], imageData);
SaveAsByteStream([parameters.savePath 'fiducialData.matb'], fiducialData);
try
    save([parameters.savePath 'parameters.mat'], 'parameters');
catch
    warning('Corrupt parameters file.....'); % Often this structure is too large to save... odd matlab bug
end
SaveAsByteStream([parameters.savePath 'bitFlipReport.matb'], bitFlipReport);
SaveAsByteStream([parameters.savePath 'moleculeStats.matb'], moleculeStats);

% Save words: splitting due to issue with object size in SaveByteStream
wordInds = unique([1 300000:300000:length(words) length(words)]);
for j=2:length(wordInds)
    SaveAsByteStream([parameters.savePath 'words' num2str(j-1) '.matb'], ...
        words( (wordInds(j-1)+1):wordInds(j)));
end
    
% Save Stripped Words
strippedWords = StripWords(words); % Remove fields that are not typically used for most analysis
SaveAsByteStream([parameters.savePath 'strippedWords.matb'], strippedWords);
    
display('------------------------------------------------------------------');
display(['Saved words, imageData, parameters, and reports to ' parameters.savePath]);

% Save Script
try
    copyfile( [mfilename('fullpath'),'.m'],[parameters.savePath,mfilename,'.m']);
    display('------------------------------------------------------------------');
    display(['Copied analysis script to ' parameters.savePath,mfilename,'.m']);
catch
    warning('Could not find mfile path');
end    

%%
error('Script complete!');