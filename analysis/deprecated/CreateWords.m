function [words, parameters] = CreateWords(imageData,varargin)
% ------------------------------------------------------------------------
% words = CreateWords(imageData,varargin)
% This generates a series of word structures from an imageData structure by
%   first finding connected objects across hyb images and then using the
%   centroids of these objects to construct actual words. 
%--------------------------------------------------------------------------
% Necessary Inputs
% imageData/A structure array with elements equal to the number of
%   images to align. For information on the fields of this structure see
%   CreateImageData
%--------------------------------------------------------------------------
% Outputs
% words/A structure array of words with the following fields
%   --
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%    'binSize' -- size of bins in nm
%    'minDotPerBin' -- min number of localizations to call a bin occupied
%    'minLocsPerDot' -- min number of localization in all bins assigned to a cluster to be called an mRNA
%    'minArea' -- min area in bins to be called a cluster of localization
%    'maxArea' -- max area in bins to be called a cluster of localization
%    'maxDtoCentroid' -- the maximum distance from a centroid to any
%       individual hybridization 
%--------------------------------------------------------------------------
% Alistair Boettiger
% boettiger@fas.harvard.edu
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 9, 2014
%--------------------------------------------------------------------------
% Based on FindMRNA.m and AssignConvSpotsToCentroids.m
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'wordConstMethod', {'commonCentroid','perLocalization'}, 'perLocalization'};
defaults(end+1,:) = {'binSize', 'positive', 0.25};
defaults(end+1,:) = {'minDotPerBin', 'nonnegative',1};
defaults(end+1,:) = {'minLocsPerDot', 'positive',1};
defaults(end+1,:) = {'minArea', 'nonnegative',0};
defaults(end+1,:) = {'maxArea', 'nonnegative',10};
defaults(end+1,:) = {'showPlots', 'boolean', false};
defaults(end+1,:) = {'verbose', 'boolean', false};
defaults(end+1,:) = {'clusterFig', 'handle',[]};
defaults(end+1,:) = {'histFig', 'handle',[]};
defaults(end+1,:) = {'imageSize', 'array', [256,256]};
defaults(end+1,:) = {'verbose', 'boolean', true};
defaults(end+1,:) = {'printedUpdates', 'boolean', true};
defaults(end+1,:) = {'numHybs', 'positive', 16};

% Report parameters
defaults(end+1,:) = {'reportsToGenerate', 'cell', []}; 
defaults(end+1,:) = {'useSubFolderForCellReport', 'boolean', true};
defaults(end+1,:) = {'overwrite', 'boolean', true};
defaults(end+1,:) = {'figFormats', 'cell', {'png', 'fig'}};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 2
    error('matlabSTORM:invalidArguments', 'two nx2 vectors of points are required');
end

% -------------------------------------------------------------------------
% Import Java utility for generating random IDs
% -------------------------------------------------------------------------
import java.util.UUID;

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Pull out mLists
% -------------------------------------------------------------------------
mLists = [imageData.mList];
mListFields = fields(mLists(1));

% -------------------------------------------------------------------------
% Clear words
% -------------------------------------------------------------------------
words = [];

% -------------------------------------------------------------------------
% Display status
% -------------------------------------------------------------------------
if parameters.printedUpdates & parameters.verbose
    display('--------------------------------------------------------------');
    display(['Creating words with method: ' parameters.wordConstMethod]);
end

% -------------------------------------------------------------------------
% Select word construction method
% -------------------------------------------------------------------------
switch parameters.wordConstMethod
    case 'commonCentroid'
        % -------------------------------------------------------------------------
        % Identify putative words by finding connected regions in all images
        % -------------------------------------------------------------------------
        % Create image parameters
        xall = [mLists.xc]';
        yall = [mLists.yc]';

        M = hist3([yall,xall],{parameters.binSize:parameters.binSize:parameters.imageSize(1),...
            parameters.binSize:parameters.binSize:parameters.imageSize(2)});

        % -------------------------------------------------------------------------
        % Find connected regions in combined images
        % -------------------------------------------------------------------------
        P = regionprops(M>=parameters.minDotPerBin,M,'PixelValues','PixelIdxList','Centroid','Area');

        % -------------------------------------------------------------------------
        % Cut regions
        % -------------------------------------------------------------------------
        clusterLocs = cellfun(@sum,{P.PixelValues});
        clusterarea = [P.Area];
        goodClusters = clusterLocs>parameters.minLocsPerDot & clusterarea > parameters.minArea & clusterarea < parameters.maxArea;

        % -------------------------------------------------------------------------
        % Identify Potential Centroids
        % -------------------------------------------------------------------------
        putativeWordCentroids = cat(1,P.Centroid)*parameters.binSize;
        putativeWordCentroids = putativeWordCentroids(goodClusters,:);

        % -------------------------------------------------------------------------
        % Calculate distances
        % -------------------------------------------------------------------------
        idx = zeros(length(putativeWordCentroids), length(imageData));
        d = inf(length(putativeWordCentroids), length(imageData)); % If empty, no distances are triggered
        for i=1:length(imageData)
            if ~isempty(imageData(i).mList.xc) % Handle lost frames
                [idx(:,i), d(:,i)] = knnsearch([imageData(i).mList.xc',imageData(i).mList.yc'],putativeWordCentroids);
            end
        end
        
        % -------------------------------------------------------------------------
        % Allocate Word Structure
        % -------------------------------------------------------------------------
        words = CreateWordsStructure(length(putativeWordCentroids), parameters.numHybs);
        
        % -------------------------------------------------------------------------
        % Build word properties specific to method
        % -------------------------------------------------------------------------
        for i=1:length(putativeWordCentroids)
            words(i).measuredCodeword = d(i,:) <= parameters.maxDtoCentroid;
            words(i).mListInds = idx(i, words(i).measuredCodeword);
            words(i).wordCentroidX = putativeWordCentroids(i,1);
            words(i).wordCentroidY = putativeWordCentroids(i,2);
        end
        
    case 'perLocalization'
        minPhotsPerStain = 1; 
        numHybes = length(imageData); 
        % record the positions of all spots in all hybes;  
        spotPostionsPerHybe = cell(numHybes,1); 
        for h=1:numHybes
            spotPostionsPerHybe{h} = [imageData(h).mList.xc',imageData(h).mList.yc'];
        end
        putativeWordCentroids = cat(1,spotPostionsPerHybe{:}); 
        
        if all(cellfun(@(x)~isempty(x), spotPostionsPerHybe)) % Don't build words for cells without all mLists
            % Assign localizations to mRNA centroids
            wordsDetected = false(length(putativeWordCentroids),numHybes);
            brightnessPerSpot = NaN(length(putativeWordCentroids),numHybes);
            xPerSpot = NaN(length(putativeWordCentroids),numHybes);
            yPerSpot = NaN(length(putativeWordCentroids),numHybes);
            idxPerSpot= NaN(length(putativeWordCentroids),numHybes);
            for h=1:numHybes
                [idx,di] = knnsearch([imageData(h).mList.xc',imageData(h).mList.yc'],putativeWordCentroids);
                brightnessPerSpot(:,h) = imageData(h).mList.a(idx);
                brightnessPerSpot(di >= parameters.maxDtoCentroid,h) = 0; 
                validIdx = brightnessPerSpot(:,h) > minPhotsPerStain;
                wordsDetected(:,h) = validIdx; 
                xPerSpot(:,h) = imageData(h).mList.xc(idx); 
                xPerSpot(~validIdx,h) = NaN; 
                yPerSpot(:,h) = imageData(h).mList.yc(idx);
                yPerSpot(~validIdx,h) = NaN; 
                idxPerSpot(:,h) = idx; 
                idxPerSpot(~validIdx,h) = NaN; 
            end

            % compute centroids of all codewords
            meanX = nanmean(xPerSpot,2);
            meanY = nanmean(yPerSpot,2);
            wordLocations = [meanX,meanY];

            % Remove redundantly recorded codewords ID'd by shared centroids.
            [wordLocations,uniqueIdx] = unique(wordLocations,'rows','stable');      
            wordsDetected = wordsDetected(uniqueIdx,:);
            idxPerSpot = idxPerSpot(uniqueIdx,:); 


            % -------------------------------------------------------------------------
            % Allocate Word Structure
            % -------------------------------------------------------------------------
            numWords = size(wordsDetected,1);
            words = CreateWordsStructure(numWords, parameters.numHybs);

            % -------------------------------------------------------------------------
            % Build word properties specific to method
            % -------------------------------------------------------------------------
            for i=1:numWords
                words(i).measuredCodeword = wordsDetected(i,:);
                words(i).mListInds = idxPerSpot(i,~isnan( idxPerSpot(i,:) ));
                words(i).wordCentroidX = wordLocations(i,1);
                words(i).wordCentroidY = wordLocations(i,2);
            end
        else
            words = CreateWordsStructure(0, parameters.numHybs);
        end
    otherwise
        error('matlabFunctions:CreateWords', 'Unknown word construction method')
end

% -------------------------------------------------------------------------
% Fill out word structure
% -------------------------------------------------------------------------
for i=1:length(words)
    % Transfer basic names and experiment info
    words(i).imageNames = {imageData.name};
    words(i).imagePaths = {imageData.filePath};
    words(i).bitOrder = parameters.bitOrder;
    words(i).numHyb = parameters.numHybs;
    words(i).cellID = imageData(1).cellNum;
    words(i).wordNumInCell =i;
    
    % Transfer image position
    words(i).imageX = imageData(1).Stage_X;
    words(i).imageY = imageData(1).Stage_Y;
    
    % Generate and transfer unique IDs
    words(i).uID = char(UUID.randomUUID());
    words(i).imageUIDs = {imageData.uID};
    words(i).fiducialUIDs = {imageData.fidUID};

    % Record properties of identification
    words(i).hasFiducialError = [imageData.hasFiducialError];
    
    % Determine measured and actual codeword
    words(i).codeword = words(i).measuredCodeword(words(i).bitOrder);
    words(i).intCodeword = bi2de(words(i).codeword,'left-msb'); %Save an integer version of the codeword
    words(i).numOnBits = sum(words(i).measuredCodeword);
    
    % Save On Bit Indicies
    words(i).measuredOnBits = find(words(i).measuredCodeword);
    words(i).onBits = find(words(i).codeword);
    words(i).paddedCellID = ones(1, words(i).numOnBits,'int32')*words(i).cellID;
    
    % Transfer mList properties
    for k = 1:length(words(i).measuredOnBits)
        listID = words(i).measuredOnBits(k);
        moleculeID = words(i).mListInds(k);
        for fieldID = 1:length(mListFields)
            words(i).(mListFields{fieldID})(listID) = mLists(listID).(mListFields{fieldID})(moleculeID);
        end
    end
    
    % Add fields for identifying words/genes
    words(i).geneName = '';
    words(i).isExactMatch = false;
    words(i).isCorrectedMatch = false;
    
    % Add focus lock quality data
    words(i).focusLockQuality = [imageData.focusLockQuality];
end

% -------------------------------------------------------------------------
% Display progress
% -------------------------------------------------------------------------
if parameters.printedUpdates
    display(['    Reconstructed ' num2str(length(words)) ' words']);
    if parameters.verbose
        [n, x] = hist([words.numOnBits], 1:parameters.numHybs);
        for j=1:parameters.numHybs
            display(['        ' num2str(n(j)) ' words with ' num2str(x(j)) ' on bits']);
        end
    end
end

