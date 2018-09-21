% Calculate Doublet Score
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 21, 2017
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2018.
% -------------------------------------------------------------------------
% Purpose: 1) Calculate and export various properties of found features
% that could be used to identify potential segmentation errors
% -------------------------------------------------------------------------

if ~exist('nDPath')
	error('A normalized data path must be provided.');
end
%% Create diary
if ~exist([nDPath 'log'], 'dir')
     mkdir([nDPath 'log']);
	 display(['Created ' nDPath 'log']); 
end
diary([nDPath 'log' filesep 'calculate.log']);
diary on;

%% Get slurm properties
PageBreak();
nodeID = getenv('SLURM_NODELIST');

display(['Running on ' nodeID]);
display(['Started at ' datestr(now)]);
scriptTimer = tic;

%% Load MERFISH Decoder
display(['Loading MERFISH Decoder from ' nDPath]);
mDecoder = MERFISHDecoder.Load(nDPath);

%% Load the found features
foundFeatures = mDecoder.GetFoundFeatures();

%% Extract the feature uids
feature_uID = {foundFeatures.uID};
feature_uID = feature_uID';
 

%% Load the barcode metadata
barcodeMetadataPath = [nDPath filesep 'reports' filesep 'barcode_metadata.csv'];
disp(['Loading barcode metadata']);
localTimer = tic;
barcodeTable = readtable(barcodeMetadataPath);
disp(['...completed in ' num2str(toc(localTimer)) ' s']);
disp(['...found ' num2str(size(barcodeTable,1)) ' barcodes']);

% Cut the barcodes to keep only those within features
barcodeTable = barcodeTable(barcodeTable.in_feature == 1, :);
disp(['...cutting to ' num2str(size(barcodeTable,1)) ' barcodes within features']);

%% Determine properties associated with these barcodes
numBarcodes = length(mDecoder.codebook);

%% Loop over features and compute numbers and mean positions
counts = zeros(length(feature_uID), numBarcodes);
comX = nan(length(feature_uID), numBarcodes);
varX= nan(length(feature_uID), numBarcodes);
comY = nan(length(feature_uID), numBarcodes);
varY = nan(length(feature_uID), numBarcodes);
for f=1:length(feature_uID)
    % Identify the index of this feature (unnecessary)
    local_feature_id = foundFeatures(strcmp(feature_uID, feature_uID{f})).feature_id;
    
    % Extract the local barcode table
    localBarcodeTable = barcodeTable(barcodeTable.feature_id == local_feature_id,:);
    
    % Compute the counts
    counts(local_feature_id,:) = accumarray(localBarcodeTable.barcode_id,1, [numBarcodes 1]);
    
    % Compute the center of mass for each barcode id
    comX(local_feature_id,:) = accumarray(localBarcodeTable.barcode_id, localBarcodeTable.abs_position_1, [numBarcodes 1], @mean, nan);
    comY(local_feature_id,:) = accumarray(localBarcodeTable.barcode_id, localBarcodeTable.abs_position_2, [numBarcodes 1], @mean, nan);
    
    % Compute variance
    varX(local_feature_id,:) = accumarray(localBarcodeTable.barcode_id, localBarcodeTable.abs_position_1, [numBarcodes 1], @var, nan);
    varY(local_feature_id,:) = accumarray(localBarcodeTable.barcode_id, localBarcodeTable.abs_position_2, [numBarcodes 1], @var, nan);

    % Display progress
    if ~mod(f, 1000)
        disp(['...completed ' num2str(f) ' of ' num2str(length(foundFeatures))]);
        disp(['...in ' num2str(toc(localTimer)) ' s']);
        localTimer = tic;
    end

end

%% Create the save path if needed
savePath = [nDPath filesep 'reports' filesep ];
if ~exist(savePath)
    mkdir(savePath);
end

%% Save these data
csvwrite([savePath 'barcode_counts.csv'], counts);
csvwrite([savePath 'barcode_center_of_mass_X.csv'], comX);
csvwrite([savePath 'barcode_center_of_mass_Y.csv'], comY);
csvwrite([savePath 'barcode_center_of_mass_var_X.csv'], varX);
csvwrite([savePath 'barcode_center_of_mass_var_Y.csv'], varY);

%% Archive analysis
PageBreak();
display(['...completed in ' num2str(toc(scriptTimer)) ' s']);
display(['Completed at ' datestr(now)]);

% Turn off diary
diary off;
