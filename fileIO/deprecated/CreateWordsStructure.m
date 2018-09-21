function words = CreateWordsStructure(numElements, numHybs)
% ------------------------------------------------------------------------
% words = CreateWordsStructure(numElements, numHybs)
% This function creates an array of empty word structures.
%--------------------------------------------------------------------------
% Necessary Inputs
%   numElements/int. The number of elements to create.
%   numHybs/int. The number of hybs in the given words (used to preallocate
%   memory). 
%--------------------------------------------------------------------------
% Outputs
%   words/structure array. This array contains the following fields:
%       %% Codeword properties
%       --uID: A unique string for each word. Useful for archival and
%       indexing
%       --codeword: A logical array specifying the word in the bit order
%       corresponding to the codebook.
%       --measuredCodeword: A logical array specifing the value of each bit
%       in the order in which they were measured.
%       --intCodeword: The unsigned integer corresponding to the codeword
%       field
%       
%       %% Decoded word properties
%       --geneName: The name of the corresponding word in the codebook
%       --isExactMatch: A boolean which is determined if the codeword
%       corresponding to the geneName is an exact match to the measured
%       codeword.
%       --isCorrectedMatch: A boolean which represents if the geneName was
%       determined using some error correction
%
%       %% Covnenient shorthand properties of the codewords
%       --numOnBits: The number of on bits in the measured codeword
%       --onBits: The indices of the on bits in the order in the codebook
%       --measuredOnBits: The indicites of the on bits in the order in
%       which they were measured.
%
%       %% Image properties
%       --imageNames: A cell array of the names of the image files for each
%       of the measured images, in the measurement order
%       --imagePaths: Paths to these images
%       --imageUIDs: A cell array of unique IDs generated for each image
%       (see file structures)
%       --imageX: The x position of the image in um
%       --imageY: The y position of the image in um
%       
%       %% Properties of the codeword in the cell
%       --wordCentroidX: X position of the word in the warp reference frame
%       in pixels
%       --wordCentroidY: Y position
%       --cellID: The cell/FOV number from the movies/images used to
%       generate this word
%       --wordNumInCell: The number of the generated word in the FOV/cell,
%       e.g. 1st, 2nd, 3rd
%       --paddedCellID: An array of the cellID equal in length to the
%       number of on bits. Useful for indexing from word arrays.
%
%       %% Fiducial alignment/warp properties
%       --fiducialUIDs: The unique ID strings for each of the images used
%       to generate the fiducial warps
%       --hasFiducialError: A logical array specifying whether or not each
%       hyb had a fiducial error.  In the order in which the hybs were
%       measured.
%       
%       %% Focus lock properties
%       --focusLockQuality: A number determining the quality of the focus
%       for each hyb image.  
%
%       %% Molecule list properties
%       Words also contain all of the fields found in molecule lists. See
%       CreateMoleculeList and the analysis software, e.g. DaoSTORM, for
%       more details. 

%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
% -- None.
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 25, 2014
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 2 || numElements < 0 || numHybs < 0
    error('matlabFunctions:invalidArguments', 'Invalid values for numElements and numHybs');
end

% -------------------------------------------------------------------------
% Define fields and defaultProperties
% -------------------------------------------------------------------------
fieldsAndValues = cell(0,2);

% Codeword values and properties
fieldsAndValues(end+1,:) = {'uID', ''};
fieldsAndValues(end+1,:) = {'codeword', false(1, numHybs)};
fieldsAndValues(end+1,:) = {'measuredCodeword', false(1, numHybs)}; % A logical array specifying the measured word in the or
fieldsAndValues(end+1,:) = {'intCodeword', nan};

fieldsAndValues(end+1,:) = {'geneName', ''};
fieldsAndValues(end+1,:) = {'isExactMatch', false};
fieldsAndValues(end+1,:) = {'isCorrectedMatch', false};
fieldsAndValues(end+1,:) = {'numOnBits', 0};
fieldsAndValues(end+1,:) = {'onBits', nan(1, numHybs)};
fieldsAndValues(end+1,:) = {'measuredOnBits', nan(1, numHybs)};

% Measurement properties
fieldsAndValues(end+1,:) = {'bitOrder', 1:numHybs};
fieldsAndValues(end+1,:) = {'numHyb', nan};
fieldsAndValues(end+1,:) = {'imageNames', cell(1, numHybs)};
fieldsAndValues(end+1,:) = {'imagePaths', cell(1, numHybs)};
fieldsAndValues(end+1,:) = {'imageUIDs', cell(1, numHybs)};
fieldsAndValues(end+1,:) = {'imageX', nan};
fieldsAndValues(end+1,:) = {'imageY', nan};

% Properties of the codeword in the cell
fieldsAndValues(end+1,:) = {'wordCentroidX', nan};
fieldsAndValues(end+1,:) = {'wordCentroidY', nan};
fieldsAndValues(end+1,:) = {'cellID', nan};
fieldsAndValues(end+1,:) = {'wordNumInCell', nan};

% Fiducial alignment/warp properties
fieldsAndValues(end+1,:) = {'fiducialUIDs', cell(1, numHybs)};
fieldsAndValues(end+1,:) = {'hasFiducialError', false(1, numHybs)};
fieldsAndValues(end+1,:) = {'paddedCellID', nan(1, numHybs)};

% Focus Lock Properties
fieldsAndValues(end+1,:) = {'focusLockQuality', nan(1, numHybs)};

% Molecule properties
fieldsAndValues(end+1,:) = {'mListInds', nan(1, numHybs)};

% -------------------------------------------------------------------------
% Generate fields and types from molecules lists
% -------------------------------------------------------------------------
defaultMList = CreateMoleculeList(numHybs);
mListFields = fields(defaultMList);

% -------------------------------------------------------------------------
% Create Words Structure
% -------------------------------------------------------------------------
for i=1:length(fieldsAndValues)
    defaultWord.(fieldsAndValues{i,1}) = fieldsAndValues{i,2};
end

% -------------------------------------------------------------------------
% Transfer Molecule list fields
% -------------------------------------------------------------------------
for i=1:length(mListFields)
    defaultWord.(mListFields{i}) = defaultMList.(mListFields{i})'; % Transpose is a kludge to produce row vectors from the default column vector
end

% -------------------------------------------------------------------------
% Create Array
% -------------------------------------------------------------------------
words = repmat(defaultWord, [1 numElements]);

