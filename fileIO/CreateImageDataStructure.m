function imageData = CreateImageDataStructure(numElements)
% ------------------------------------------------------------------------
% imageData = CreateImageDataStructure(numElements)
% This function creates an array of empty imageData structures.
%--------------------------------------------------------------------------
% Necessary Inputs
%   numElements/int. The number of elements to create.
%--------------------------------------------------------------------------
% Outputs
% imageData/A structure array of length specified by numElements. Each
%   element contains the following fields
%
%       %% File information
%       --name: The name of the file
%       --filePath: The path to the file
%       --infFilePath: The path to the corresponding .inf file
%       --uID: A string unique to this instance of this element
%
%       %% Movie/Image information
%       --movieType: A string specifying the type of movie, e.g. STORM,
%       bleach
%       --hybNum: The number of the hyb
%       --cellNum: The number of the FOV or cell
%       --isFiducial: A boolean specifying whether the image is of fiducial
%       markers
%       --binType: A string defining the type of bin file to use for
%       analysis, e.g. alist or med300_alist
%       --delimiters: A cell array of the delimiters used to parse the
%       filename
%
%       %% Movie/Image information from .inf file
%       --imageH: The height of the image in pixels
%       --imageW: The width of the image in pixels
%       --Stage_X: The x position of the stage in um
%       --Stage_Y: The y position of the stage in um
%
%       %% Focus lock
%       --focusLockQuality: A scaler specifying the quality of the focus
%       lock for this image
%
%       %% Molecules
%       --mList: A structure containing information on all of the molecules
%       identified in this image.  See CreateMoleculeList for field
%       information.
%
%       %% Warp/Alignment
%       --tform: An affine transformation structure used to align this
%       image to a common coordinate system
%       --warpErrors: A 1x5 vector containing errors associated with
%       aligning fiducial markers
%       --hasFiducialError: A boolean specifying whether there was an error
%       in warping this image
%       --fiducialErrorMessage: The provided error message
%       --fiducialUID: The unique string (ID) of the imageData structure
%       corresponding to the fiducial image used to warp this image.
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
if nargin < 1 || numElements < 0
    error('matlabFunctions:invalidArguments', 'Invalid values for numElements and numHybs');
end

% -------------------------------------------------------------------------
% Define fields and defaultProperties
% -------------------------------------------------------------------------
fieldsAndValues = cell(0,2);

% File name and unique ID
fieldsAndValues(end+1,:) = {'name', ''};
fieldsAndValues(end+1,:) = {'filePath', ''};
fieldsAndValues(end+1,:) = {'infFilePath', ''};
fieldsAndValues(end+1,:) = {'uID', ''};

% Movie type and details
fieldsAndValues(end+1,:) = {'movieType', ''};
fieldsAndValues(end+1,:) = {'hybNum', -1};
fieldsAndValues(end+1,:) = {'cellNum', -1};
fieldsAndValues(end+1,:) = {'isFiducial', false};
fieldsAndValues(end+1,:) = {'binType', ''};
fieldsAndValues(end+1,:) = {'delimiters', {}};

% Movie properties from inf file
fieldsAndValues(end+1,:) = {'imageH', 0};
fieldsAndValues(end+1,:) = {'imageW', 0};
fieldsAndValues(end+1,:) = {'Stage_X', 0};
fieldsAndValues(end+1,:) = {'Stage_Y', 0};

% Movie focus quality
fieldsAndValues(end+1,:) = {'focusLockQuality', 0};

% Molecule lists
fieldsAndValues(end+1,:) = {'mList', CreateMoleculeList(0)};

% Warp properties
fieldsAndValues(end+1,:) = {'tform', maketform('affine',[1 0 0; 0 1 0; 0 0 1])};
fieldsAndValues(end+1,:) = {'warpErrors', zeros(1,5)};
fieldsAndValues(end+1,:) = {'hasFiducialError', false};
fieldsAndValues(end+1,:) = {'fiducialErrorMessage', []};
fieldsAndValues(end+1,:) = {'fidUID', ''};

% -------------------------------------------------------------------------
% Create imageData Structure
% -------------------------------------------------------------------------
for i=1:length(fieldsAndValues)
    defaultImageData.(fieldsAndValues{i,1}) = fieldsAndValues{i,2};
end

imageData = repmat(defaultImageData, [1 numElements]);

