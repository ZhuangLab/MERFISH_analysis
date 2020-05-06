function [parameters] = WriteCodebook(codebookPath, barcodes, readouts, names, ids, varargin)
% ------------------------------------------------------------------------
% [codebook, header] = WriteCodebook(codebookPath, barcodes, readouts, geneNames, ids, varargin)
% This function writes a codebook file given the specified barcodes, names,
% and ids
%--------------------------------------------------------------------------
% Necessary Inputs
% --codebookPath/A path to the codebook that will be written
% --barcodes/A MxN binary matrix in specifying M barcodes of N bits in
% length
% -- names/A cell array of names for each barcode, e.g. gene names.
% -- ids/A cell array of ids for each barcode, e.g. transcript ids
% 
%--------------------------------------------------------------------------
% Outputs
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 21, 2017
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2018.
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true};           % Display progress?

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 5 
    error('matlabFunctions:invalidArguments', 'A valid path to a codebook is required');
end

if size(barcodes,1) < length(names)
    error('matlabFunctions:invalidArguments', 'Insufficient barcodes provided for the provided names');
end

if size(barcodes,2) < length(readouts)
    error('matlabFunctions:invalidArguments', 'Insufficient readouts for the number of bits in the barcodes');
end

if ~isstruct(readouts) || ~isfield(readouts, 'Header')
    error('matlabFunctions:invalidArguments', 'The provided readouts are not a structure with a Header field');
end

% Check file extension of provided file
[~, ~, fileExt] = fileparts(codebookPath);
if ~strcmp(fileExt, '.csv')
    warning('matlabFunctions:invalidExtension', 'It is recommended to save codebooks as csv files');
end

% -------------------------------------------------------------------------
% Display progress
% -------------------------------------------------------------------------
if parameters.verbose
    PageBreak();
    display(['Writing codebook: ' codebookPath]);
end

% -------------------------------------------------------------------------
% Open file
% -------------------------------------------------------------------------
fid = fopen(codebookPath, 'W');

% -------------------------------------------------------------------------
% Define and write header
% -------------------------------------------------------------------------
% Define bit names/labels
bitNames = {readouts(1:size(finalBarcodes,2)).Header};
bitLabels = cell(1, 2*length(bitNames) - 1);
bitLabels(1:2:end) = bitNames;
bitLabels(2:2:end) = repmat({', '}, [1 length(bitNames)-1]);
bitLabels = cat(2, bitLabels{:});

% Write flat header information
fprintf(fid, 'version, %s\n', '1.0');    
fprintf(fid, 'codebook_name, %s\n', 'M3E1'); 
fprintf(fid, 'bit_names, %s\n', bitLabels);
fprintf(fid, '%s, %s, %s\n', 'name', 'id', 'barcode');


    % Identify ids using target regions
    localIds = {}; % IDs
    for i=1:length(finalGenes)
        tRegions = finalTargetRegions(strcmp({finalTargetRegions.geneName}, finalGenes{i}));
        if ~isempty(tRegions)
            localIds{i} = tRegions.id;
            localIds{i} = strrep(localIds{i}, ',', '_');
        else
            localIds{i} = '';
        end
    end

    % Open codebook file
    

    % Write codebook data
    for i=1:length(finalGenes)
        fprintf(fid, '%s, %s, %s\n', finalGenes{i}, localIds{i}, num2str(finalBarcodes(i,:)));
    end
    
    % Write codebook
    PageBreak();
    display(['Wrote: ' codebookPath])
    
    % Close fid
    fclose(fid);
    
% else
%     error('Found existing codebook!');
end
