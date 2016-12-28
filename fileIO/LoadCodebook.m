function [codebook, header, parameters] = LoadCodebook(codebookPath, varargin)
% ------------------------------------------------------------------------
% [codebook, header] = LoadCodebook(codebookPath, varargin)
% This function creates a codebook structure from codebook file
%--------------------------------------------------------------------------
% Necessary Inputs
% --codebookPath/A path to a valid codebook
% 
%--------------------------------------------------------------------------
% Outputs
% --codebook/A structure array with the fields provided in the codebook
% file. These must include name and barcode fields. 
% --header/A structure with fields containing the initial header information
%  of the codebook. Required fields are 
%     --version/A string specifying the format version of the codebook
%     --bit_names/A cell array of the names of the bits in the order in which 
%        they occur in the barcodes.
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% March 13, 2016
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true};           % Display progress?
defaults(end+1,:) = {'barcodeConvFunc', 'function', @char}; % How to set the format of the barcode?
                                                            % e.g. @(x)x=='1'
% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1 || ~exist(codebookPath)
    error('matlabFunctions:invalidArguments', 'A valid path to a codebook is required');
end

% -------------------------------------------------------------------------
% Display progress
% -------------------------------------------------------------------------
if parameters.verbose
    PageBreak();
    display(['Loading codebook from: ' codebookPath]);
end

% -------------------------------------------------------------------------
% Open the file
% -------------------------------------------------------------------------
fid = fopen(codebookPath, 'r');
if fid < 0
    error('matlabFunctions:invalidArguments', 'A valid path to a codebook is required');
end

% -------------------------------------------------------------------------
% Load header
% -------------------------------------------------------------------------
done = false;
while ~done
    % Load line and split on comma
    line = fgetl(fid);
    stringParts = strsplit(line, ','); 
    
    % Remove whitespace
    stringParts = cellfun(@strtrim, stringParts, 'UniformOutput', false);
    
    % If the line is name, id, barcode, then the header portion of the
    % codebook has been read
    if length(stringParts) == 3 && isempty(setdiff(stringParts, {'name', 'id', 'barcode'}));
        done = true;
        break;
    end
    
    % Assign name value pairs
    if length(stringParts)==2
        header.(stringParts{1}) = stringParts{2};
    else
        header.(stringParts{1}) = stringParts(2:end);
    end
end

% -------------------------------------------------------------------------
% Check header
% -------------------------------------------------------------------------
if ~(isfield(header, 'version') & isfield(header, 'bit_names'))
    error('matlabFunctions:invalidArguments', 'The codebook is corrupt. Both a version and bit_names flag must be present.');
end

% -------------------------------------------------------------------------
% Display progress
% -------------------------------------------------------------------------
if parameters.verbose
    headerFields = fields(header);
    for f=1:length(headerFields)
        data = header.(headerFields{f});
        if iscell(data)
            displayString = cell(1, 2*length(data)-1);
            displayString(1:2:end) = data;
            displayString(2:2:end) = {', '};
            displayString = [displayString{:}];
        else
            displayString = data;
        end
        display(['...' headerFields{f} ': ' displayString]);
    end
end


% -------------------------------------------------------------------------
% Build codebook
% -------------------------------------------------------------------------
% Define function to remove whitespace
rmspace = @(x)x(~isspace(x));
% Switch based on version
switch header.version
    case {'1.0', '1'}
        codebook = repmat(struct('name', '', 'id', '', 'barcode', ''), [0 1]);
        
        done = false;
        while ~done
            % Check for end of file
            if feof(fid)
                done = true;
                break;
            end
            
            % Read line, split, and assign values
            line = fgetl(fid);
            stringParts = strsplit(line, ',');
            stringParts = cellfun(@strtrim, stringParts, 'UniformOutput', false);
            
            codebook(end+1).name = stringParts{1};
            codebook(end).id = stringParts{2};
            codebook(end).barcode = parameters.barcodeConvFunc(rmspace(stringParts{3}));
        end
end

% -------------------------------------------------------------------------
% Close the file
% -------------------------------------------------------------------------
fclose(fid);

% -------------------------------------------------------------------------
% Display progress
% -------------------------------------------------------------------------
if parameters.verbose
    display(['...loaded ' num2str(length(codebook)) ' barcodes']);
end


