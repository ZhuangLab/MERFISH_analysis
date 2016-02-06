function [map, geneNames, codewords, parameters] = CodebookToMap(codebook, varargin)
% ------------------------------------------------------------------------
% [map, geneNames, codeWords] = CodebookToMap(codebook, varargin)
% This function generates a container.Map object from a either the path to
%   a valid codebook or a codebook structure. 
%--------------------------------------------------------------------------
% Necessary Inputs
% --codebook/Either a path to a valid codebook or a codebook structure. 
% 
%--------------------------------------------------------------------------
% Outputs
% --map/A container.Map object with keys corresponding to the valid 
%   codewords described in the codebook (and error correctable codewords) 
%   and the names of the corresponding objects (genes). 
% --geneNames/A cell array of the value entries for the map
% --codewords/A cell array of the key entries for the map
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 7, 2014
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'errCorrFunc', 'function', []};
defaults(end+1,:) = {'keyType', {'int', 'binStr'}, 'int'};
defaults(end+1,:) = {'mapContents', {'exact', 'all', 'correctable'}, 'all'};

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1 
    error('matlabFunctions:invalidArguments', 'A codebook is required');
end

% -------------------------------------------------------------------------
% Parse provided data
% -------------------------------------------------------------------------
if ischar(codebook)
    if ~exist(codebook, 'file')
        error('matlabFunctions:invalidArguments', 'The provided path is not to a valid file.');
    else
        try
            codebook = fastaread(codebook);
        catch
            error('matlabFunctions:invalidArguments', 'The provided path is not to a valid if !fasta file.');
        end
    end
end

if isstruct(codebook)
    if ~all(ismember({'Header', 'Sequence'}, fields(codebook)))
        error('matlabFunctions:invalidArguments', 'The provided structure does not have the correct fields.');
    end
end

% -------------------------------------------------------------------------
% Determine Key Conversion Function
% -------------------------------------------------------------------------
removeWs = @(x) x(~isspace(x)); % Useful short hand for removing whitespace
switch parameters.keyType
    case 'int'
        % keyConv= @(x) uint16(bin2dec(removeWs(num2str(x))));
        keyConv= @(x) uint16( bi2de(logical(str2num(x)),'left-msb') ); 
    case 'binStr'
        keyConv = removeWs;
end

% -------------------------------------------------------------------------
% Determine values for map object
% -------------------------------------------------------------------------
keys = {};
values = {};
for i=1:length(codebook)
    keys{i} = keyConv(codebook(i).Header);
    wsInds = find(isspace(codebook(i).Sequence));
    geneName = codebook(i).Sequence(1:(wsInds(1)-1));
    %geneName = strsplit(codebook(i).Sequence);
    values{i} = geneName;
end

% -------------------------------------------------------------------------
% Return geneNames and codewords in the order in the codebook
% -------------------------------------------------------------------------
geneNames = values;
codewords = keys;

% -------------------------------------------------------------------------
% Create map
% -------------------------------------------------------------------------
exactMap = containers.Map(keys, values);

% -------------------------------------------------------------------------
% Add correctable keys: Pass as logical arrays
% -------------------------------------------------------------------------
keys = {}; values = {};
if ~isempty(parameters.errCorrFunc)
    newKeyConv = @(x) keyConv(num2str(x));
    for i=1:length(codebook)
        codewordLogical = removeWs(codebook(i).Header) == '1';
        newKeys = cellfun(newKeyConv, parameters.errCorrFunc(codewordLogical), ...
            'UniformOutput', false);
        if ~isempty(newKeys)
            wsInds = find(isspace(codebook(i).Sequence));
            geneName = codebook(i).Sequence(1:(wsInds(1)-1));
            newValues = repmat({geneName}, [1 length(newKeys)]);
            keys((end+1):(end+length(newKeys))) = newKeys;
            values((end+1):(end+length(newKeys))) = newValues;
        end
    end
    correctableMap = containers.Map(keys, values);
else
    map = exactMap;
end

% -------------------------------------------------------------------------
% Return map
% -------------------------------------------------------------------------
if ~isempty(parameters.errCorrFunc)
    switch parameters.mapContents
        case 'exact'
            map = exactMap;
        case 'correctable'
            map = correctableMap;
        case 'all'
            map = [exactMap; correctableMap];
    end
end


