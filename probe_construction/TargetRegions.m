classdef TargetRegions < handle
% ------------------------------------------------------------------------
% [targetRegions] = TargetRegions(varargin)
% This class stores target region information. 
%
%--------------------------------------------------------------------------
% Necessary Inputs
%--------------------------------------------------------------------------
% Methods
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% May 3, 2015
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Define properties
% -------------------------------------------------------------------------
properties (SetAccess=protected)
    geneName = ''           % Common name of the gene
    id = ''                 % Accession for the gene
    sequence = {}           % Sequences of the target regions
    startPos = []           % Starting position for each target region
    regionLength = []       % Length for each target region
    GC = []                 % GC for each target region
    Tm = []                 % Melting temperature for each target region
    specificity = []        % Specificty of each target region wrt to all other genes (not including isoforms)
    isoSpecificity = []     % Specificity of each target region wrt to other isoforms
    penalties = []          % Penalty values for each target region (by row)
    penaltyNames = {}       % Name of each penalty (i.e. each row in penalities)
    numRegions = 0          % Number of target regions in class
end

properties (SetAccess=protected, Hidden=true)
    map = '*ACGTRYKMSWBDHVN-'; 
end

% -------------------------------------------------------------------------
% Public Methods
% -------------------------------------------------------------------------
methods
    % -------------------------------------------------------------------------
    % Constructor
    % -------------------------------------------------------------------------
    function obj = TargetRegions(varargin)
        % Create the TargetRegions object
        % obh = TargetRegions() % Create empty class
        % obj = TargetRegions(..., 'geneName', name)
        % obj = TargetRegions(..., 'id', geneAccession)
        % obj = TargetRegions(..., 'geneSequence', seq)
        % obj = TargetRegions(..., 'startPos', posVec)
        % obj = TargetRegions(..., 'regionLength', lengthVec)
        % obj = TargetRegions(..., 'GC', gcVec)
        % obj = TargetRegions(..., 'Tm', tmVec)
        % obj = TargetRegions(..., 'specificity', specVec)
        % obj = TargetRegions(..., 'penalties', penMat)
        % obj = TargetRegions(..., 'penaltyNames', {name1, name2, ...})
        % obj = TargetRegions(..., 'sequence', parsedSeqs)
        
        % -------------------------------------------------------------------------
        % Handle empty object request
        % -------------------------------------------------------------------------
        if nargin < 1
            return;
        end
        
        % -------------------------------------------------------------------------
        % Parse variable inputs
        % -------------------------------------------------------------------------
        defaults = cell(0,3); 
        
        defaults(end+1,:) = {'geneName', 'string', ''}; % Display progress of construction
        defaults(end+1,:) = {'id', 'string', ''};
        defaults(end+1,:) = {'geneSequence', 'freeType', []};
        defaults(end+1,:) = {'startPos', 'positive', {}};
        defaults(end+1,:) = {'regionLength', 'positive', []};
        defaults(end+1,:) = {'GC', 'nonnegative', []};
        defaults(end+1,:) = {'Tm', 'positive', []};
        defaults(end+1,:) = {'specificity', 'nonnegative', []};
        defaults(end+1,:) = {'isoSpecificity', 'nonnegative', []};
        defaults(end+1,:) = {'penalties', 'nonnegative', []};
        defaults(end+1,:) = {'penaltyNames', 'cell', {}};
        defaults(end+1,:) = {'sequence', 'cell', {}};
        
        parameters = ParseVariableArguments(varargin, defaults, mfilename);
        
        % -------------------------------------------------------------------------
        % Store data internally
        % -------------------------------------------------------------------------
        fieldsToTransfer = setdiff(properties(obj), {'geneSequence', 'numRegions', 'map'});
        for i=1:length(fieldsToTransfer)
            obj.(fieldsToTransfer{i}) = parameters.(fieldsToTransfer{i});
        end
 
        obj.numRegions = length(obj.startPos);
        % -------------------------------------------------------------------------
        % Parse sequences if a gene sequence is provided as contiguous
        % sequence
        % -------------------------------------------------------------------------

        if ~iscell(parameters.geneSequence)
            obj.numRegions = length(obj.startPos);
            if ~isempty(parameters.geneSequence) && ischar(parameters.geneSequence)
                for i=1:obj.numRegions
                    obj.sequence{i} = parameters.geneSequence(...
                        obj.startPos(i):(obj.startPos(i) + obj.regionLength(i) - 1) );
                end
            end
            if ~isempty(parameters.geneSequence) && ~ischar(parameters.geneSequence)
                for i=1:obj.numRegions
                    obj.sequence{i} = obj.map(parameters.geneSequence(...
                        obj.startPos(i):(obj.startPos(i) + obj.regionLength(i) - 1) ) + 2);
                    % Note: 2 is added to map the unknown characters to 1. 
                end
            end
        else
            % Gene sequence is already cell array
            % No need to parse
            obj.sequence = parameters.geneSequence;
        end
    end    

    % -------------------------------------------------------------------------
    % Fasta write
    % -------------------------------------------------------------------------
    function fastawrite(obj, filePath, varargin)
        % Write the targetRegions to an existing or new fasta file
        % obj.fastawrite(filePath)
        % obj.fastawrite(filePath, 'overwrite', boolean)
        
        % -------------------------------------------------------------------------
        % Parse variable inputs
        % -------------------------------------------------------------------------
        defaults = cell(0,3); 
        defaults(end+1,:) = {'overwrite', 'boolean', false}; 
        
        parameters = ParseVariableArguments(varargin, defaults, mfilename);
        
        % -------------------------------------------------------------------------
        % Prepare file for ovewrite
        % -------------------------------------------------------------------------
        if exist(filePath)==2 && parameters.overwrite
            display(['Deleting existing file: ' filePath]);
            delete(filePath);
        end
        
        % -------------------------------------------------------------------------
        % Determine validity of path 
        % -------------------------------------------------------------------------
        fastawrite(filePath, {}, {});
        
        % -------------------------------------------------------------------------
        % Disable warnings 
        % -------------------------------------------------------------------------
        if length(obj) > 1
            warnState = warning;
            warning('off', 'bioinfo:fastawrite:AppendToFile');
        end
        
        % -------------------------------------------------------------------------
        % Format headers 
        % -------------------------------------------------------------------------
        fieldsToTransfer = setdiff(properties(obj), {'sequence', 'penalties', 'penaltyNames', 'numRegions'}); 
       
        for o=1:length(obj)
            for i=1:obj(o).numRegions
                localHeader = ['id=' obj(o).id ...
                    ' geneName=' obj(o).geneName ...
                    ' startPos=' num2str(obj(o).startPos(i)) ...
                    ' regionLength=' num2str(obj(o).regionLength(i)) ...
                    ' GC=' num2str(obj(o).GC(i)) ...
                    ' Tm=' num2str(obj(o).Tm(i)) ...
                    ' specificity=' num2str(obj(o).specificity(i)) ...
                    ' isoSpecificity=' num2str(obj(o).isoSpecificity(i)) ...
                    ];
                % Add penalties
                for j=1:length(obj(o).penaltyNames)
                    localHeader = [localHeader ' p_' obj(o).penaltyNames{j} '=' num2str(obj(o).penalties(j,i))];
                end
                headers{i} = strtrim(localHeader);
            end

            % -------------------------------------------------------------------------
            % Write fasta 
            % -------------------------------------------------------------------------
            fastawrite(filePath, headers, obj(o).sequence);
        end
        
        % -------------------------------------------------------------------------
        % Reenable warnings 
        % -------------------------------------------------------------------------
        if length(obj) > 1
            warning(warnState); % Restore previous warning state
        end
        
    end
    
    % -------------------------------------------------------------------------
    % Save Function
    % -------------------------------------------------------------------------
    function Save(obj, dirPath)
        % Save the transcriptome object in a directory specified by dirPath
        % obj.Save(dirPath)
        
        % -------------------------------------------------------------------------
        % Check directory validity
        % -------------------------------------------------------------------------
        if dirPath(end) ~= filesep
            dirPath(end+1) = filesep;
        end
        [status, message] = mkdir(dirPath);
        if ~status
            error('matlabFunctions:invalidArguments', 'Invalid directory path');
            return;
        end
        
        % -------------------------------------------------------------------------
        % Check if array
        % -------------------------------------------------------------------------
        numElement = length(obj);
        if numElement > 1
            display(['Saving array of ' num2str(numElement)]);
            SaveAsByteStream([dirPath 'numElement.matb'], numElement, 'verbose', false);
        end
        
        % -------------------------------------------------------------------------
        % Save fields
        % -------------------------------------------------------------------------
        fieldsToSave = properties(obj);
        for i=1:length(fieldsToSave)
            switch fieldsToSave{i}
                case {'GC', 'Tm', 'penalties', 'startPos', 'regionLength', 'numRegions', 'specificity', 'isoSpecificity'}
                    data = [obj.(fieldsToSave{i})];
                case {'id', 'geneName'}
                    data = {obj.(fieldsToSave{i})};
                case {'sequence'}
                    data = [obj.(fieldsToSave{i})];
                case {'penaltyNames'} % All elements in the array must have the same penalty fields
                    data = [obj(1).(fieldsToSave{i})];
            end
            if ~isempty(data)
                SaveAsByteStream([dirPath fieldsToSave{i} '.matb'], ...
                    data, 'verbose', false);
            end
        end
    end
end

% -------------------------------------------------------------------------
% Static methods
% -------------------------------------------------------------------------
methods (Static)
    % -------------------------------------------------------------------------
    % Build a TargetRegions object from a saved version
    % -------------------------------------------------------------------------
    function obj = Load(dirPath, varargin)
        % obj = TargetRegions.Load(dirPath)
        % obj = TargetRegions.Load(..., 'verbose', boolean)
        
        % -------------------------------------------------------------------------
        % Check provided path
        % -------------------------------------------------------------------------
        if dirPath(end) ~= filesep
            dirPath(end+1) = filesep;
        end
        if ~isdir(dirPath)
            error('matlabFunctions:invalidArguments', 'The provided path is not valid');
        end
        
        % -------------------------------------------------------------------------
        % Parse variable inputs
        % -------------------------------------------------------------------------
        defaults = cell(0,3); 
        defaults(end+1,:) = {'verbose', 'boolean', false}; % Display progress of construction
        parameters = ParseVariableArguments(varargin, defaults, mfilename);
        
        % -------------------------------------------------------------------------
        % Determine if it is multiple or single objects
        % -------------------------------------------------------------------------
        if exist([dirPath 'numElement.matb'])
            numElement = LoadByteStream([dirPath 'numElement.matb'], ...
                'verbose', parameters.verbose);
        end
        
        % -------------------------------------------------------------------------
        % Create empty object (to define fields to load)
        % -------------------------------------------------------------------------
        fieldsToLoad = properties(TargetRegions());
        
        % -------------------------------------------------------------------------
        % Load flat objects
        % -------------------------------------------------------------------------
        for i=1:length(fieldsToLoad)
            if exist([dirPath fieldsToLoad{i} '.matb'])
                data.(fieldsToLoad{i}) = LoadByteStream([dirPath fieldsToLoad{i} '.matb'], ...
                    'verbose', parameters.verbose);
            end
        end
        
        % -------------------------------------------------------------------------
        % Parse into objects
        % -------------------------------------------------------------------------
        index = [0 cumsum(data.numRegions)];
        foundFields = fields(data);
        for o=1:numElement
            % Compose input for TargetRegions
            variableInput = cell(0,0);
            for j=1:length(foundFields)
                % Create flag
                variableInput{end+1} = foundFields{j};
                
                % Parse loaded data
                switch foundFields{j}
                    case {'numRegions'} % Remove this flag, it is calculated by the constructor
                        variableInput = variableInput(1:(end-1));
                    case {'GC', 'Tm', 'penalties', 'startPos', 'regionLength', 'numRegions', 'sequence', 'specificity', 'isoSpecificity'}
                        variableInput{end+1} = data.(foundFields{j})((index(o)+1):index(o+1));
                    case {'id', 'geneName'}
                        variableInput{end+1} = data.(foundFields{j}){o};
                    case {'penaltyNames'}
                        variableInput{end+1} = data.(foundFields{j}); % All entries share the same names
                end
            end
            obj(o) = TargetRegions(variableInput{:});
        end
    end
    
end % Static methods

end
