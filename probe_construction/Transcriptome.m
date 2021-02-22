classdef Transcriptome < handle
% ------------------------------------------------------------------------
% [transcriptomeObj] = Transcriptome(sequenceData, varargin)
% This class contains information about a transcriptome, including
% sequences, gene names, gene IDs, expression data, and isoform
% information.
%--------------------------------------------------------------------------
% Necessary Inputs
% targetSequences -- Either the path to a fasta file that contains target
% names or a structure array that has the same organization as the output
% of fastaread. 
%--------------------------------------------------------------------------
% Methods
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
% 'abundPath' -- A path to a cufflinks isoforms.fpkm_tracking file.
% 'headerType' -- Reserved for future expansion of abundance data types
% 'IDType' -- Reserved for future expansion of isoform/gene ID processing
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% April 27, 2015
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Define properties
% -------------------------------------------------------------------------
properties
    numTranscripts      % The number of loaded transcripts
    numGenes            % The number of loaded genes
    verbose             % The verbose status of the class
    headerType          % The type of the abundance data file -- reserved for future use
    IDType              % The type of the transcript ID -- reserved for future use
    abundPath           % The path to the abundance data file
    abundLoaded         % A boolean that determines if abundance data are loaded
    transPath           % The path to the transcriptome file, if provided
end

properties (Hidden=true, SetAccess=protected)
    ids                 % Accession values for each entry
    geneNames           % Common name for each entry
    intSequences        % Integer form for the sequence of each entry
    abundance           % Abundance value for each entry, if loaded
    id2Ind              % A map for conversion of ID to internal index
    name2Ind            % A map for converstion of gene name to internal index
end

% -------------------------------------------------------------------------
% Define methods
% -------------------------------------------------------------------------
methods
    % -------------------------------------------------------------------------
    % Define constructor
    % -------------------------------------------------------------------------
    function obj = Transcriptome(transcriptome, varargin)
        % Create the transcriptome object
        % obj = Transcriptome(fastaread(...))
        % obj = Transcriptome(transcriptomePath)
        % obj = Transcriptome(..., 'abundPath', pathToAbundanceFile)
        % obj = Transcriptome(..., 'verbose', true/false);
        % obj = Transcriptome({ids, geneNames, sequences, abund}, ...) %
       
        % -------------------------------------------------------------------------
        % Default variables
        % -------------------------------------------------------------------------
        defaults = cell(0,3);
        
        defaults(end+1,:) = {'verbose', 'boolean', false}; % Display progress of construction
        defaults(end+1,:) = {'headerType', {'cufflinks', 'ensembl', 'custom', 'refseq'}, 'cufflinks'}; % Assume a cufflinks format
        defaults(end+1,:) = {'IDType', 'string', 'NCBI'}; % The type of transcript ID, e.g. ENSEMBL or NCBI. 
        defaults(end+1,:) = {'abundPath', 'filePath', []}; % Path to abundance data
        
        % -------------------------------------------------------------------------
        % Parse variable input for constructor class
        % -------------------------------------------------------------------------
        parameters = ParseVariableArguments(varargin, defaults, mfilename);
        
        % Transfer parameters
        f = intersect(fields(parameters), properties(obj));        
        for i=1:length(f)
            obj.(f{i}) = parameters.(f{i});
        end
        
        % -------------------------------------------------------------------------
        % Initialize empty values for everything
        % -------------------------------------------------------------------------
        obj.numTranscripts = 0;
        obj.numGenes = 0;
        obj.abundLoaded = false;
        obj.transPath = '';
        
        % -------------------------------------------------------------------------
        % Parse necessary input
        % -------------------------------------------------------------------------
        if nargin < 1 % Define empty class
            return;
        elseif ischar(transcriptome) % A path was provided
            if exist(transcriptome) == 2 % Does this path exist
                if obj.verbose
                    display('----------------------------------------------------------')
                    display(['Loading transcriptome from ' transcriptome]);
                    tic;
                end
                obj.transPath = transcriptome;
                transcriptome = fastaread(transcriptome);
                if obj.verbose
                    display(['... completed in ' num2str(toc) ' s']);
                    display(['Found ' num2str(length(transcriptome)) ' sequences']);
                end
            else
                error('matlabFunctions:invalidArguments', 'Invalid path to target sequences');
            end
        elseif isstruct(transcriptome) % A fasta structure was provided
            if ~isempty(setdiff({'Header', 'Sequence'}, fields(transcriptome)))
                error('matlabFunctions:invalidArguments', 'Provided structure does not have all required fields');
            end
        elseif iscell(transcriptome) && length(transcriptome) == 4 % All data were provided in a cell
            % Insert everything (but sequences)
            
%             assignin('base', 'transObj', transcriptome);
            
            obj.ids = transcriptome{1};
            obj.geneNames = transcriptome{2};
            obj.abundance = transcriptome{4};
            obj.abundLoaded = false;
            % Update abundance loaded flag
            if ~isempty(obj.abundance)
                obj.abundLoaded = true;
            end
           
            % Check sequences for integer versus character
            if ~isa(transcriptome{3}{1}, 'int8')
                obj.intSequences = cellfun(@(x) int8(nt2int(x)) - 1, ...
                    transcriptome{3}, 'UniformOutput', false);
            else
                obj.intSequences = transcriptome{3};
            end
        else
            error('matlabFunctions:invalidArguments', 'Invalid input');
        end

        % -------------------------------------------------------------------------
        % Parse Headers into IDs and gene names
        % -------------------------------------------------------------------------
        if isempty(obj.ids) % Handle the direct construction
            if obj.verbose
                PageBreak();
                display(['Parsing transcriptome']);
                tic;
            end
            switch obj.headerType
                case 'cufflinks'
                    parseFunc = @(x)regexp(x, '(?<id>\S*) (?<name>gene=\S*)', 'names');
                    for i=1:length(transcriptome)
                        parsedData = parseFunc(transcriptome(i).Header);
                        obj.ids{i} = parsedData.id;
                        obj.geneNames{i} = parsedData.name(6:end); % Strip off gene=
                        obj.intSequences{i} = int8(nt2int(transcriptome(i).Sequence)) - 1; % 0=A, 1=C, ... >3 = not valid
                    end
                case 'ensembl'
                    parseFunc = @(x)regexp(x, '(?<id>^ENS\D+T\d+[.]\d).*(?<geneID>gene:ENS\D+G\d*[.]\d).*(?<name>gene_symbol:\S+)', 'names');
                    for i=1:length(transcriptome)
                        
                        try
                        
                            parsedData = parseFunc(transcriptome(i).Header);
                            obj.ids{i} = parsedData.id;
    %                         obj.geneID{i} = parsedData.geneID(6:end); % Strip
    %                         off gene 
                            % Might be unused elsewhere, so remove here.
                            obj.geneNames{i} = parsedData.name((numel('gene_symbol:')+1):end); 
                            obj.intSequences{i} = int8(nt2int(transcriptome(i).Sequence)) - 1;
                            
                        catch 
                            
                            fprintf(1, 'Failed on transcriptome i = %d\n', i);
                            assignin('base', 'parsedData', parsedData);
                            assignin('base', 'header', transcriptome(i).Header);
                            
                        end
                    end
                    
                case 'refseq'
                    
   
                    % Brute-force name convention cases here
                    parseFuncA = @(x)regexp(x, '(?<id>^NM_\d+[.]\d).*(?<name>\(\S+\))', 'names');
                    parseFuncB = @(x)regexp(x, '(?<id>^XM_\d+[.]\d).*(?<name>\(\S+\))', 'names');
                    parseFuncC = @(x)regexp(x, '(?<id>^NR_\d+[.]\d).*(?<name>\(\S+\))', 'names');
                    parseFuncD = @(x)regexp(x, '(?<id>^XR_\d+[.]\d).*(?<name>\(\S+\))', 'names');
                    for i = 1:length(transcriptome)
                        
                        if transcriptome(i).Header(1) == 'X'
                            if transcriptome(i).Header(2) == 'M'
                                parseFunc = parseFuncB;
                            
                            elseif transcriptome(i).Header(2) == 'R'
                                parseFunc = parseFuncD;
                                
                            else
                                fprintf(1, 'Transcriptome i = %d has unexpected header\n', i);
                                assignin('base', 'parsedData', parsedData);
                                assignin('base', 'header', transcriptome(i).Header);
                            end
                            
                        elseif transcriptome(i).Header(1) == 'N'
                            if transcriptome(i).Header(2) == 'M'
                                parseFunc = parseFuncA;
                            
                            elseif transcriptome(i).Header(2) == 'R'
                                parseFunc = parseFuncC;
                                
                            else
                                fprintf(1, 'Transcriptome i = %d has unexpected header\n', i);
                                assignin('base', 'parsedData', parsedData);
                                assignin('base', 'header', transcriptome(i).Header);
                            end
                            
                        else
                            fprintf(1, 'Transcriptome i = %d has unexpected header\n', i);
%                             assignin('base', 'parsedData', parsedData);
                            assignin('base', 'header', transcriptome(i).Header);
                                
                        end
                        
                        try
                            parsedData = parseFunc(transcriptome(i).Header);
                            obj.ids{i} = parsedData.id;
                            obj.geneNames{i} = parsedData.name((1 + strfind(parsedData.name, '(')):(strfind(parsedData.name, ')') - 1));
                            obj.intSequences{i} = int8(nt2int(transcriptome(i).Sequence)) - 1;
                        catch 
                            fprintf(1, 'Failed on transcriptome i = %d\n', i);
%                             assignin('base', 'parsedData', parsedData);
                            assignin('base', 'header', transcriptome(i).Header);
                        end
                    end
            end
            if obj.verbose
                display(['... completed in ' num2str(toc) ' s']);
            end
        end

        % -------------------------------------------------------------------------
        % Build internal index maps: ID to index
        % -------------------------------------------------------------------------
        obj.numTranscripts = length(obj.ids);
        obj.id2Ind = containers.Map(obj.ids, 1:obj.numTranscripts);
        
        % -------------------------------------------------------------------------
        % Build internal index maps: gene name to index
        % -------------------------------------------------------------------------
        [geneNames, ~, ic] = unique(obj.geneNames);
        isoFormInds = {};
        for i=1:length(geneNames)
            isoFormInds{i} = find(ic==i);
        end
        obj.name2Ind = containers.Map(geneNames, isoFormInds);
        obj.numGenes = length(geneNames);
        
        % -------------------------------------------------------------------------
        % Set abundances
        % -------------------------------------------------------------------------
        if isempty(obj.abundance)
            if ~isempty(parameters.abundPath)
                obj.AddAbundances(parameters.abundPath)
            else
                obj.abundance = zeros(1, obj.numTranscripts);
            end
        end
    end

    % -------------------------------------------------------------------------
    % AddAbundances
    % -------------------------------------------------------------------------
    function [notIn, notIncluded] = AddAbundances(obj, varargin)
        % Add abundances to transcriptome
        % obj.AddAbundances(pathToAbundanceFile)
        % obj.AddAbundances(idNames, abundVec)
        % Currently only supports cufflinks .fpkm_tracking files
        
        % -------------------------------------------------------------------------
        % Check file type and path
        % -------------------------------------------------------------------------
        if ischar(varargin{1})
            abundancePath = varargin{1};
            if ~exist(abundancePath) == 2
                error('matlabFunctions:invalidArguments', 'Invalid path to abundance data');
            end
            [~, ~, ext] = fileparts(abundancePath);
            if ~strcmp(ext,'.fpkm_tracking')
                warning('Only fpkm_tracking files are currently supported.');
            end

            % -------------------------------------------------------------------------
            % Load file and extract ids and fpkm
            % -------------------------------------------------------------------------
            if obj.verbose
                display(['Loading abundances from ' abundancePath]);
                tic;
            end
            fid = fopen(abundancePath);
            fileData = textscan(fid, '%s%s%s%s%s%s%s%d%f%f%f%f%s', 'HeaderLines', 1);
            fclose(fid);

            foundIds = fileData{1};
            fpkm = fileData{10};
        else nargin == 2 & all(varargin{2}>=0) & iscell(varargin{1}) 
            % -------------------------------------------------------------------------
            % Handle the case of direct loading of abundance names
            % -------------------------------------------------------------------------
            foundIds = varargin{1};
            fpkm = varargin{2};
        end
        
        % -------------------------------------------------------------------------
        % Clear abundances and load new data 
        % -------------------------------------------------------------------------
        obj.abundance = zeros(1, obj.numTranscripts);
        
        [commonIds, ia, ib] = intersect(foundIds, obj.ids);
        obj.abundance(ib) = fpkm(ia);
        
        obj.abundLoaded = true;
        
        % -------------------------------------------------------------------------
        % Return ids that were not included or provided
        % -------------------------------------------------------------------------
        if nargout > 1
            notIncluded = setdiff(obj.ids, commonIds);
        end
        if nargout > 0
            notIn = setdiff(foundIds, commonIds);
        end
        
        if obj.verbose
            display(['... completed in ' num2str(toc) ' s']);
        end
    end
    
    % -------------------------------------------------------------------------
    % GetAbundancesById -- Return the abundance of a set of ids 
    % -------------------------------------------------------------------------
    function abund = GetAbundanceByID(obj, ids)
        % Return a list of abundances corresponding to the ids provided
        % abund = obj.GetAbundanceByID(ids)
        % ids not in the transcriptome are returned as Nan

        validKeys = isKey(obj.id2Ind, ids);
        abund = nan(1, length(ids));
        keys = obj.id2Ind.values(ids(validKeys));
        keys = cat(1, keys{:});
        abund(validKeys) = obj.abundance(keys);
    end

    % -------------------------------------------------------------------------
    % GetAbundanceByName -- Return the abundance of a gene 
    % -------------------------------------------------------------------------
    function abund = GetAbundanceByName(obj, names, varargin)
        % Return a list of abundances corresponding to provided names
        % abund = obj.GetAbundanceByName(names)
        % abund = obj.GetAbundanceByName(names, 'isoformFunc', functionHandle)
        % The abundance of multiple isoforms is added by default, but a
        % different function, e.g. @mean/@min/@max, can be passed with
        % 'isoformFunc'
        % abund = obj.GetAbundanceByName(names, 'returnType',
        % 'allIsoforms'/'isoformFunction') If 'allIsoforms' is specified
        % abundances are returned as a cellarray of isoform abundances.
        
        % -------------------------------------------------------------------------
        % Handle variable input
        % -------------------------------------------------------------------------
        defaults = cell(0,3);
        defaults(end+1,:) = {'returnType', {'allIsoforms', 'isoformFunction'}, 'isoformFunction'};
        defaults(end+1,:) = {'isoformFunc', 'function', @sum}; % How to report abundance with different isoforms        
        parameters = ParseVariableArguments(varargin, defaults, mfilename);
        
        % -------------------------------------------------------------------------
        % Find isoforms and associated keys
        % -------------------------------------------------------------------------
        if ~iscell(names) % Handle single input
            names = {names};
        end
        
        % Find valid keys
        validKeys = isKey(obj.name2Ind, names);
        inds = cell(1, length(validKeys));
        inds(validKeys) = obj.name2Ind.values(names(validKeys));
        
        % -------------------------------------------------------------------------
        % Accumulate isoform abundance information
        % -------------------------------------------------------------------------
        switch parameters.returnType
            case 'isoformFunction'
                % Create index for accum array
                flatInds = cat(1, inds{:});
                l = cellfun(@length, inds);
                accumInds = repelem(1:length(inds), l);

                % Get flat abundances
                abund = obj.abundance(flatInds);

                % Accumulate abundances for different isoforms
                abund = accumarray(accumInds', abund, [], parameters.isoformFunc);
            case 'allIsoforms'
                % Prepare abundance cell array
                abund = cell(0, length(inds));
                
                % Fill
                for i=1:length(inds)
                    abund{i} = obj.abundance(inds{i});
                end
        end
    end
    
    % -------------------------------------------------------------------------
    % GetSequencesByName: Return sequences by gene name
    % -------------------------------------------------------------------------
    function seqs = GetSequencesByName(obj, names)
        % Return all sequences for the specified name
        % seqs = obj.GetSequencesByName(names)
        % Broken. Requires fix to return correct sequence.
        % Does not support names with more than 1 transcript in object.
        % If multiple transcripts assigned to name, multiple values
        % returned by values(obj.name2Ind, names(isValid)); 
        % generating error.
        
        % Handle a single input
        wasSingleInput = false;
        if ~iscell(names)
            names = {names}; 
            wasSingleInput = true;
        end
        
        % Create sequences cell array
        % Must be corrected to allow case where multiple sequences can share same
        % name.
        
        seqs = cell(1, length(names));

        % Find valid names
        isValid = isKey(obj.name2Ind, names);
        
        % Assign into base workspace for troubleshooting.
        assignin('base', 'name2Ind', obj.name2Ind);
        assignin('base', 'isValid', isValid);
        assignin('base', 'intSequences', obj.intSequences);
        assignin('base', 'values', values);
        
        seqs(isValid) = obj.intSequences(cell2mat(values(obj.name2Ind, names(isValid))));
        
        % Convert sequences
        for i=1:length(seqs)
            if isValid(i)
                seqs{i} = int2nt(seqs{i}+1);
            end
        end
        
        % If was a single input remove cell wrapper
        if wasSingleInput
            seqs = seqs{1};
        end
    end
    
    % -------------------------------------------------------------------------
    % GetSequenceByID: Return sequences by gene id
    % -------------------------------------------------------------------------
    function seqs = GetSequenceByID(obj, ids)
        % Return sequence for the specified id
        % seq = obj.GetSequenceByID(ids)
        
        % Handle a single input
        wasSingleInput = false;
        if ~iscell(ids)
            ids = {ids}; 
            wasSingleInput = true;
        end

        % Create sequences cell array
        seqs = cell(1, length(ids));

        % Find valid names
        isValid = isKey(obj.id2Ind, ids);
        
        seqs(isValid) = obj.intSequences(cell2mat(values(obj.id2Ind, ids(isValid))));

        % Convert to nt
        seqs = cellfun(@(x)int2nt(x+1), seqs, 'UniformOutput', false);
                
        % If was a single input remove cell wrapper
        if wasSingleInput
            seqs = seqs{1};
        end
    end

    % -------------------------------------------------------------------------
    % GetIDsByName: Return gene ids by gene name(s)
    % -------------------------------------------------------------------------
    function ids = GetIDsByName(obj, names)
        % Return all transcript ids corresponding to a gene name or names
        % ids = obj.GetIDsByName(names)
        
        % Coerce input to cell
        if ~iscell(names)
            names = {names};
        end
        
        % Prepare output: empty cells indicates invalid names
        ids = cell(0, length(names));

        % Identify valid names
        validKeys = isKey(obj.name2Ind, names);
        
        if any(validKeys)
            validKeys = find(validKeys);
            inds = obj.name2Ind.values(names(validKeys));
            for i=1:length(inds)
                ids{validKeys(i)} = obj.ids(inds{i});
            end
        end
    end
    
    % -------------------------------------------------------------------------
    % GetNameById
    % -------------------------------------------------------------------------
    function name = GetNameById(obj, id)
        % Return gene name for a given id
        % ids = obj.GetIDsByName(name)

        % FUTURE WORK: Update to allow multiple id inputs
        if ~isKey(obj.id2Ind, id)
            name = '';
            return;
        end
        name = obj.geneNames(obj.id2Ind(id));
    end
    
    % -------------------------------------------------------------------------
    % GetNameById
    % -------------------------------------------------------------------------
    function names = GetNames(obj)
        % Return all gene names
        % names = obj.GetNames()

        names = keys(obj.name2Ind);
    end
    
    %-------------------------------------------------------------------------
    % Slice the transcriptome
    % -------------------------------------------------------------------------
    function newTranscriptome = Slice(obj, varargin)
        % Generate a new transcriptome object with a subset of entries
        % newObj = obj.Slice('geneID', {id1, id2, ...})
        
        % -------------------------------------------------------------------------
        % Handle variable input
        % -------------------------------------------------------------------------
        defaults = cell(0,3);
        defaults(end+1,:) = {'geneID', 'cell', {}};
        defaults(end+1,:) = {'geneName', 'cell', {}};       
        parameters = ParseVariableArguments(varargin, defaults, mfilename);

        % -------------------------------------------------------------------------
        % Get internal inds to slice
        % -------------------------------------------------------------------------
        if ~isempty(parameters.geneID)
            inds = obj.GetInternalInds('geneID', parameters.geneID);
        elseif ~isempty(parameters.geneName)
            inds = obj.GetInternalInds('geneName', parameters.geneName);
        else
            inds = 1:obj.numTranscripts;
        end
        
        % -------------------------------------------------------------------------
        % Build new transcriptome object
        % -------------------------------------------------------------------------
        newTranscriptome = Transcriptome(...
            {obj.ids(inds), obj.geneNames(inds), obj.intSequences(inds), obj.abundance(inds)}, ...
            'verbose', obj.verbose, ...
            'headerType', obj.headerType, ...
            'IDType', obj.IDType, ...
            'abundPath', obj.abundPath);        
    end
    
    % -------------------------------------------------------------------------
    % Save Function
    % -------------------------------------------------------------------------
    function Save(obj, dirPath)
        % Save the Transcriptome object in a directory specified by dirPath
        % obj.Save(dirPath)
        
        % -------------------------------------------------------------------------
        % Check directory validity
        % -------------------------------------------------------------------------
        if dirPath(end) ~= filesep
            dirPath(end+1) = filesep;
        end
        if ~isdir(dirPath)
            [status, message] = mkdir(dirPath);
            if ~status
                error('matlabFunctions:invalidArguments', 'Invalid directory path');
                return;
            end
        end

        % -------------------------------------------------------------------------
        % Define fields to save
        % -------------------------------------------------------------------------
        fieldsToSave = properties(obj);
        fieldsToSave = union(fieldsToSave, {'ids', 'geneNames', ...
            'intSequences', 'abundance', 'id2Ind', 'name2Ind'});
    
        % -------------------------------------------------------------------------
        % Save fields
        % -------------------------------------------------------------------------
        for i=1:length(fieldsToSave)
            SaveAsByteStream([dirPath fieldsToSave{i} '.matb'], ...
                obj.(fieldsToSave{i}), 'verbose', obj.verbose);
        end
    end
end

% -------------------------------------------------------------------------
% Hidden methods
% -------------------------------------------------------------------------
methods (Hidden=true)
    % -------------------------------------------------------------------------
    % Return internal index for gene names or gene ids 
    % -------------------------------------------------------------------------
    function [inds, ids, names, validKeys] = GetInternalInds(obj, varargin)
        % Return the internal order of elements in the transcriptome as
        % specified by a list of gene names or gene ids
        % inds = obj.GetInternalInds('geneName', names)
        % inds = obj.GetInternalInds('geneID', ids)
        
        % -------------------------------------------------------------------------
        % If nothing is requested return everything 
        % -------------------------------------------------------------------------
        if nargin < 1
            inds = 1:obj.numTranscripts;
            names = obj.geneNames;
            ids = obj.ids;
            validKeys = true(1, length(inds));
        end
        
        % -------------------------------------------------------------------------
        % Parse variable inputs
        % -------------------------------------------------------------------------
        defaults = cell(0,3); 
        defaults(end+1,:) = {'geneName', 'cell', {}}; % Display progress of construction
        defaults(end+1,:) = {'geneID', 'cell', {}};
        defaults(end+1,:) = {'inds', 'positive',[]};
        
        parameters = ParseVariableArguments(varargin, defaults, mfilename);
        
        % -------------------------------------------------------------------------
        % Return nothing if no valid inputs
        % -------------------------------------------------------------------------
        inds = [];
        names = {};
        ids = {};

        % -------------------------------------------------------------------------
        % Return desired properties
        % -------------------------------------------------------------------------
        if ~isempty(parameters.inds)
            inds = parameters.inds;
            names = obj.geneNames(inds);
            ids = obj.ids(inds);
        elseif ~isempty(parameters.geneName)
            validKeys = obj.name2Ind.isKey(parameters.geneName);
            if any(validKeys)
                inds = obj.name2Ind.values(parameters.geneName(validKeys));
                inds = cat(1, inds{:}); % Flatten
                names = obj.geneNames(inds);
                ids = obj.ids(inds);
            end
        elseif ~isempty(parameters.geneID)
            validKeys = obj.id2Ind.isKey(parameters.geneID);
            if any(validKeys)
                inds = obj.id2Ind.values(parameters.geneID(validKeys));
                inds = cat(1, inds{:});
                names = obj.geneNames(inds);
                ids = obj.ids(inds);
            end
        else % If nothing is requested, return everything
            inds = 1:obj.numTranscripts;
            names = obj.geneNames;
            ids = obj.ids;
            validKeys = true(1, length(inds));
        end
    end
end

% -------------------------------------------------------------------------
% Static methods
% -------------------------------------------------------------------------
methods (Static)
    % -------------------------------------------------------------------------
    % Build a Transcriptome object from a saved version
    % -------------------------------------------------------------------------
    function obj = Load(dirPath, varargin)
        % obj = Transcriptome.Load(dirPath, 'verbose', boolean)
        
        % -------------------------------------------------------------------------
        % Handle variable arguments
        % -------------------------------------------------------------------------
        defaults = cell(0,3);
        defaults(end+1,:) = {'verbose', 'boolean', false};
        parameters = ParseVariableArguments(varargin, defaults, mfilename);
        
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
        % Create empty object (to define fields to load)
        % -------------------------------------------------------------------------
        obj = Transcriptome();
        
        % -------------------------------------------------------------------------
        % Define fields to load
        % -------------------------------------------------------------------------
        fieldsToLoad = properties(obj);
        fieldsToLoad = union(fieldsToLoad, {'ids', 'geneNames', ...
            'intSequences', 'abundance', 'id2Ind', 'name2Ind'});
        % -------------------------------------------------------------------------
        % Load properties/data
        % -------------------------------------------------------------------------
        for i=1:length(fieldsToLoad)
            obj.(fieldsToLoad{i}) = LoadByteStream([dirPath fieldsToLoad{i} '.matb'], ...
            'verbose', parameters.verbose);
        end
    end
end

end
