classdef TRDesigner < handle
% ------------------------------------------------------------------------
% [trDesignerObj] = TRDesigner(varargin)
% This class designs target regions given a transcriptome object. 
%
% See Transcriptome and PrimerDesigner
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% April 2018
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2018.
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Define properties
% -------------------------------------------------------------------------
properties 
    verbose 	% A boolean that determines whether the class displays progress
end

properties (SetAccess=protected)
    transcriptome           % A transcriptome object that defines the sequences to target
    parallel                % A parallel.Pool object
    numPar                  % The number of parallel workers
    
    dG                      % The free energies associated with all nearest neighbours in all sequences
    
    gc                      % The GC content of all sequences
    isValid                 % The validity of all sequences, i.e. is each base A,C,T,U, or G
    
    OTTableNames            % Names of the off-target tables
    OTTables                % An array of off-target tables
    penalties               % The penalties associated with each sequence for each table
    
    forbiddenSeqs           % The list of forbidden sequences
    isForbiddenSeq          % The presence of forbidden sequences
   
    specificity             % The fraction of each sequence that is unique to that sequence
    specificityTable        % The OTTable used to calculate region specificity
    
    isoSpecificity          % The fraction of each sequence that is unique to the isoforms of that gene
    isoSpecificityTables    % An array of OTTables used to calculate isoform specificity
end

% -------------------------------------------------------------------------
% Public Methods
% -------------------------------------------------------------------------
methods
    % -------------------------------------------------------------------------
    % Define constructor
    % -------------------------------------------------------------------------
    function obj = TRDesigner(varargin)
        % Create the TRDesigner object
        % obj = TRDesigner('transcriptome', transcriptomeObject)
        % obj = TRDesigner(..., 'OTTables', {OTTableObj1, OTTableObj2, ...})
        % obj = TRDesigner(..., 'OTTableNames', {tableName1, tableName2, ...})
        % obj = TRDesigner(..., 'parallel', parallel.PoolObj)
        % obj = TRDesigner(..., 'verbose', boolean)
        % obj = TRDesigner(..., 'forbiddenSeqs', {seq1, seq2, ...});
        % obj = TRDesigner(..., 'specificityTable', OTTableObj);
        % obj = TRDesigner(..., 'isoSpecificityTables', OTTableObjArray);
        
        % -------------------------------------------------------------------------
        % Parse variable inputs
        % -------------------------------------------------------------------------
        % Define defaults
        defaults = cell(0,3); 
        defaults(end+1,:) = {'verbose', 'boolean', true}; % Display progress of construction
        defaults(end+1,:) = {'transcriptome', 'freeType', []};
        defaults(end+1,:) = {'OTTables', 'freeType', []};
        defaults(end+1,:) = {'OTTableNames', 'cell', {}};
        defaults(end+1,:) = {'parallel', 'parallel', []};
        defaults(end+1,:) = {'forbiddenSeqs', 'cell', {}};
        defaults(end+1,:) = {'specificityTable', 'freeType', []};
        defaults(end+1,:) = {'isoSpecificityTables', 'freeType', []};

        % Parse variable input
        parameters = ParseVariableArguments(varargin, defaults, mfilename);
        
        % Transfer values to object
        f = setdiff(fields(parameters), {'parallel', 'OTTableNames', 'OTTables', 'forbiddenSeqs'});
        for i=1:length(f)
            obj.(f{i}) = parameters.(f{i});
        end
        % Set parallel pool object
        SetParallel(obj, parameters.parallel); 
        
        % -------------------------------------------------------------------------
        % Initialize all parameters for empty constructor case
        % -------------------------------------------------------------------------
        obj.gc = {};
        obj.dG = {};
        obj.isValid = {};
        obj.isForbiddenSeq = {};
        obj.specificity = {};
        obj.isoSpecificity = {};
        if nargin < 1
            return;
        end
        
        % -------------------------------------------------------------------------
        % Check input values
        % -------------------------------------------------------------------------
        if ~isempty(obj.transcriptome) && ~strcmp(class(obj.transcriptome), 'Transcriptome')
            error('matlabFunctions:invalidArguments', 'Provided transcriptome must be a transcriptome object');
        end
        if ~isempty(parameters.OTTables) && length(parameters.OTTables)~=length(parameters.OTTableNames)
            error('matlabFunctions:invalidArguments', 'An equal number of names as OT tables must be provided.');
        end

        % -------------------------------------------------------------------------
        % Create annotations for penalty table 
        % -------------------------------------------------------------------------
        for i=1:length(parameters.OTTables)
            obj.AddOTTable(parameters.OTTables(i), parameters.OTTableNames{i});
        end
        
        % -------------------------------------------------------------------------
        % Add specificity table
        % -------------------------------------------------------------------------
        if ~isempty(parameters.specificityTable) || ~isempty(parameters.isoSpecificityTables)
            obj.AddSpecificityTable(parameters.specificityTable, parameters.isoSpecificityTables);
        end
        
        % -------------------------------------------------------------------------
        % Create annotations for GC and Tm calculations
        % -------------------------------------------------------------------------
        seqs = obj.transcriptome.intSequences;
        gc = {};
        valid = {};
        parfor (s=1:length(seqs), obj.numPar) % Loop over sequences in parallel
            gc{s} = seqs{s} == 1 | seqs{s} == 2;        % Determine G/C
            valid{s} = seqs{s} <= 3 & seqs{s} >=0;      % Determine validity
            dG{s} = TRDesigner.SantaLuciaNearestNeighbor(seqs{s});  % Calculate thermodynamic properties
        end
        % Store values
        obj.gc = gc;
        obj.dG = dG;
        obj.isValid = valid;

        % -------------------------------------------------------------------------
        % Create forbidden sequence annotations  
        % -------------------------------------------------------------------------
        for i=1:length(parameters.forbiddenSeqs)
            obj.AddForbiddenSeq(parameters.forbiddenSeqs{i});
        end        
    end
    
    % -------------------------------------------------------------------------
    % Add a forbidden sequence
    % -------------------------------------------------------------------------
    function AddForbiddenSeq(obj, seq, varargin)
        % Add a forbidden sequence or replace existing sequences
        % obj.AddForbiddenSeq(sequence)
        % obj.AddForbbidenSeq(..., 'replace', boolean) -- Replace the
        %   existing sequences
        
        % -------------------------------------------------------------------------
        % Parse necessary input
        % -------------------------------------------------------------------------
        if nargin < 1 || ~(ischar(seq) || all(seq >= 0 & seq <= 3))
            error('matlabFunctions:invalidArguments', 'The provided sequence is invalid.');
        end

        % -------------------------------------------------------------------------
        % Handle variable arguments
        % -------------------------------------------------------------------------
        defaults = cell(0,3);
        defaults(1,:) = {'replace', 'boolean', false}; % Either add to existing tables (false) or replace (true)
        parameters = ParseVariableArguments(varargin, defaults, mfilename);
     
        % -------------------------------------------------------------------------
        % Reset sequences if requested
        % -------------------------------------------------------------------------
        if parameters.replace
            obj.forbiddenSeqs = {};        
            obj.isForbiddenSeq = {};          
        end
        
        % -------------------------------------------------------------------------
        % Add sequence
        % -------------------------------------------------------------------------
        if ischar(seq)
            seq = int8(nt2int(seq, 'ACGTOnly', true)) - 1;
        end
        obj.forbiddenSeqs{end+1} = seq;
        
        % -------------------------------------------------------------------------
        % Display progress
        % -------------------------------------------------------------------------
        if obj.verbose
            PageBreak();
            tic;
            display(['Finding all occurrences of ' int2nt(seq+1)]);
        end
        
        % -------------------------------------------------------------------------
        % Calculate penalty
        % -------------------------------------------------------------------------
        hashBase = fliplr(4.^[0:(length(seq)-1)]);
        forbiddenHash = sum(hashBase.*double(seq));
        seqs = obj.transcriptome.intSequences; 

        parfor (s=1:length(seqs), obj.numPar) % Loop over transcriptome
            % Hash sequence
            seqHash = filter(hashBase, 1, double(seqs{s}));
            seqHash = seqHash(length(hashBase):end);
                        
            % Find forbidden sequences
            hasForbiddenSeq{s} = seqHash == forbiddenHash;
        end
        
        obj.isForbiddenSeq{end+1} = hasForbiddenSeq;
        
        % -------------------------------------------------------------------------
        % Display progress
        % -------------------------------------------------------------------------
        if obj.verbose
            display(['... completed in ' num2str(toc) ' at ' datestr(now)]);
        end
    end
    
    % -------------------------------------------------------------------------
    % Add a penalty table and calculate the appropriate penalties
    % -------------------------------------------------------------------------
    function AddOTTable(obj, otTable, tableName, varargin)
        % Add or modify an existing penalty table
        % obj.AddOTTable(obj, OTTable, OTTableName)
        % obj.AddOTTable(..., 'replace', boolean) -- replace existing
        %  tables
        
        % -------------------------------------------------------------------------
        % Parse necessary input
        % -------------------------------------------------------------------------
        if nargin < 2 || ~isa(otTable, 'OTTable') || ~ischar(tableName)
            error('matlabFunctions:invalidArguments', 'Both an OTTable and a name must be provided');
        end
        
        % -------------------------------------------------------------------------
        % Handle variable arguments
        % -------------------------------------------------------------------------
        defaults = cell(0,3);
        defaults(1,:) = {'replace', 'boolean', false}; % Either add to existing tables (false) or replace (true)
        parameters = ParseVariableArguments(varargin, defaults, mfilename);
        
        % -------------------------------------------------------------------------
        % Reset tables if requested
        % -------------------------------------------------------------------------
        if parameters.replace
            obj.OTTableNames = {};           
            obj.OTTables = [];
            obj.penalties = {};
        end
        
        % -------------------------------------------------------------------------
        % Add table and name
        % -------------------------------------------------------------------------
        if isempty(obj.OTTables) % Handle initial table/name pair
            obj.OTTables = otTable;
            obj.OTTableNames = {tableName};
        else
            obj.OTTables(end+1) = otTable;
            obj.OTTableNames{end+1} = tableName;
        end
        
        % -------------------------------------------------------------------------
        % Display progress
        % -------------------------------------------------------------------------
        if obj.verbose
            PageBreak();
            tic;
            display(['Calculating penalty values for ' obj.OTTableNames{end} ': ' datestr(now)]);
        end
        
        % -------------------------------------------------------------------------
        % Build penalty
        % -------------------------------------------------------------------------
        % Make a local copy of the sequences
        seqs = obj.transcriptome.intSequences; 
                        
        % Make reference to table
        table = obj.OTTables(end);
        
        % Loop over all sequences
        for s=1:length(seqs)
            penalty{s} = table.CalculatePenalty(seqs{s}); % Use table to calculate penalty
        end

        % Save penalty
        obj.penalties{end+1} = penalty; 
        
        % -------------------------------------------------------------------------
        % Display progress
        % -------------------------------------------------------------------------
        if obj.verbose
            display(['... completed in ' num2str(toc) ' at ' datestr(now)]);
        end
        
    end
        
    % -------------------------------------------------------------------------
    % Add a table for determining region specificity
    % -------------------------------------------------------------------------
    function AddSpecificityTable(obj, specificityTable, isoSpecificityTables)
        % Add or replace a specificity table
        % obj.AddSpecificityTable(obj, specificityTable, isoSpecificityTables)
        
        % -------------------------------------------------------------------------
        % Parse necessary input
        % -------------------------------------------------------------------------
        % Handle case that isoSpecificityTable is not provided
        if nargin < 3
            isoSpecificityTables = [];
        end
        if nargin < 2 || ~isa(specificityTable, 'OTTable')
            error('matlabFunctions:invalidArguments', 'The provided specificity table is not a valid OTTable');
        end
        if nargin == 3 && ~isa(isoSpecificityTables, 'OTTable')
            error('matlabFunctions:invalidArguments', 'The provided isoform specificity table array does not contain OTTables');
        end
        
        % -------------------------------------------------------------------------
        % Update table
        % -------------------------------------------------------------------------
        obj.specificityTable = specificityTable;
        obj.isoSpecificityTables = isoSpecificityTables;
        
        % -------------------------------------------------------------------------
        % Check seed lengths
        % -------------------------------------------------------------------------
        if ~isempty(isoSpecificityTables)
            isoSeedLength = unique([obj.isoSpecificityTables.seedLength]);
            if length(isoSeedLength) ~= 1
                error('matlabFunctions:invalidArguments', 'The OTTable list contains elements with different seed lengths!');
            end
        else
            isoSeedLength = unique([obj.isoSpecificityTables.seedLength]);
            if isoSeedLength ~= obj.specificityTable.seedLength
                error('matlabFunctions:invalidArguments', 'The seed lengths must be equal for the specificity table and the isoform specificity tables!');
            end
        end

        % -------------------------------------------------------------------------
        % Display progress
        % -------------------------------------------------------------------------
        if obj.verbose
            PageBreak();
            tic;
            display(['Calculating region specificity for a seed length of ' ...
                num2str(obj.specificityTable.seedLength) ': ' datestr(now)]);
        end
        
        % -------------------------------------------------------------------------
        % Build specificity score
        % -------------------------------------------------------------------------
        
        % Determine normalization based on the weighting of the table
        if obj.specificityTable.uniformWeight 
            normalization = ones(1, length(obj.transcriptome.intSequences));
        else
            normalization = obj.transcriptome.abundance;
        end
        
        % Determine if isoform specificity should be calculated
        if isempty(obj.isoSpecificityTables) % If no isoform specificity tables are provided, calculate specificity directly
            display(['Calculating specificity without isoform information...']);
            for s=1:length(obj.transcriptome.intSequences)
                obj.specificity{s} = normalization(s)./(obj.specificityTable.CalculatePenalty(obj.transcriptome.intSequences{s}));
                if obj.verbose && ~mod(s,1000)
                    display(['... completed ' num2str(s) ' seqs']);
                end
            end
        else
            display(['Utilizing isoform information']);
            for s=1:length(obj.transcriptome.intSequences)
                % Get sequence name and find OTTables
                localGeneName = obj.transcriptome.geneNames{s};
                localTable = obj.isoSpecificityTables(strcmp({obj.isoSpecificityTables.name}, localGeneName));
                
                % Confirm that a table was found
                if length(localTable) > 1
                    warning('matlab:invalidArguments', ['Found more than one OTTable for ' localGeneName]);
                    continue;
                elseif isempty(localTable)
                    warning('matlab:invalidArguments', ['Did not find an OTTable for ' localGeneName]);
                    continue;
                end
                
                % Calculate isoform adjusted specificity
                isoCounts = localTable.CalculatePenalty(obj.transcriptome.intSequences{s});
                obj.isoSpecificity{s} = normalization(s)./isoCounts;
                obj.specificity{s} = isoCounts./ ... 
                    obj.specificityTable.CalculatePenalty(obj.transcriptome.intSequences{s});
                
                % Update progress
                if obj.verbose && ~mod(s,1000)
                    display(['... completed ' num2str(s) ' seqs']);
                end
            end
        end
        
        % -------------------------------------------------------------------------
        % Display progress
        % -------------------------------------------------------------------------
        if obj.verbose
            display(['... completed in ' num2str(toc) ' at ' datestr(now)]);
        end
        
    end
    
    % -------------------------------------------------------------------------
    % SetParallel
    % -------------------------------------------------------------------------
    function SetParallel(obj, p)
        % Set or update the parallel.Pool object
        % obj.setParallel(p) 
        % obj.setParallel([]) removes the existing pool
        
        % -------------------------------------------------------------------------
        % Check validity
        % -------------------------------------------------------------------------
        if ~isempty(p) && ~isa(p, 'parallel.Pool')
            error('matlabFunctions:invalidArgument', 'Must provide a valid parallel.Pool object');
        end
        
        obj.parallel = p;
        if isempty(obj.parallel)
            obj.numPar = 0;
        else
            obj.numPar = obj.parallel.NumWorkers;
        end
    end
    
    % -------------------------------------------------------------------------
    % Get the forbidden seq penalty for all putative regions
    % -------------------------------------------------------------------------
    function [noForbiddenSeqs, ids, names] = GetRegionForbiddenSeqs(obj, len, varargin)
        % Return a boolean for each valid probe
        % [noForbiddenSeqs, ids, names] = obj.GetRegionForbiddenSeqs(len)
        % [noForbiddenSeqs, ids, names] = obj.GetRegionForbiddenSeqs(len, 'geneName', names)
        % [noForbiddenSeqs, ids, names] = obj.GetRegionForbiddenSeqs(len, 'geneID', ids)
        % noForbiddenSeqs specifies if the target region starting at each position
        % does not contain any of the forbidden sequences

        % -------------------------------------------------------------------------
        % Get internal indices for the requested transcripts
        % -------------------------------------------------------------------------
        [inds, ids, names] = obj.transcriptome.GetInternalInds(varargin{:});
        
        % -------------------------------------------------------------------------
        % Display Progress
        % -------------------------------------------------------------------------
        if obj.verbose
            PageBreak();
            timer = tic;
            display(['Finding forbidden sequences for ' num2str(length(inds)) ' transcripts']);
        end
        
        % -------------------------------------------------------------------------
        % Search for forbidden sequences
        % -------------------------------------------------------------------------
        for s=1:length(obj.forbiddenSeqs) % Loop over forbidden sequences
            % Display progress
            if obj.verbose
                display(['Finding ' int2nt(obj.forbiddenSeqs{s}+1)]);
            end
            
            % Calculate filter window
            windowLen = len - length(obj.forbiddenSeqs{s}) + 1;
            
            % Handle first forbidden sequence
            if s==1
                for l=1:length(inds) % Loop over transcriptome
                    hasSeq = filter(ones(1, windowLen)/windowLen, 1, obj.isForbiddenSeq{s}{inds(l)});
                    noForbiddenSeqs{l} = ~hasSeq(windowLen:end);
                end
            else
                for l=1:length(inds) % Loop over transcriptome
                    hasSeq = filter(ones(1, windowLen)/windowLen, 1, obj.isForbiddenSeq{s}{inds(l)});
                    noForbiddenSeqs{l} = ~hasSeq(windowLen:end) & noForbiddenSeqs{l};
                end
            end
        end
        
        % -------------------------------------------------------------------------
        % Display Progress
        % -------------------------------------------------------------------------
        if obj.verbose
            display(['... completed in ' num2str(toc(timer))]);
        end
    end

    % -------------------------------------------------------------------------
    % Get all valid regions
    % -------------------------------------------------------------------------
    function [isValid, ids, names] = GetRegionValidity(obj, len, varargin)
        % Return a boolean for each valid probe
        % [isValid, ids, names] = obj.GetRegionValidity(len)
        % [isValid, ids, names] = obj.GetRegionValidity(len, 'geneName', names)
        % [isValid, ids, names] = obj.GetRegionValidity(len, 'geneID', ids)
        % isValid specifies if the target region starting at each position
        % is valid

        % -------------------------------------------------------------------------
        % Get internal indices for the requested transcripts
        % -------------------------------------------------------------------------
        [inds, ids, names] = obj.transcriptome.GetInternalInds(varargin{:});
        
        % -------------------------------------------------------------------------
        % Compute sliding GC window
        % -------------------------------------------------------------------------
        if obj.verbose
            PageBreak();
            timer = tic;
            display(['Calculating validity for ' num2str(length(inds)) ' transcripts']);
        end
        for l=1:length(inds)
            valid = filter(ones(1, len)/len, 1, ...
                obj.transcriptome.intSequences{inds(l)} > 3 | ...
                obj.transcriptome.intSequences{inds(l)} < 0);
            isValid{l} = ~valid(len:end);
        end
        if obj.verbose
            display(['... completed in ' num2str(toc(timer))]);
        end
    end
    
    % -------------------------------------------------------------------------
    % Get GC values for all putative regions
    % -------------------------------------------------------------------------
    function [GC, ids, names] = GetRegionGC(obj, len, varargin)
        % Return the GC content of a set of probes of specified length
        % [GC, ids, names, validRequest] = obj.GetRegionGC(len)
        % [GC, ids, names, validRequest] = obj.GetRegionGC(len, 'geneName', names)
        % [GC, ids, names, validRequest] = obj.GetRegionGC(len, 'geneID', ids)
        % GC is a cell array of all GC values 
        % ids are the gene ids associated with each GC curve
        % names are the gene names associated with each GC curve
        % validRequest is a boolean that can be used to confirm that a
        % requested gene/id is actually in the transcriptome. 
        
        % -------------------------------------------------------------------------
        % Get internal indices for the requested transcripts
        % -------------------------------------------------------------------------
        [inds, ids, names] = obj.transcriptome.GetInternalInds(varargin{:});
        
        % -------------------------------------------------------------------------
        % Compute sliding GC window
        % -------------------------------------------------------------------------
        if obj.verbose
            PageBreak();
            timer = tic;
            display(['Calculating GC content for ' num2str(length(inds)) ' transcripts']);
        end
        for l=1:length(inds)
            gc = filter(ones(1, len)/len, 1, obj.gc{inds(l)});
            GC{l} = gc(len:end);
        end
        if obj.verbose
            display(['... completed in ' num2str(toc(timer))]);
        end
    end
           
    % -------------------------------------------------------------------------
    % GetRegionPenalty
    % -------------------------------------------------------------------------
    function [penalty, ids, names] = GetRegionPenalty(obj, len, OTtableName, varargin)
        % Return the penalty associated with all possible probes of specified length
        % [penalty, ids, names] = obj.GetRegionPenalty(len, OTtableName)
        % [penalty, ids, names] = obj.GetRegionPenalty(len, OTtableName, 'geneName', names)
        % [penalty, ids, names] = obj.GetRegionPenalty(len, OTtableName, 'geneID', ids)
        % penalty is a cell array of all penalties 
        % ids are the gene ids associated with each curve
        % names are the gene names associated with each curve
        
        % -------------------------------------------------------------------------
        % Get internal indices for the requested transcripts
        % -------------------------------------------------------------------------
        [inds, ids, names] = obj.transcriptome.GetInternalInds(varargin{:});
        
        % -------------------------------------------------------------------------
        % Find table ind
        % -------------------------------------------------------------------------
        id = find(strcmp(obj.OTTableNames, OTtableName));
        if isempty(id)
            error('matlabFunctions:invalidArgument', 'The specified table does not exist');
        end
        if length(id) > 1
            warning('Multiple tables matching the specified name were found');
        end
        
        localPenalties = obj.penalties{id};
        
        % -------------------------------------------------------------------------
        % Compute appropriate window length
        % -------------------------------------------------------------------------
        windowLen = len - obj.OTTables(id).seedLength + 1;
        
        % -------------------------------------------------------------------------
        % Compute sliding GC window
        % -------------------------------------------------------------------------
        if obj.verbose
            PageBreak();
            timer = tic;
            display(['Calculating penalty for ' num2str(length(inds)) ' transcripts' ...
                ' using table: ' OTtableName]);
        end
        for l=1:length(inds)
            pen = filter(ones(1, windowLen), 1, localPenalties{inds(l)}); % Return sum, not average
            penalty{l} = pen(windowLen:end);
        end
        if obj.verbose
            display(['... completed in ' num2str(toc(timer))]);
        end
    end
    
    % -------------------------------------------------------------------------
    % Get the specificity of putative target regions
    % -------------------------------------------------------------------------
    function [specificity, ids, names] = GetRegionSpecificity(obj, len, varargin)
        % Return the penalty associated with all possible probes of specified length
        % [specificity, ids, names] = obj.GetRegionSpecificity(len)
        % ... = obj.GetRegionSpecificity(len, 'geneName', names)
        % ... = obj.GetRegionSpecificity(len, 'geneID', ids)
        
        % -------------------------------------------------------------------------
        % Get internal indices for the requested transcripts
        % -------------------------------------------------------------------------
        [inds, ids, names] = obj.transcriptome.GetInternalInds(varargin{:});
        
        % -------------------------------------------------------------------------
        % Compute appropriate window length
        % -------------------------------------------------------------------------
        windowLen = len - obj.specificityTable.seedLength + 1;
        
        % -------------------------------------------------------------------------
        % Compute sliding window
        % -------------------------------------------------------------------------
        if obj.verbose
            PageBreak();
            timer = tic;
            display(['Calculating specificity for ' num2str(length(inds)) ' transcripts']);
        end
        for l=1:length(inds)
            spec = filter(ones(1, windowLen)/windowLen, 1, obj.specificity{inds(l)}); % Return average
            specificity{l} = spec(windowLen:end);
        end
        if obj.verbose
            display(['... completed in ' num2str(toc(timer))]);
        end
    end
    
    % -------------------------------------------------------------------------
    % Get the isoform specificity of putative target regions
    % -------------------------------------------------------------------------
    function [isoSpecificity, ids, names] = GetRegionIsoSpecificity(obj, len, varargin)
        % Return the penalty associated with all possible probes of specified length
        % [specificity, ids, names] = obj.GetRegionIsoSpecificity(len)
        % ... = obj.GetRegionIsoSpecificity(len, 'geneName', names)
        % ... = obj.GetRegionIsoSpecificity(len, 'geneID', ids)
        
        % -------------------------------------------------------------------------
        % Get internal indices for the requested transcripts
        % -------------------------------------------------------------------------
        [inds, ids, names] = obj.transcriptome.GetInternalInds(varargin{:});
        
        % -------------------------------------------------------------------------
        % Compute appropriate window length
        % -------------------------------------------------------------------------
        windowLen = len - obj.isoSpecificityTables(1).seedLength + 1;
        
        % -------------------------------------------------------------------------
        % Compute sliding window
        % -------------------------------------------------------------------------
        if obj.verbose
            PageBreak();
            timer = tic;
            display(['Calculating isoform specificity for ' num2str(length(inds)) ' transcripts']);
        end
        for l=1:length(inds)
            spec = filter(ones(1, windowLen)/windowLen, 1, obj.isoSpecificity{inds(l)}); % Return average
            isoSpecificity{l} = spec(windowLen:end);
        end
        if obj.verbose
            display(['... completed in ' num2str(toc(timer))]);
        end
    end
    
    % -------------------------------------------------------------------------
    % Return the Tm for all regions of the desired transcripts
    % -------------------------------------------------------------------------
    function [Tm, ids, names] = GetRegionTm(obj, len, varargin)
        % Return the Tm for all specified probes of length, len
        % [Tm, ids, names, validRequest] = obj.GetRegionTm(len)
        % [Tm, ids, names, validRequest] = obj.GetRegionTm(len, 'geneName', names)
        % [Tm, ids, names, validRequest] = obj.GetRegionTm(len, 'geneID', ids)
        % [Tm, ids, names, validRequest] = obj.GetRegionTm(len, 'monovalentSalt', saltConc)
        % [Tm, ids, names, validRequest] = obj.GetRegionTm(len, 'probeConc', probeConc)
        % Tm is a cell array of all Tm values 
        
        % -------------------------------------------------------------------------
        % Parse variable inputs
        % -------------------------------------------------------------------------
        defaults = cell(0,3); 
        defaults(end+1,:) = {'geneName', 'cell', {}}; % Display progress of construction
        defaults(end+1,:) = {'geneID', 'cell', {}};
        defaults(end+1,:) = {'monovalentSalt', 'positive', 0.3};
        defaults(end+1,:) = {'probeConc', 'positive', 5e-9};
        defaults(end+1,:) = {'inds', 'positive',[]};
        
        parameters = ParseVariableArguments(varargin, defaults, mfilename);
        
        % -------------------------------------------------------------------------
        % Get internal indices for the requested transcripts
        % -------------------------------------------------------------------------
        [inds, ids, names] = obj.transcriptome.GetInternalInds('parameters', parameters);
        
        % -------------------------------------------------------------------------
        % Display Progress
        % -------------------------------------------------------------------------
        if obj.verbose
            PageBreak();
            timer = tic;
            display(['Calculating Tm for ' num2str(length(inds)) ' transcripts' ...
                ' with ' num2str(parameters.monovalentSalt) ' M salt and ' ...
                num2str(parameters.probeConc) ' M probe']);
        end
        
        % -------------------------------------------------------------------------
        % Calculate Tm
        % -------------------------------------------------------------------------
        localDG = obj.dG(inds); % Local thermodynamic properties
        seqs = obj.transcriptome.intSequences(inds); % Local copy of sequences
        for l=1:length(inds) % Loop over requested sequences
            % Get local sequence
            intSeq = seqs{l};
            
            % Calculate total H and S per putative probe
            H = filter(ones(1, len-1), 1, localDG{l}(1,:)); % len-1 nn per sequence of len
            S = filter(ones(1, len-1), 1, localDG{l}(2,:));
            H = H((len-1):end);
            S = S((len-1):end);
            
            % Determine ends
            fivePrimeAT = intSeq(1:(end-len+1)) == 0 | intSeq(1:(end-len+1)) == 3;
            threePrimeAT = intSeq(len:end) == 0 | intSeq(len:end) == 3;
            
            % Add end corrections
            H = H + 0.2 + 2.2*fivePrimeAT + 2.2*threePrimeAT;
            S = S + -5.7 + 6.9*fivePrimeAT + 6.9*threePrimeAT;
            
            % Apply salt correction
            S = S + 0.368*(len-1)*log(parameters.monovalentSalt);
            
            % Calculate Tm in C
            Tm{l} = H*1000 ./ (S + 1.9872 * log(parameters.probeConc)) - 273.15; 
            
            % NOTE: For the future, I should provide two concentrations,
            % probe and target, and this should be log(probeC - targetC/2)
            % where probeC > targetC. They are switched in the
            % concentrations are switched.
        end
        
        % -------------------------------------------------------------------------
        % Display Progress
        % -------------------------------------------------------------------------
        if obj.verbose
            display(['... completed in ' num2str(toc(timer))]);
        end
    end
    
    % -------------------------------------------------------------------------
    % Design target regions
    % -------------------------------------------------------------------------
    function targetRegions = DesignTargetRegions(obj,varargin)
        % This method returns target regions, tiled to be non-overlapping,
        % and subject to a variety of constraints.  
        % targetRegions = obj.DesignTargetRegions(..., 'regionLength', [allPossibleLengthValues])
        % targetRegions = obj.DesignTargetRegions(..., 'Tm', [low,up])
        % targetRegions = obj.DesignTargetRegions(..., 'GC', [low,up])
        % targetRegions = obj.DesignTargetRegions(..., 'OTTables', {'name', [low, up], 'name', range, ...})
        % targetRegions = obj.DesignTargetRegions(..., 'geneName', names)
        % targetRegions = obj.DesignTargetRegions(..., 'geneID', ids)
        % targetRegions = obj.DesignTargetRegions(..., 'specificity', [low, up])
        % targetRegions = obj.DesignTargetRegions(..., 'isoSpecificity', [low, up])

        % -------------------------------------------------------------------------
        % Parse variable inputs
        % -------------------------------------------------------------------------
        % Define defaults
        defaults = cell(0,3); 
        defaults(end+1,:) = {'geneName', 'cell', {}}; % Display progress of construction
        defaults(end+1,:) = {'geneID', 'cell', {}};
        defaults(end+1,:) = {'regionLength', 'positive', 30};
        defaults(end+1,:) = {'Tm', 'positive', []};
        defaults(end+1,:) = {'GC', 'positive', []};
        defaults(end+1,:) = {'specificity', 'positive', []};
        defaults(end+1,:) = {'isoSpecificity', 'positive', []}; 
        defaults(end+1,:) = {'monovalentSalt', 'positive', 0.3};
        defaults(end+1,:) = {'probeConc', 'positive', 5e-9};
        defaults(end+1,:) = {'OTTables', 'cell', {}};
        defaults(end+1,:) = {'includeSequence', 'boolean', true};
        defaults(end+1,:) = {'threePrimeSpace', 'integer', 0};
        defaults(end+1,:) = {'removeForbiddenSeqs', 'boolean', false};
        
        % Parse variable inputs
        parameters = ParseVariableArguments(varargin, defaults, mfilename);
                
        % -------------------------------------------------------------------------
        % Get internal indices for the requested transcripts
        % -------------------------------------------------------------------------
        inds = obj.transcriptome.GetInternalInds('parameters', parameters);
        
        % -------------------------------------------------------------------------
        % Handle case of no valid requested transcripts
        % -------------------------------------------------------------------------
        if isempty(inds)
            warning('matlabFunctions:noValidEntries', ...
                'None of the requested entries were found');
            targetRegions = TargetRegions();
            return;
        end
        
        % -------------------------------------------------------------------------
        % Initialize variables for loops
        % -------------------------------------------------------------------------
        targetRegions = [];
        indsToKeep = cell(length(parameters.regionLength), length(inds));
        
        % -------------------------------------------------------------------------
        % Loop over probe lengths
        % -------------------------------------------------------------------------
        for l=1:length(parameters.regionLength)
            % -------------------------------------------------------------------------
            % Display progress
            % -------------------------------------------------------------------------
            if obj.verbose
                PageBreak();
                timer1= tic;
                display(['Designing target regions of length ' num2str(parameters.regionLength(l))]);
            end
            
            % -------------------------------------------------------------------------
            % Calculate all region properties
            % -------------------------------------------------------------------------
            GC(l,:) = obj.GetRegionGC(parameters.regionLength(l), 'inds', inds);
            Tm(l,:) = obj.GetRegionTm(parameters.regionLength(l), 'inds', inds);
            for p=1:length(obj.OTTables)
                penalties(l,p,:) = obj.GetRegionPenalty(parameters.regionLength(l), ...
                    obj.OTTableNames{p}, 'inds', inds);
            end
            specificity(l,:) = obj.GetRegionSpecificity(parameters.regionLength(l), 'inds', inds);
            isoSpecificity(l,:) = obj.GetRegionIsoSpecificity(parameters.regionLength(l), 'inds', inds);
            
            % -------------------------------------------------------------------------
            % Define numerical pad: to address round off error in filter
            % -------------------------------------------------------------------------
            numPad = 10*eps;
            
            % -------------------------------------------------------------------------
            % Initialize Inds to Keep
            % -------------------------------------------------------------------------
            indsToKeep(l,:) = obj.GetRegionValidity(parameters.regionLength(l), 'inds', inds);
            
            % -------------------------------------------------------------------------
            % Determine probe properties (only if needed) and cut
            % -------------------------------------------------------------------------
            if parameters.removeForbiddenSeqs
                noForbiddenSeqs = obj.GetRegionForbiddenSeqs(parameters.regionLength(l), 'inds', inds);
                for i=1:length(noForbiddenSeqs)
                    indsToKeep{l,i} = indsToKeep{l,i} & noForbiddenSeqs{i};
                end
            end
            if ~isempty(parameters.GC)
                for i=1:length(GC)
                    indsToKeep{l,i} = indsToKeep{l,i} & GC{l,i} >= (parameters.GC(1)-numPad) & GC{l,i} <= (parameters.GC(2)+numPad);
                end
            end
            if ~isempty(parameters.Tm)
                for i=1:length(Tm)
                    indsToKeep{l,i} = indsToKeep{l,i} & Tm{l,i} >= (parameters.Tm(1)-numPad) & Tm{l,i} <= (parameters.Tm(2)+numPad);
                end
            end
            if ~isempty(parameters.specificity)
                for i=1:length(specificity)
                    indsToKeep{l,i} = indsToKeep{l,i} & specificity{l,i} >= (parameters.specificity(1)-numPad) ...
                        & specificity{l,i} <= (parameters.specificity(2)+numPad);
                end
            end
            if ~isempty(parameters.isoSpecificity)
                for i=1:length(isoSpecificity)
                    indsToKeep{l,i} = indsToKeep{l,i} & isoSpecificity{l,i} >= (parameters.isoSpecificity(1)-numPad) ...
                        & isoSpecificity{l,i} <= (parameters.isoSpecificity(2)+numPad);
                end
            end
            if ~isempty(parameters.OTTables)
                for t=1:(length(parameters.OTTables)/2)
                    pid = find(strcmp(obj.OTTableNames, parameters.OTTables{2*(t-1)+1}), 1);
                    if isempty(pid)
                        error('matlabFunctions:invalidArguments', 'Unrecognized OTTable name');
                    end
                    for i=1:length(inds)
                        indsToKeep{l,i} = indsToKeep{l,i} & penalties{l,pid,i} >= (parameters.OTTables{2*t}(1)-numPad) ...
                            & penalties{l,pid,i} <= (parameters.OTTables{2*t}(2)+numPad);
                    end
                end
            end
            
            % -------------------------------------------------------------------------
            % Display progress
            % -------------------------------------------------------------------------
            if obj.verbose
                display(['... completed identification of possible targets of length ' ...
                    num2str(parameters.regionLength(l)) ' nt in ' num2str(toc(timer1)) ' s']);
            end
        end
        
        % -------------------------------------------------------------------------
        % Display Progress
        % -------------------------------------------------------------------------
        % Display progress
        if obj.verbose
            PageBreak();
            timer2 = tic;
            display(['Tiling and compiling target regions: ' datestr(now)]);
        end
          
        % -------------------------------------------------------------------------
        % Prepare variables for parfor
        % -------------------------------------------------------------------------
        targetRegions = cell(1, length(inds));
        regionLength = parameters.regionLength;
        threePrimeSpace = parameters.threePrimeSpace;
        ids = obj.transcriptome.ids(inds);
        geneNames = obj.transcriptome.geneNames(inds);
        intSequences = obj.transcriptome.intSequences(inds);
        
        % -------------------------------------------------------------------------
        % Loop over objects in parallel
        % -------------------------------------------------------------------------
        parfor (i=1:length(inds), obj.numPar)   
            % Compile region properties
            regionProps = zeros(6,0);
            for l=1:length(regionLength)
                startPos = find(indsToKeep{l,i});
                localProps = zeros(6,length(startPos));
                if ~isempty(startPos) % Only build object if there are any valid starting positions
                    localProps(1,:) = startPos;
                    localProps(2,:) = regionLength(l)*ones(1, sum(indsToKeep{l,i}));
                    localProps(3,:) = Tm{l,i}(indsToKeep{l,i});
                    localProps(4,:) = GC{l,i}(indsToKeep{l,i});
                    localProps(5,:) = specificity{l,i}(indsToKeep{l,i});
                    localProps(6,:) = isoSpecificity{l,i}(indsToKeep{l,i});
                    regionProps = [regionProps localProps];
                end
            end
            
            % Tile regions and return properties of selected regions
            selectedRegionProps = TRDesigner.TileRegions(regionProps, ...
                threePrimeSpace);
            
            % Build a new target region object
            newTR = TargetRegions('id', ids{i}, ...
                'geneName', geneNames{i}, ...
                'geneSequence', intSequences{i}, ...
                'startPos', selectedRegionProps(1,:), ...
                'regionLength', selectedRegionProps(2,:), ...
                'Tm', selectedRegionProps(3,:), ...
                'GC', selectedRegionProps(4,:), ...
                'specificity',selectedRegionProps(5,:), ...
                'isoSpecificity', selectedRegionProps(6,:));
            
            % Append to growing list of target regions
            targetRegions{i} = newTR;
        end
        
        % -------------------------------------------------------------------------
        % Flatten objects
        % -------------------------------------------------------------------------
        targetRegions = [targetRegions{:}];
        
        % -------------------------------------------------------------------------
        % Display progress
        % -------------------------------------------------------------------------
        if obj.verbose
            display(['... completed in ' num2str(toc(timer2)) ' s']);
        end        
    end
        
    % -------------------------------------------------------------------------
    % Save Function
    % -------------------------------------------------------------------------
    function Save(obj, dirPath)
        % Save the TRDesigner object in a directory specified by dirPath
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
        % Define fields to save
        % -------------------------------------------------------------------------
        fieldsToSave = properties(obj);
        fieldsToSave = setdiff(fieldsToSave, {'parallel', 'numPar'});
        
        % -------------------------------------------------------------------------
        % Save fields
        % -------------------------------------------------------------------------
        for i=1:length(fieldsToSave)
            switch fieldsToSave{i}
                case 'OTTables'
                    for j=1:length(obj.OTTables)
                        table = obj.OTTables(j);
                        table.Save([dirPath 'OTTable_' num2str(j)]);
                    end
                case 'specificityTable'
                    obj.specificityTable.Save([dirPath 'specificityTable']);
                case 'isoSpecificityTables'
                    obj.isoSpecificityTables.Save([dirPath 'isoSpecificityTables']);
                case 'transcriptome'
                    obj.transcriptome.Save([dirPath 'transcriptome']);
                otherwise
                    SaveAsByteStream([dirPath fieldsToSave{i} '.matb'], ...
                        obj.(fieldsToSave{i}), 'verbose', obj.verbose);
            end
        end
    end
end

% -------------------------------------------------------------------------
% Static methods
% -------------------------------------------------------------------------
methods (Static)
    % -------------------------------------------------------------------------
    % Build a TRDesigner object from a saved version
    % -------------------------------------------------------------------------
    function obj = Load(dirPath)
        % obj = TRDesigner.Load(dirPath)
        
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
        obj = TRDesigner();
        
        % -------------------------------------------------------------------------
        % Define fields to load
        % -------------------------------------------------------------------------
        fieldsToLoad = setdiff(properties(obj), {'parallel', 'numPar'});
        
        % -------------------------------------------------------------------------
        % Load properties/data
        % -------------------------------------------------------------------------
        for i=1:length(fieldsToLoad)
            switch fieldsToLoad{i}
                case 'OTTables'
                    contents = dir([dirPath 'OTTable_*']);
                    contents = contents([contents.isdir]);
                    for j=1:length(contents)
                        try
                            obj.OTTables(j) = OTTable.Load([dirPath contents(j).name]);
                        catch
                            obj.OTTables = OTTable.Load([dirPath contents(j).name]);
                        end
                    end
                case 'specificityTable'
                    obj.specificityTable = OTTable.Load([dirPath 'specificityTable']);
                case 'isoSpecificityTables'
                    obj.isoSpecificityTables = OTTable.Load([dirPath 'isoSpecificityTables']);
                case 'transcriptome'
                    obj.transcriptome = Transcriptome.Load([dirPath 'transcriptome']);
                otherwise
                    obj.(fieldsToLoad{i}) = LoadByteStream([dirPath fieldsToLoad{i} '.matb'], ...
                    'verbose', true);
            end
        end
    end
    
    % -------------------------------------------------------------------------
    % SantaLucia Nearest Neighor Calculations  
    % -------------------------------------------------------------------------
    function dG = SantaLuciaNearestNeighbor(intSeq)
        % Return the enthalpy and entropy associated with each nearest
        % neighbor pair in the sequence
        % dG = TRDesigner.SantaLuciaNearestNeighbor(intSeq)
        
        % -------------------------------------------------------------------------
        % Coerce sequence format  
        % -------------------------------------------------------------------------
        if ischar(intSeq)
            intSeq = int8(nt2int(intSeq)) - 1;
        end
        
        % -------------------------------------------------------------------------
        % Initialize dG: dG(1,:) is the enthalpy (kcal/mol); dG(2,:) is the entropy (cal/(mol K) 
        % -------------------------------------------------------------------------
        dG = nan(2, length(intSeq)-1); % nan is a not valid flag
        
        % -------------------------------------------------------------------------
        % Assign terms: AA = 1, AC = 2, ... CA = 5 ... TT = 16
        % -------------------------------------------------------------------------
        nnID = 4*intSeq(1:(end-1)) + intSeq(2:end) + 1; % Convert pairs to index
        % AC -> 4*0 + 1 + 1= 2; CA = 4*1 + 0 + 1 = 5
        isValid = intSeq(1:(end-1)) <= 3 & intSeq(2:end) <= 3 & ...
            intSeq(1:(end-1)) >= 0 & intSeq(2:end) >=0 ; % Check for ambiguious nucleotides

        % -------------------------------------------------------------------------
        % Nearest neighbor H and S from Santa Lucia Jr and Hicks, Annu.
        % Rev. Biomol. Struct. 2004 33:415-40. Units of kcal/mol (H) or
        % cal/(mol K) (S).
        % -------------------------------------------------------------------------
        H = [-7.6 -8.4 -7.8 -7.2 -8.5 -8.0 -10.6 -7.8 ...   % AA/TT GT/CA CT/GA AT/TA CA/GT GG/CC CG/GC CT/GA
            -8.2 -9.8 -8.0 -8.4 -7.2 -8.2 -8.5 -7.6];       % GA/CT GC/CG GG/CC GT/CA TA/AT GA/CT CA/GT AA/TT
        S = [-21.3 -22.4 -21.0 -20.4 -22.7 -19.9 -27.2 -21.0 ...    % AA/TT GT/CA CT/GA AT/TA CA/GT GG/CC CG/GC CT/GA
            -22.2 -24.4 -19.9 -22.4 -21.3 -22.2 -22.7 -21.3];       % GA/CT GC/CG GG/CC GT/CA TA/AT GA/CT CA/GT AA/TT
        % Note the different values for AA/TT from Santa Lucia Jr PNAS 1998

        % -------------------------------------------------------------------------
        % Nearest neighbour H and S from Santa Lucia Jr, PNAS 95, 1460-65
        % (1998)
        % -------------------------------------------------------------------------
%         H = [-7.9 -8.4 -7.8 -7.2 -8.5 -8.0 -10.6 -7.8 ...   % AA/TT GT/CA CT/GA AT/TA CA/GT GG/CC CG/GC CT/GA
%             -8.2 -9.8 -8.0 -8.4 -7.2 -8.2 -8.5 -7.6];       % GA/CT GC/CG GG/CC GT/CA TA/AT GA/CT CA/GT AA/TT
%         S = [-22.2 -22.4 -21.0 -20.4 -22.7 -19.9 -27.2 -21.0 ...    % AA/TT GT/CA CT/GA AT/TA CA/GT GG/CC CG/GC CT/GA
%             -22.2 -24.4 -19.9 -22.4 -21.3 -22.2 -22.7 -21.3];       % GA/CT GC/CG GG/CC GT/CA TA/AT GA/CT CA/GT AA/TT
        
        % -------------------------------------------------------------------------
        % Define dG
        % -------------------------------------------------------------------------
        dG(1,isValid) = H(nnID(isValid));
        dG(2,isValid) = S(nnID(isValid));      
    end
    
    % -------------------------------------------------------------------------
    % Tile Regions  
    % -------------------------------------------------------------------------
    function selectedRegionData = TileRegions(regionProps, padLength)
        % Identify a non-overlapping tiling of regions separated by at
        % least padLength
        % selectedRegionData = TRDesigner.TileRegions(regionProps, padLength)
        
        % -------------------------------------------------------------------------
        % Handle empty
        % -------------------------------------------------------------------------
        if isempty(regionProps)
            selectedRegionData = regionProps;
            return;
        end
                  
        
        % -------------------------------------------------------------------------
        % Find start positions, sort, and find next available positions
        % -------------------------------------------------------------------------
        startPos = regionProps(1,:);
        [startPos, sind] = sort(startPos);
        
        regionProps = regionProps(:,sind); % Sort data
        nextAvailablePos = startPos + regionProps(2,:) + padLength;
        
        % -------------------------------------------------------------------------
        % Tile probes
        % -------------------------------------------------------------------------
        done = false;
        indsToKeep = 1;
        while ~done
            minNextPos = nextAvailablePos(indsToKeep(end)); % Identify the minimum starting position
            newInd = find(startPos >= minNextPos, 1);
            if isempty(newInd)
                done = true;
            else
                indsToKeep(end+1) = newInd;
            end
        end
        
        % -------------------------------------------------------------------------
        % Return probe data
        % -------------------------------------------------------------------------
        selectedRegionData = regionProps(:,indsToKeep);
    end
        
    
end % Static methods

end
