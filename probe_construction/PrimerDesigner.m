classdef PrimerDesigner < handle
% ------------------------------------------------------------------------
% [pDesignerObj] = PrimerDesigner(varargin)
% This class designs orthogonal primers. 
%--------------------------------------------------------------------------
% Necessary Inputs
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% May 4, 2015
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Define properties
% -------------------------------------------------------------------------
properties 
    verbose             % Determines the verbosity of the class
    ntComposition       % A 4x1 vector that controls the composition of the primers (A,C,G,T)
    OTTables            % OTTables an array of Off-Target Tables for calculating penalties
    OTTableNames        % A cell array of string names for each Off-Target table
    monovalentSalt      % The monovalent salt concentration (M) for Tm calculations
    primerConc          % The concentration of primer (M) for Tm calculations
    seqsToRemove        % A cell array of sequences that are not allowed, e.g. GGGG
end
properties (SetAccess=protected)
    numPrimers          % The number of current primers in the class
    primerLength        % The length of the primers
    seqs                % The sequences of the primers in integer format, a NxL matrix 
    parallel            % A parallel.Pool object
    numPar              % The number of current parallel workers
    homologyMax         % The maximum homology length accepted
    gc                  % The GC content for each sequence
    Tm                  % The Tm (C) for each sequence
    penalties           % The penalty values associated with each sequence for each OTTable
end

properties (GetAccess=private)
end

properties (Hidden=true)
    homologyMat         % A sparse matrix specifying the primer pairs that share homology
    seqHash             % The hash values for all sequences
    seqRCHash           % The hash values for the reverse complement of all sequences
end

% -------------------------------------------------------------------------
% Public Methods
% -------------------------------------------------------------------------
methods
    
    % -------------------------------------------------------------------------
    % Define constructor
    % -------------------------------------------------------------------------
    function obj = PrimerDesigner(varargin)
        % Create the PrimerDesigner object

        % -------------------------------------------------------------------------
        % Parse variable inputs
        % -------------------------------------------------------------------------
        defaults = cell(0,3); 
        
        defaults(end+1,:) = {'verbose', 'boolean', true}; % Display progress of construction
        defaults(end+1,:) = {'ntComposition', 'positive', [0.25 0.25 0.25 0.25]};
        defaults(end+1,:) = {'OTTables', 'freeType', []};
        defaults(end+1,:) = {'OTTableNames', 'cell', {}};
        defaults(end+1,:) = {'parallel', 'parallel', []};
        defaults(end+1,:) = {'seqs', 'array', []};
        defaults(end+1,:) = {'primerLength', 'positive', 20};
        defaults(end+1,:) = {'numPrimersToGenerate', 'positive', 1e6};
        defaults(end+1,:) = {'homologyMax', 'positive', 8};
        defaults(end+1,:) = {'monovalentSalt', 'positive', 0.3};
        defaults(end+1,:) = {'primerConc', 'positive', 0.5e-6};
        defaults(end+1,:) = {'seqsToRemove', 'cell', {'AAAA', 'TTTT', 'GGGG', 'CCCC'}};
        
        parameters = ParseVariableArguments(varargin, defaults, mfilename);
        
        % Transfer values to object
        fieldsToTransfer = setdiff(fields(parameters), {'parallel', ...
            'numPrimersToGenerate', 'primerLength'});
        for f=1:length(fieldsToTransfer)
            obj.(fieldsToTransfer{f}) = parameters.(fieldsToTransfer{f});
        end
        SetParallel(obj, parameters.parallel); %Handles transfer of parallel.Pool and sets obj.numPar
        
        % -------------------------------------------------------------------------
        % Handle empty class request
        % -------------------------------------------------------------------------
        if nargin < 1
            return;
        end
        
        % -------------------------------------------------------------------------
        % Generate random seed for randseq
        % -------------------------------------------------------------------------
        rng('shuffle');
        
        % -------------------------------------------------------------------------
        % Check validity of input sequences if provided
        % -------------------------------------------------------------------------
        if ~isempty(obj.seqs)
            if ~isa(obj.seqs, 'int8') || any(obj.seqs(:) > 3 | obj.seqs(:) < 0)
                error('matlabFunctions:invalidArgument', 'The sequences must be provided as a NxL array of unit8 from 0 to 3');
            end
            obj.primerLength = size(obj.seqs,2);
            obj.numPrimers = size(obj.seqs,1);
        else % Generate randome sequences
            obj.AddRandomSequences('parameters', parameters)
        end
            
    end
    
    % -------------------------------------------------------------------------
    % Generate Random Sequences 
    % -------------------------------------------------------------------------
    function AddRandomSequences(obj, varargin)
		% Add random sequences to the class
	
        % -------------------------------------------------------------------------
        % Handle variable input
        % -------------------------------------------------------------------------
        defaults = cell(0,3);
        defaults(end+1,:) = {'ntComposition', 'positive', [0.25 0.25 0.25 0.25]};
        defaults(end+1,:) = {'primerLength', 'positive', 20};
        defaults(end+1,:) = {'numPrimersToGenerate', 'positive', 1e6};
        
        parameters = ParseVariableArguments(varargin, defaults, mfilename);
        
        % -------------------------------------------------------------------------
        % Parse nucleotide composition 
        % -------------------------------------------------------------------------
        if ~isempty(parameters.ntComposition)
            if length(parameters.ntComposition) ~= 4
                error('matlabFunctions:invalidArguments', 'nt composition must have four entries');
            end
            parameters.ntComposition = parameters.ntComposition/sum(parameters.ntComposition);
        else
            parameters.ntComposition = 0.25*ones(1,4);
        end
        
        % -------------------------------------------------------------------------
        % Display Progress  
        % -------------------------------------------------------------------------
        if obj.verbose
            PageBreak();
            display(['Creating ' num2str(parameters.numPrimersToGenerate) ' new primers of ' ...
                num2str(parameters.primerLength) '-nt length']);
        end
        
        % -------------------------------------------------------------------------
        % Create random sequences: Following randseq
        % -------------------------------------------------------------------------
        rseq = rand(parameters.numPrimersToGenerate,parameters.primerLength);
        edges = [0, cumsum(parameters.ntComposition)];
        edges(end) = 1;
        
        [~, seq] = histc(rseq, edges);
        seq = int8(seq-1); %Shift 
        
        % -------------------------------------------------------------------------
        % Update sequences
        % -------------------------------------------------------------------------
        obj.seqs = cat(1,obj.seqs,seq);
        obj.primerLength = size(obj.seqs,2);
        obj.numPrimers = size(obj.seqs,1);
        
        % -------------------------------------------------------------------------
        % Recalculate Tm and GC and penalties
        % -------------------------------------------------------------------------
        obj.CalculatePrimerProperties();
        
    end
    
    % -------------------------------------------------------------------------
    % AddPrimer  
    % -------------------------------------------------------------------------
    function AddPrimer(obj, seq)
        % Add a specific primer sequence to the list of sequences
        % obj.AddPrimer('ACTG....')
        % obj.AddPrimer([0 0 1 2 3 1 0 2 ...])
        
        % -------------------------------------------------------------------------
        % Check input sequence  
        % -------------------------------------------------------------------------
        if nargin < 1 | ~ismember(class(seq), {'char', 'int8'});
            error('matlabFunctions:invalidArguments', 'A valid char or int8 sequence must be provided');
        end
        if length(seq) ~= primerLength
            error('matlabFunctions:invalidArguments', 'The provided sequence must match the length of all primers');
        end
        
        % -------------------------------------------------------------------------
        % Convert if necessary and check the sequence space  
        % -------------------------------------------------------------------------
        switch class(seq)
            case 'char'
                % Convert
                intSeq = nt2int(obj.seqsToRemove{s}, 'ACGTOnly', true);
                if ~all(intSeq > 0)
                    error('matlabFunctions:invalidArguments', 'The sequence can only contain A, C, T, or G');
                end
                intSeq = int8(intSeq - 1); % Change basis and coerce to type

            case 'int8'
                % Confirm that only 0, 1, 2, 3 are provided
                if ~all(ismember(intSeq, [0 1 2 3]))
                    error('matlabFunctions:invalidArguments', 'All elements in integer sequences must be 0,1,2,or 3');
                end
                
            otherwise
                error('matlabFunctions:invalidArguments', 'Unrecognized sequence class');
        end
        
        % -------------------------------------------------------------------------
        % Add to the primer list and update the number of primers
        % -------------------------------------------------------------------------
        obj.seqs(end+1,:) = intSeq;
        obj.numPrimers = size(obj.seqs,1);
        
        % -------------------------------------------------------------------------
        % Recalculate Tm and GC and penalties
        % -------------------------------------------------------------------------
        obj.CalculatePrimerProperties(); %% IN FUTURE VERSIONS THIS COULD BE DONE MORE EFFICIENTLY
    end
        
    % -------------------------------------------------------------------------
    % Calculate Primer Properties  
    % -------------------------------------------------------------------------
    function CalculatePrimerProperties(obj, varargin)
        % Calculate the Tm, GC, and penalties of all sequences
        % obj.CalculatePrimerProperties()
        % obj.CalculatePrimerProperties(..., 'monovalentSalt',saltConcInM);
        % obj.CalculatePrimerProperties(..., 'primerConc',primerConcInM);

        % -------------------------------------------------------------------------
        % Handle variable input  
        % -------------------------------------------------------------------------
        defaults = cell(0,3); 
        defaults(end+1,:) = {'monovalentSalt', 'positive', []};
        defaults(end+1,:) = {'primerConc', 'positive', []};

        parameters = ParseVariableArguments(varargin, defaults, mfilename);
        
        % -------------------------------------------------------------------------
        % Update object salt and probe concentrations  
        % -------------------------------------------------------------------------
        if ~isempty(parameters.monovalentSalt)
            obj.monovalentSalt = parameters.monovalentSalt;
        end
        if ~isempty(parameters.primerConc)
            obj.primerConc = parameters.primerConc;
        end
        
        % -------------------------------------------------------------------------
        % Display Progress  
        % -------------------------------------------------------------------------
        if obj.verbose
            PageBreak();
            display(['Calculating Tm and GC for ' num2str(obj.numPrimers) ' primers with ' ...
                num2str(obj.monovalentSalt) ' M salt and ' num2str(obj.primerConc) ...
                ' M primer concentration']);
            timer = tic;
        end
        
        % -------------------------------------------------------------------------
        % Calculate primer GC  
        % -------------------------------------------------------------------------
        filterLength = obj.primerLength;
        
        obj.gc = filter(ones(1,filterLength)/filterLength, 1, double(obj.seqs == 1 | obj.seqs == 2), [], 2);
        obj.gc = obj.gc(:, filterLength:end);
        
        % -------------------------------------------------------------------------
        % Calculate Tm  
        % -------------------------------------------------------------------------
        obj.Tm = nan(obj.numPrimers,1);
        for i=1:obj.numPrimers
            % Derive enthalpy and entropy
            dG = TRDesigner.SantaLuciaNearestNeighbor(obj.seqs(i,:));
            dG = sum(dG,2);
            
            % Calculate 5' and 3' corrections
            fivePrimeAT = obj.seqs(i,1) == 0 | obj.seqs(i,1) == 3;
            threePrimeAT = obj.seqs(i,end) == 0 | obj.seqs(i,end) == 3;
            
            dG(1) = dG(1) + 0.2 + 2.2*fivePrimeAT + 2.2*threePrimeAT;
            dG(2) = dG(2) + -5.7 + 6.9*fivePrimeAT + 6.9*threePrimeAT;

            % Apply salt corrections
            dG(2) = dG(2) + 0.368*(obj.primerLength-1)*log(obj.monovalentSalt);
            
            % Calculate Tm
            obj.Tm(i) = dG(1)*1000 ./ (dG(2) + 1.9872 * log(obj.primerConc)) - 273.15;
            
            % NOTE: For the future, I should provide two concentrations,
            % probe and target, and this should be log(probeC - targetC/2)
            % where probeC > targetC. 
        end
        
        % -------------------------------------------------------------------------
        % Display Progress  
        % -------------------------------------------------------------------------
        if obj.verbose
            display(['... completed in ' num2str(toc(timer)) ' s']);
        end
        
        % -------------------------------------------------------------------------
        % Calculate Penalty Values  
        % -------------------------------------------------------------------------
        obj.penalties = nan(obj.numPrimers, length(obj.OTTables));
        for o=1:length(obj.OTTables)
            if obj.verbose
                PageBreak();
                display(['Calculating penalty for the ' obj.OTTableNames{o} ' table with ' ...
                  'seed length ' num2str(obj.OTTables(o).seedLength)]);
                timer = tic;
            end
            
            for s=1:obj.numPrimers 
                % Calculate total penalty for the sequence and for its
                % reverse complement
                obj.penalties(s,o) = sum(obj.OTTables(o).CalculatePenalty(obj.seqs(s,:))) + ...
                    sum(obj.OTTables(o).CalculatePenalty(fliplr(3-obj.seqs(s,:))));
            end
            
            if obj.verbose
                display(['... completed in ' num2str(toc(timer)) ' s']);
            end
        end
    end
    
    % -------------------------------------------------------------------------
    % Cut primers on GC, Tm, or penalty  
    % -------------------------------------------------------------------------
    function indsToKeep = CutPrimers(obj, varargin)
        % Cut primers based on their GC, Tm, or penalty
        % obj.CutPrimers(..., 'Tm', [low,up])
        % obj.CutPrimers(..., 'GC', [low,up])
        % obj.CutPrimers(..., 'OTTables', {'name', [low, up], 'name', range'})

        % -------------------------------------------------------------------------
        % Parse variable inputs
        % -------------------------------------------------------------------------
        defaults = cell(0,3); 
        defaults(end+1,:) = {'Tm', 'positive', []};
        defaults(end+1,:) = {'GC', 'positive', []};
        defaults(end+1,:) = {'OTTables', 'cell', {}};
        
        parameters = ParseVariableArguments(varargin, defaults, mfilename);

        % -------------------------------------------------------------------------
        % Display progress  
        % -------------------------------------------------------------------------
        if obj.verbose
            PageBreak();
            display('Keeping primers with ')
            if ~isempty(parameters.GC)
                display(['   GC in [' num2str(parameters.GC(1)) ...
                    ', ' num2str(parameters.GC(2)) ']']);
            end
            if ~isempty(parameters.Tm)
                display(['   Tm in [' num2str(parameters.Tm(1)) ...
                    ', ' num2str(parameters.Tm(2)) ']']);
            end
            if ~isempty(parameters.OTTables)
                for t=1:(length(parameters.OTTables)/2)
                    pid = find(strcmp(obj.OTTableNames, parameters.OTTables{2*(t-1)+1}), 1);
                    if isempty(pid)
                        error('matlabFunctions:invalidArguments', 'Unrecognized OTTable name');
                    end
                    display(['   ' obj.OTTableNames{pid} ' penalty in [' num2str(parameters.OTTables{2*t}(1)) ...
                    ', ' num2str(parameters.OTTables{2*t}(2)) ']']);
                end
            end
            timer = tic;
        end
        
        % -------------------------------------------------------------------------
        % Cut properties if needed
        % -------------------------------------------------------------------------
        indsToKeep = true(obj.numPrimers,1);
        if ~isempty(parameters.GC)
            indsToKeep = indsToKeep & obj.gc >= parameters.GC(1) & obj.gc <= parameters.GC(2);
        end
        if ~isempty(parameters.Tm)
            indsToKeep = indsToKeep & obj.Tm >= parameters.Tm(1) & obj.Tm <= parameters.Tm(2);
        end
        if ~isempty(parameters.OTTables)
            for t=1:(length(parameters.OTTables)/2)
                pid = find(strcmp(obj.OTTableNames, parameters.OTTables{2*(t-1)+1}), 1);
                if isempty(pid)
                    error('matlabFunctions:invalidArguments', 'Unrecognized OTTable name');
                end
                indsToKeep = indsToKeep & ...
                    obj.penalties(:,pid) >= parameters.OTTables{2*t}(1) & ...
                    obj.penalties(:,pid) <= parameters.OTTables{2*t}(2);
            end
        end
        
        % -------------------------------------------------------------------------
        % Cut oligos
        % -------------------------------------------------------------------------
        obj.seqs = obj.seqs(indsToKeep,:);
        obj.numPrimers = sum(indsToKeep);
        obj.gc = obj.gc(indsToKeep);
        obj.Tm = obj.Tm(indsToKeep);
        obj.penalties = obj.penalties(indsToKeep,:);

        % -------------------------------------------------------------------------
        % Display progress  
        % -------------------------------------------------------------------------
        if obj.verbose
            display(['... completed in ' num2str(toc(timer)) ' s']);
            display(['Removed ' num2str(sum(~indsToKeep)) ' primers']);
            display(['Kept ' num2str(sum(indsToKeep)) ' primers']);
        end
    end
    
    % -------------------------------------------------------------------------
    % Remove sequences that are not permitted
    % -------------------------------------------------------------------------
    function indsToKeep = RemoveForbiddenSeqs(obj, varargin)
        % Remove primers based on self complementarity
        % obj.RemoveForbiddenSeqs()
        % obj.RemoveForbiddenSeqs('seqsToRemove', {'seq1','seq2',...})
        
        % -------------------------------------------------------------------------
        % Handle variable input  
        % -------------------------------------------------------------------------
        defaults = cell(0,3); 
        defaults(end+1,:) = {'seqsToRemove', 'cell', {}};

        parameters = ParseVariableArguments(varargin, defaults, mfilename);

        % -------------------------------------------------------------------------
        % Update object properties  
        % -------------------------------------------------------------------------
        if ~isempty(parameters.seqsToRemove)
            obj.seqsToRemove = parameters.seqsToRemove;
        end
        
        % -------------------------------------------------------------------------
        % Display progress  
        % -------------------------------------------------------------------------
        if obj.verbose
            PageBreak();
            display(['Removing forbidden sequences']);
            timer = tic;
        end
        
        
        % -------------------------------------------------------------------------
        % Loop through forbidden sequences  
        % -------------------------------------------------------------------------
        indsToKeep = true(obj.numPrimers,1);
        for s=1:length(obj.seqsToRemove)
            % Display progress
            if obj.verbose
                display(['... finding ' obj.seqsToRemove{s}]);
            end
            
            % Convert forbidden sequence  
            intSeq = double(int8(nt2int(obj.seqsToRemove{s}, 'ACGTOnly', true))-1);
            if any(intSeq < 0)
                error('matlabFunctions:invalidArguments', 'Invalid forbidden sequence');
            end
            
            % Hash forbidden sequence
            hashBase = fliplr(4.^[0:(length(intSeq)-1)]);
            forbiddenHash = sum(hashBase.*intSeq) + 1;
            
            % Hash sequence
            obj.seqHash = filter(hashBase, 1, double(obj.seqs), [], 2) + 1;
            obj.seqHash = obj.seqHash(:,length(intSeq):end);
            
            % Find matches
            indsToKeep = indsToKeep & ~any(obj.seqHash == forbiddenHash,2);
        end
        
        % -------------------------------------------------------------------------
        % Update primers  
        % -------------------------------------------------------------------------
        obj.seqs = obj.seqs(indsToKeep,:);
        obj.numPrimers = size(obj.seqs,1);
        obj.gc = obj.gc(indsToKeep);
        obj.Tm = obj.Tm(indsToKeep);
        obj.penalties = obj.penalties(indsToKeep,:);

        % -------------------------------------------------------------------------
        % Display progress  
        % -------------------------------------------------------------------------
        if obj.verbose
            display(['... completed in ' num2str(toc(timer)) ' s']);
            display(['Removed ' num2str(sum(~indsToKeep)) ' primers']);
            display(['Kept ' num2str(sum(indsToKeep)) ' primers']);
        end
        
    end
    
    % -------------------------------------------------------------------------
    % Remove primers based on internal self complementarity
    % -------------------------------------------------------------------------
    function indsToKeep = RemoveSelfCompPrimers(obj, varargin)
        % Remove primers based on self complementarity
        % obj.RemoveSelfCompPrimers()
        % obj.RemoveSelfCompPrimers('maxHomology', homologyRegionLength)
        
        % -------------------------------------------------------------------------
        % Handle variable input  
        % -------------------------------------------------------------------------
        defaults = cell(0,3); 
        defaults(end+1,:) = {'homologyMax', 'positive', []};

        parameters = ParseVariableArguments(varargin, defaults, mfilename);
        
        % -------------------------------------------------------------------------
        % Update homology max property  
        % -------------------------------------------------------------------------
        if ~isempty(parameters.homologyMax)
            obj.homologyMax = parameters.homologyMax;
        end

        % -------------------------------------------------------------------------
        % Display progress  
        % -------------------------------------------------------------------------
        if obj.verbose
            PageBreak();
            display(['Identifying internal homology within each primer']);
            timer = tic;
        end

        % -------------------------------------------------------------------------
        % Hash the sequences
        % -------------------------------------------------------------------------
        hashBase = fliplr(4.^[0:(obj.homologyMax-1)]);
        obj.seqHash = filter(hashBase, 1, double(obj.seqs), [], 2) + 1;
        obj.seqHash = obj.seqHash(:,obj.homologyMax:end);
        
        % Hash the reverse complement: 3-seq = A=0->3=T, C=1->2=G
        obj.seqRCHash = filter(hashBase,1, double(fliplr(3-obj.seqs)), [], 2) + 1;
        obj.seqRCHash = obj.seqRCHash(:,obj.homologyMax:end);

        % -------------------------------------------------------------------------
        % Scan and remove primers based on internal homology
        % -------------------------------------------------------------------------
        indsToKeep = true(obj.numPrimers,1);
        for s=1:obj.numPrimers
            indsToKeep(s) = ~any(ismember(obj.seqHash(s,:), obj.seqRCHash(s,:)));
        end
        
        % -------------------------------------------------------------------------
        % Update primers  
        % -------------------------------------------------------------------------
        obj.seqs = obj.seqs(indsToKeep,:);
        obj.numPrimers = size(obj.seqs,1);
        obj.gc = obj.gc(indsToKeep);
        obj.Tm = obj.Tm(indsToKeep);
        obj.penalties = obj.penalties(indsToKeep,:);

        % -------------------------------------------------------------------------
        % Display progress  
        % -------------------------------------------------------------------------
        if obj.verbose
            display(['... completed in ' num2str(toc(timer)) ' s']);
            display(['Removed ' num2str(sum(~indsToKeep)) ' primers']);
            display(['Kept ' num2str(sum(indsToKeep)) ' primers']);
        end
        
    end
    
    % -------------------------------------------------------------------------
    % Remove primers based on cross homology  
    % -------------------------------------------------------------------------
    function indsToKeep = RemoveHomologousPrimers(obj, varargin)
        % Remove primers which share homologous regions
        % obj.RemoveHomologousPrimers()
        % obj.RemoveHomologousPrimers('maxHomology', homologyRegionLength)
        
        % -------------------------------------------------------------------------
        % Handle variable input  
        % -------------------------------------------------------------------------
        defaults = cell(0,3); 
        defaults(end+1,:) = {'homologyMax', 'positive', []};

        parameters = ParseVariableArguments(varargin, defaults, mfilename);
        
        % -------------------------------------------------------------------------
        % Update object salt and probe concentrations  
        % -------------------------------------------------------------------------
        if ~isempty(parameters.homologyMax)
            obj.homologyMax = parameters.homologyMax;
        end
        
        % -------------------------------------------------------------------------
        % Display progress  
        % -------------------------------------------------------------------------
        if obj.verbose
            PageBreak();
            display(['Identifying cross homology within primers']);
            timer = tic;
        end
        
        % -------------------------------------------------------------------------
        % Hash the sequences
        % -------------------------------------------------------------------------
        hashBase = fliplr(4.^[0:(obj.homologyMax-1)]);
        obj.seqHash = filter(hashBase, 1, double(obj.seqs), [], 2) + 1;
        obj.seqHash = obj.seqHash(:,obj.homologyMax:end);
        
        % Hash the reverse complement: 3-seq = A=0->3=T, C=1->2=G
        obj.seqRCHash = filter(hashBase,1, double(fliplr(3-obj.seqs)), [], 2) + 1;
        obj.seqRCHash = obj.seqRCHash(:,obj.homologyMax:end);

        % -------------------------------------------------------------------------
        % Identify unique hashes
        % -------------------------------------------------------------------------
        [uniqueHash, ~, ic] = unique([obj.seqHash(:) obj.seqRCHash(:)]);

        % -------------------------------------------------------------------------
        % Display progress  
        % -------------------------------------------------------------------------
        if obj.verbose
            display(['Identified ' num2str(length(uniqueHash)) ' unique ' ...
                num2str(obj.homologyMax) '-mers']);
        end
        
        % -------------------------------------------------------------------------
        % Create homology mat: non-zero entries indicate a shared sequence
        % -------------------------------------------------------------------------
        obj.homologyMat = sparse(size(obj.seqHash,1), size(obj.seqHash,1));
        
        % Loop over unique hash values
        for i=1:length(uniqueHash)
            localIDs = find(ic == i); % find inds of sequences that hold these values
            seqIDs = mod(localIDs-1, size(obj.seqHash,1))+1; % Map to sequence id
            % Loop over all sequence ids
            for k=1:length(seqIDs)
                for l=1:length(seqIDs)
                    if k~=l % Don't mark self homology
                        obj.homologyMat(seqIDs(k), seqIDs(l)) = 1;
                    end
                end
            end
            
            if obj.verbose && ~mod(i,1000)
                display(['... completed ' num2str(i) ' regions']);
            end
        end
        
        % -------------------------------------------------------------------------
        % Display progress  
        % -------------------------------------------------------------------------
        if obj.verbose
            display(['... completed in ' num2str(toc(timer)) ' s']);
        end
        
        % -------------------------------------------------------------------------
        % Display progress  
        % -------------------------------------------------------------------------
        if obj.verbose
            PageBreak();
            display(['Removing primers with cross homology']);
            timer = tic;
        end
        
        % -------------------------------------------------------------------------
        % Prepare row sum  
        % -------------------------------------------------------------------------
        rowSum = sum(obj.homologyMat,1);
        
        % -------------------------------------------------------------------------
        % Iterate until no cross homology  
        % -------------------------------------------------------------------------
        count = 1;
        while any(rowSum>0)
            [~, maxID] = max(rowSum); % Find primer with the most connections
            rowSum = rowSum - obj.homologyMat(maxID,:); % Delete its connections
            rowSum(maxID) = -Inf; % Flag as removed
            count = count + 1;
            if ~mod(count, 1000) && obj.verbose
                display(['... completed ' num2str(count) ' iterations']);
            end
        end
        
        % -------------------------------------------------------------------------
        % Update primers  
        % -------------------------------------------------------------------------
        indsToKeep = rowSum == 0;

        obj.seqs = obj.seqs(indsToKeep,:);
        obj.numPrimers = size(obj.seqs,1);
        obj.gc = obj.gc(indsToKeep);
        obj.Tm = obj.Tm(indsToKeep);
        obj.penalties = obj.penalties(indsToKeep,:);
        
        % -------------------------------------------------------------------------
        % Display progress  
        % -------------------------------------------------------------------------
        if obj.verbose
            display(['... completed in ' num2str(toc(timer)) ' s']);
            display(['Identified ' num2str(obj.numPrimers) ' orthogonal primers']);
        end
    end

    % -------------------------------------------------------------------------
    % Write fasta files
    % -------------------------------------------------------------------------
    function WriteFasta(obj, filePath, varargin)
        % Write a fasta file
        % obj.WriteFasta(filePath)
        % obj.WriteFasta(..., 'namePrefix', namePrefixValue)
        
        % -------------------------------------------------------------------------
        % Handle variable input  
        % -------------------------------------------------------------------------
        defaults = cell(0,3); 
        defaults(end+1,:) = {'namePrefix', 'string', ''};
        defaults(end+1,:) = {'fieldPad', 'positive', []};

        parameters = ParseVariableArguments(varargin, defaults, mfilename);

        % -------------------------------------------------------------------------
        % Generate random prefix  
        % -------------------------------------------------------------------------
        if isempty(parameters.namePrefix)
            import java.util.UUID;
            parameters.namePrefix = char(UUID.randomUUID());
            parameters.namePrefix = parameters.namePrefix(1:8);
        end
        
        % -------------------------------------------------------------------------
        % Determine field padding  
        % -------------------------------------------------------------------------
        if isempty(parameters.fieldPad)
            parameters.fieldPad = ceil(log10(obj.numPrimers));
        end
        
        % -------------------------------------------------------------------------
        % Check file path  
        % -------------------------------------------------------------------------
        if nargin < 1
            error('matlabFunctions:invalidArguments', 'A valid file path must be provided');
        end
        
        % -------------------------------------------------------------------------
        % Display progress  
        % -------------------------------------------------------------------------
        if obj.verbose
            PageBreak();
            display(['Writing fasta: ' filePath]);
            timer = tic;
        end

        % -------------------------------------------------------------------------
        % Delete existing files  
        % -------------------------------------------------------------------------
        if exist(filePath) == 2
            delete(filePath);
            if obj.verbose
                display('... Deleted existing file');
            end
        end
        
        % -------------------------------------------------------------------------
        % Build object  
        % -------------------------------------------------------------------------
        headers = {};
        seqs = {};
        for s=1:obj.numPrimers
            headers{s} = [parameters.namePrefix '-' num2str(s, ['%0' num2str(parameters.fieldPad) 'u']) ...
                ' Tm=' num2str(obj.Tm(s)) ...
                ' GC=' num2str(obj.gc(s))];
            seqs{s} = int2nt(obj.seqs(s,:)+1);
        end
        
        % -------------------------------------------------------------------------
        % Write fasta  
        % -------------------------------------------------------------------------
        fastawrite(filePath, headers, seqs);
        
        % -------------------------------------------------------------------------
        % Display progress  
        % -------------------------------------------------------------------------
        if obj.verbose
            display(['... completed in ' num2str(toc(timer)) ' s']);
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
        if ~isempty(p) & ~isa(p, 'parallel.Pool')
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
    % Save Function
    % -------------------------------------------------------------------------
    function Save(obj, dirPath)
        % Save the primer designer object in a directory specified by dirPath
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
    % Build a PrimerDesigner object from a saved version
    % -------------------------------------------------------------------------
    function obj = Load(dirPath)
        % obj = PrimerDesigner.Load(dirPath)
        
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
        obj = PrimerDesigner();
        
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
                otherwise
                    obj.(fieldsToLoad{i}) = LoadByteStream([dirPath fieldsToLoad{i} '.matb'], ...
                    'verbose', true);
            end
        end
    end
        
end % Static methods

end % Class definition
