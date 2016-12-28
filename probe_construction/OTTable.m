classdef OTTable < handle
% ------------------------------------------------------------------------
% [otTable, parameters] = OTTable(targetSequences, seedLength, varargin)
% This class returns a look up table that contains the penalty
% assigned to each seed sequence based on the frequency with which it appears
% in a set of off-target sequences. A seed sequence is every unique n-mer 
% (where n is defined by the seedLength).
% 
%--------------------------------------------------------------------------
% Necessary Inputs
% targetSequences -- A array of sequences in the fasta structure format,
%   i.e. each element needs a Header and a Sequence entry. Alternatively, 
%	this can be a transcriptome object or a cell array of sequences.
% seedLength -- the length of the seed sequences
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%  -- weights: An array with a weight to apply to each of the sequences in
%     targetSequences. The default is 1. These values can adjust for the
%     relative abundance of each of the targetSequences.
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% April 15, 2015; April 20, 2015
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Define properties
% -------------------------------------------------------------------------
properties
    name            % A string which defines the name of the OT Table, useful for indexing arrays of OTTables
    verbose         % A boolean that determines whether or not the classes displays progress
end

properties (SetAccess=protected)
    seedLength      			% The length of the seed region
    numPar          			% The number of parallel workers
    mapType         			% The type of map
    numEntries = 0  			% The number of entries in the map
    uniformWeight = false 		% Whether every sequence is weighted equally or not
    parallel        			% A parallel.Pool object
end

properties (Hidden=true, SetAccess=protected)
    data 						% Either a OTMap or OTMap2 instance, stores the key for each n-mer and the associated penalty
    hashBase					% The hash base used to quickly convert a n-mer sequence to an integer
end

% -------------------------------------------------------------------------
% Define methods
% -------------------------------------------------------------------------
methods
    
    % -------------------------------------------------------------------------
    % Define constructor
    % -------------------------------------------------------------------------
    function obj = OTTable(targetSequences, seedLength, varargin)
        % obj = OTTable(targetSequences, seedLength) creates an OTTable
        % using the targetSequences and a seedLength. targetSequences can
        % be a Transcriptome object, a structure with Header and Sequence
        % fields, or a cell array of sequences.
        %
        % obj = OTTable(..., 'weights', weights) % Specify the penalty
        % weight for each object
        %
        % obj = OTTable(..., 'mapType', mapType) % Specify the type of map
        %
        % obj = OTTable(..., 'parallel', parallel.pool instance)
        %
        % obj = OTTable(..., 'transferAbund', boolean) % Are the abundances
        % in the transcriptome object used as the penalty weights?
        %
        % obj = OTTable(..., 'name', string) % A name for the table
        
        % -------------------------------------------------------------------------
        % Default variables
        % -------------------------------------------------------------------------
        defaults = cell(0,3);

        % Parameters for parsing file names
        defaults(end+1,:) = {'weights', 'array', []};           % Weight to apply to each target sequence
        defaults(end+1,:) = {'verbose', 'boolean', false};      % Display progress of construction
        defaults(end+1,:) = {'mapType', 'function', @OTMap};    % The type of map that needs to be created.
        defaults(end+1,:) = {'convertMapType', 'boolean', true};% Whether to convert from OTMap to OTMap2
        defaults(end+1,:) = {'parallel', 'parallel', []};       % A parallel.pool object
        defaults(end+1,:) = {'transferAbund', 'boolean', false};% A boolean that specifies whether or not to add abundances from transcriptome obj
        defaults(end+1,:) = {'name', 'string', ''};             % A name for the OTTable instance

        % -------------------------------------------------------------------------
        % Parse variable input
        % -------------------------------------------------------------------------
        parameters = ParseVariableArguments(varargin, defaults, mfilename);
        
        % Transfer parameter values to object
        f = setdiff(fields(parameters), {'parallel', 'transferAbund', 'convertMapType', 'weights'});
        for i=1:length(f)
            obj.(f{i}) = parameters.(f{i});
        end        
        SetParallel(obj, parameters.parallel); % Set parallel pool object
 
        % -------------------------------------------------------------------------
        % Prepare empty object
        % -------------------------------------------------------------------------
        if nargin == 0
            return;
        end
        
        % -------------------------------------------------------------------------
        % Parse necessary input
        % -------------------------------------------------------------------------
        % Check number of required inputs
        if nargin < 2
            error('matlabFunctions:invalidArguments', ...
                ['A valid set of target sequences in fasta format or a transcriptome object as well as a valid seed length must be provided.']);
        end
        
        % Check properties of seed length
        if ~length(seedLength) == 1 && seedLength > 0
            error('matlabFunctions:invalidArgument', ...
                'Seed length must be a positive integer');
        end
		
        % Archive seed length and create hash base
        obj.seedLength = seedLength;
        obj.hashBase = fliplr(4.^[0:(obj.seedLength-1)]);
        obj.data = parameters.mapType(); % Create empty map of the desired class

        % Check for valid targetSequence arguments
        if isempty(targetSequences) % Handle empty class request
            return;
        else
            switch class(targetSequences)
                case 'Transcriptome'
                    % Do nothing this is a valid input and the class
                    % definition guarantees that all required attributes
                    % are present
                case 'struct' % This must be a structure with the output fields of fastaread
                    if ~isempty(setdiff({'Sequence'}, fields(targetSequences))) 
                        error('matlabFunctions:invalidArguments', ...
                    	'targetSequences must have a Sequence field.');
                    end
                case 'cell'
                    if ~all(strcmp(cellfun(@class, targetSequences, 'UniformOutput', false), 'char')) % Check that all entries are character arrays
                        error('matlabFunctions:invalidArguments', ...
                            'targetSequences provided as a cell array must all be character arrays.');
                    end
                otherwise
                    error('matlabFunctions:invalidArguments', ...
                        'targetSequences must be either a Transcriptome object, a structure with Header and Sequence fields, or a cell array of sequences.');e
            end
        end

        % -------------------------------------------------------------------------
        % Transfer abundances/weights if needed
        % -------------------------------------------------------------------------
        if isa(targetSequences, 'Transcriptome') % A transcriptome object was provided
            if parameters.transferAbund && targetSequences.abundLoaded % transfer abundances to weights
                parameters.weights = targetSequences.abundance;
            else
                parameters.weigths = [];
            end
            targetSequences = struct('Sequence', targetSequences.intSequences);
        end

        % -------------------------------------------------------------------------
        % Change format if provided sequences are a cell array
        % -------------------------------------------------------------------------
        if isa(targetSequences, 'cell')
            targetSequences = struct('Sequence', targetSequences);
        end
        
        % -------------------------------------------------------------------------
        % Check provided weights
        % -------------------------------------------------------------------------
        if isempty(parameters.weights)
            parameters.weights = ones(1, length(targetSequences));
            obj.uniformWeight = true;
        elseif length(parameters.weights) ~= length(targetSequences);
            error('matlabFunctions:invalidArguments', 'Weights must be equal in length to targetSequences');
        else
            obj.uniformWeight = false;
        end

        % -------------------------------------------------------------------------
        % Update internal number of parallel workers
        % -------------------------------------------------------------------------
        if isempty(obj.parallel)
            obj.numPar = 0;
        else
            obj.numPar = obj.parallel.NumWorkers;
        end
        
        % -------------------------------------------------------------------------
        % Display information on table construction
        % -------------------------------------------------------------------------
        if obj.verbose
            display('-------------------------------------------------------------------------');
            display(['Creating OT Table for ' num2str(length(targetSequences)) ' sequences '...
                'and a seed length of ' num2str(obj.seedLength)]);
            display(['Started on ' datestr(now)]);
            timer = tic;
        end
        
        % -------------------------------------------------------------------------
        % Allocate and build hash tables
        % -------------------------------------------------------------------------
        % Create separate key/value matrix for each worker
        data = Composite(obj.numPar);
        for p=1:length(data)
            data{p} = obj.mapType();
        end
        % Display progress
        if obj.verbose
            display(['Utilizing ' num2str(obj.numPar) ' parallel workers']);
        end
            
        % Create local variables to handle problem passing object
        % attributes into a parallel process
        seedLength = obj.seedLength;
        weights = parameters.weights;
        
        % Define properties for fast hash function
        hashBase = obj.hashBase;
        
        % Determine if nt2int conversion is needed
        performConversion = ischar(targetSequences(1).Sequence);
        
        % Loop through transcriptome in parallel
        spmd (obj.numPar)
            for t=labindex:numlabs:length(targetSequences) % Start at each workers id
                % Convert sequence to numbers
                if performConversion
                    localSeq = int8(nt2int(targetSequences(t).Sequence))-1;
                else
                    localSeq = targetSequences(t).Sequence; % Already int
                end

                % Hash sequence
                hash = filter(hashBase, 1, double(localSeq));
                hash = hash(seedLength:end);
                
                % Find >3 which corresponds to ambiguous nucleotides
                isValid = filter(ones(1, seedLength)/seedLength, ...
                    1, (localSeq > 3 | localSeq < 0)); 
                isValid = ~isValid(seedLength:end); % Kludge...
                
                % Add valid hashes with penalties to map
                AddToMap(data, cat(1,hash(isValid), zeros(1, sum(isValid)) + weights(t))); 
            end
        end
            
        % Add maps from parallel workers
        obj.data = data{1};
        data{1} = [];
        for j=2:length(data)
            obj.data.AddToMap(GetTable(data{j}));
        end

        % -------------------------------------------------------------------------
        % Convert MapType: OTMap is faster to build and OTMap2 is faster
        % for lookup
        % -------------------------------------------------------------------------
        if parameters.convertMapType && isa(obj.data, 'OTMap');
            obj.data = OTMap2(GetTable(obj.data));
            obj.mapType = @OTMap2;
        end
        
        % -------------------------------------------------------------------------
        % Display information on table construction
        % -------------------------------------------------------------------------
        if obj.verbose
            display(['Completed in  ' num2str(toc(timer)) ' s']);
        end
        
        % -------------------------------------------------------------------------
        % Collect properties of table
        % -------------------------------------------------------------------------
        obj.numEntries = length(obj.data);
        
    end

    % -------------------------------------------------------------------------
    % CheckValidity
    % -------------------------------------------------------------------------
    function isValid = IsValidSequence(obj,localSeq)
        if length(localSeq) < obj.seedLength
            isValid = 0;
            return;
        end
        if ischar(localSeq) % Coerce to integer if need be
            localSeq = int8(nt2int(localSeq, 'ACGTOnly', true))-1;
        end
        isValid = filter(ones(1, obj.seedLength)/obj.seedLength, ...
            1, (localSeq > 3 | localSeq < 0)); 
        isValid = ~isValid(obj.seedLength:end); % Kludge...
    end

    % -------------------------------------------------------------------------
    % Calculate Penalty
    % -------------------------------------------------------------------------
    function [penalty, hash] = CalculatePenalty(obj,seq)
		% Calculate the penalty associated with a given sequence. 
		% penalty = obj.CalculatePenalty(seq); 
	
        % Calculate number of seeds within the specified sequence
        numSeed = length(seq) - obj.seedLength + 1;
        
        % Check for short sequences
        if numSeed < 1
            penalty = [];
            hash = [];
            return;
        end
        
        % Convert only if needed
        if ischar(seq) 
            seq = int8(nt2int(seq, 'ACGTOnly', true))-1;
        end
        
        % Hash sequence
        hash = filter(obj.hashBase, 1, double(seq)); 
        hash = hash(obj.seedLength:end);
        
        % Find >3 which corresponds to ambiguous nucleotides
        isValid = filter(ones(1, obj.seedLength)/obj.seedLength, ...
            1, (seq >3 | seq < 0)); 
        isValid = ~isValid(obj.seedLength:end); % Kludge...
        
        % Calculate penalty associated with hash
        penalty = nan(1, numSeed); % Initialize as nan
        penalty(isValid) = obj.data.GetValues(hash(isValid));
    end
    
    % -------------------------------------------------------------------------
    % Return Keys
    % -------------------------------------------------------------------------
    function keyValues = keys(obj)
        keyValues = keys(obj.data);
    end
   
    % -------------------------------------------------------------------------
    % Return Penalty Values
    % -------------------------------------------------------------------------
    function penaltyValues = values(obj)
        penaltyValues = values(obj.data);
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
    % Overload of plus
    % -------------------------------------------------------------------------
    function obj3 = plus(obj1, obj2)
        % Add two OTTables
        % obj3 = obj1 + obj2
        % obj3 = plus(obj1, obj2)
        % obj3 = obj1.plus(obj2)
        
        % -------------------------------------------------------------------------
        % Check validity of objects
        % -------------------------------------------------------------------------
        if ~strcmp(class(obj1), 'OTTable') || ~strcmp(class(obj2), 'OTTable')
            error('matlabFunctions:invalidClass', ...
                'Cannot add OTTable with an object of a different class');
        end
        
        % -------------------------------------------------------------------------
        % Check compatibility of objects
        % -------------------------------------------------------------------------
        if obj1.seedLength ~= obj2.seedLength
            error('matlabFunctions:invalidAddition', ...
                'Cannot add two OTTables with different seed lengths.');
        end
        if ~strcmp(func2str(obj1.mapType), func2str(obj2.mapType))
            warning('matlabFunctions:invalidAddition', ...
                ['OTTables have different map types. Using the mapType of the first table.']);
        end
        
        % -------------------------------------------------------------------------
        % Create empty OTTable
        % -------------------------------------------------------------------------
        obj3 = OTTable([], obj1.seedLength, ...
            'mapType', obj1.mapType);
        
        % -------------------------------------------------------------------------
        % Combine data from each object
        % -------------------------------------------------------------------------
        % Get data tables
        data1 = GetTable(obj1.data);
        data2 = GetTable(obj2.data);
        
        % Find unique keys
        [keys, ~, inds] = unique(cat(2, data1(1,:), data2(1,:)));
        
        % Accumulate values for unique keys
        values = cat(2,data1(2,:), data2(2,:));
        values = accumarray(inds, values, [])';
        
        % -------------------------------------------------------------------------
        % Create new map
        % -------------------------------------------------------------------------
        obj3.data = obj1.mapType();
        obj3.data.AddToMap(cat(1,keys, values));
        obj3.numEntries = length(keys);
    end
    
    % -------------------------------------------------------------------------
    % Overload of sum
    % -------------------------------------------------------------------------
    function obj = sum(otTables)
        % obj = sum(otTables)
        
        % -------------------------------------------------------------------------
        % Check otTable array for validity
        % -------------------------------------------------------------------------
        if length(unique([otTables.seedLength])) ~= 1
            error('matlabFunctions:invalidArguments', ...
                'otTable arrays must have the same seed length for all elements');
        end
        if length(unique([otTables.uniformWeight])) ~= 1
            warning('matlabFunctions:invalidArguments', ...
                'Adding OTTables with different weighting');
        end

        % -------------------------------------------------------------------------
        % Compile all data tables
        % -------------------------------------------------------------------------
        % Display progress
        if otTables(1).verbose
            PageBreak();
            display(['Compiling all data tables']);
            tic;
        end
        
        % Allocate memory
        totalEntries = sum([otTables.numEntries]);
        allData = zeros(2, totalEntries);
        
        % Add tables
        count = 0;
        for i=1:length(otTables)
            localData = GetTable(otTables(i).data);
            if ~isempty(localData)
                allData(:, count + (1:size(localData,2))) = localData;
                count = count + size(localData,2);
            end
            
            % Display progress
            if ~mod(i,1000) && otTables(1).verbose
                display(['.... completed ' num2str(i) ' of ' num2str(length(otTables))]);
            end
            
        end
        
        % Display completion message
        if otTables(1).verbose
            display(['.... completed in ' num2str(toc) ' s']);
        end
        
        % -------------------------------------------------------------------------
        % Identify and sum redundant entries
        % -------------------------------------------------------------------------
        % Display progress
        if otTables(1).verbose
            display(['Adding redundant entries']);
            tic;
        end
        
        % Identify unique keys and sum values that share the same key
        [keys, ~, ic] = unique(allData(1,:));
        values = accumarray(ic, allData(2,:), [])';
        
        % Display completion message
        if otTables(1).verbose
            display(['.... completed in ' num2str(toc) ' s']);
        end
        
        % -------------------------------------------------------------------------
        % Create new otTable
        % -------------------------------------------------------------------------
        % Define seed length and uniformWeight
        seedLength = otTables(1).seedLength;
        uniformWeight = otTables(1).uniformWeight;
        
        % Create empty OTTable
        obj = OTTable([], seedLength, ...
            'verbose', otTables(1).verbose, ...
            'mapType', otTables(1).mapType);
        obj.uniformWeight = uniformWeight; % Transfer uniformWeight flag
        
        % Use the data table to create a new map
        obj.data = otTables(1).mapType([keys; values]);
        obj.numEntries = length(keys);
        
    end 
    
    % -------------------------------------------------------------------------
    % Save Function
    % -------------------------------------------------------------------------
    function Save(obj, dirPath)
        % Save the OTTable object (or array of OTTable objects) in a directory specified by dirPath
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
        end
        
        % -------------------------------------------------------------------------
        % Define fields to save
        % -------------------------------------------------------------------------
        fieldsToSave = properties(obj);
        fieldsToSave = setdiff(fieldsToSave, {'parallel', 'numPar', ...
            'mapType', 'verbose'}); % 'mapType' is not used again
        %fieldsToSave = union(fieldsToSave, {'data', 'hashBase'}); %hashBase is never used once the table is built
        fieldsToSave = union(fieldsToSave, {'data'});

        % -------------------------------------------------------------------------
        % Determine number of entries in object array (if it is an array)
        % -------------------------------------------------------------------------
        numObj = length(obj);
        
        % -------------------------------------------------------------------------
        % Save fields
        % -------------------------------------------------------------------------
        for i=1:length(fieldsToSave)
            switch fieldsToSave{i}
                case 'data' % Write data as a flat binary file/ numEntries will be the index for this flat file
                    fid = fopen([dirPath 'data.bin'], 'W');
                    if fid < 1
                        error('error opening file');
                    end
                    
                    % Record progress
                    if obj(1).verbose
                        display(['Saving ' dirPath 'data.bin']);
                        tic;
                    end
                    
                    % Loop over all entries
                    for i=1:numObj
                        fwrite(fid, GetTable(obj(i).data), 'double');
                    end
                    fclose(fid);
                    
                    % Record progress
                    if obj(1).verbose
                        display(['.... finished in ' num2str(toc)]);
                    end

                case {'name'} % Separate out fields that need to be saved as cells versus arrays
                    SaveAsByteStream([dirPath fieldsToSave{i} '.matb'], ...
                        {obj.(fieldsToSave{i})}, 'verbose', obj(1).verbose);
                otherwise % All other fields are entries that can be entered into an array
                    SaveAsByteStream([dirPath fieldsToSave{i} '.matb'], ...
                        [obj.(fieldsToSave{i})], 'verbose', obj(1).verbose);
            end
        end

    end
        
end

% -------------------------------------------------------------------------
% Static methods
% -------------------------------------------------------------------------
methods (Static)
    % -------------------------------------------------------------------------
    % Build a OTTable or an OTTable array from a saved version
    % -------------------------------------------------------------------------
    function obj = Load(dirPath, varargin)
        % obj = OTTable.Load(dirPath)
        % obj = OTTable.Load(..., 'verbose', verbose) % Determine verbosity of class
        % obj = OTTable.Load(..., 'mapType', mapType) % Determine mapType
        
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
        % Default variables
        % -------------------------------------------------------------------------
        defaults = cell(0,3);
        
        defaults(end+1,:) = {'verbose', 'boolean', true};       % Display progress of construction
        defaults(end+1,:) = {'mapType', 'function', @OTMap2};   % The type of map that needs to be created.
        
        % -------------------------------------------------------------------------
        % Parse variable input for constructor class
        % -------------------------------------------------------------------------
        parameters = ParseVariableArguments(varargin, defaults, mfilename);

        % -------------------------------------------------------------------------
        % Load number of entries and create empty OTTable array
        % -------------------------------------------------------------------------
        numEntries = LoadByteStream([dirPath 'numEntries.matb'], ...
                        'verbose', parameters.verbose);
                            
        % Parse out numEntries
        for j=1:length(numEntries)
            obj(j) = OTTable(); % Must create this way to avoid copy by reference issues
            obj(j).numEntries = numEntries(j);
        end

        % -------------------------------------------------------------------------
        % Load and deal generic fields
        % -------------------------------------------------------------------------
        fieldsToLoad = {'name', 'seedLength', 'uniformWeight'};
        for f=1:length(fieldsToLoad)
            % Define path to load
            fileToLoad = [dirPath fieldsToLoad{f} '.matb'];
            
            % Check for the possibility that no name field is defined --
            % backwards compatability
            if ~exist(fileToLoad) & strcmp(fieldsToLoad{f}, 'name')
                warning('matlabFunctions::missingInput', 'Did not find the name of the OTTable. Leaving blank');
                obj(j).name = '';
                continue;
            end
            
            % Load field 
            data = LoadByteStream(fileToLoad, ...
                        'verbose', parameters.verbose);
            % Parse out to object fields
            switch class(data)
                case 'cell'
                    for j=1:length(data)
                        obj(j).(fieldsToLoad{f}) = data{j};
                    end
                otherwise
                    for j=1:length(data)
                        obj(j).(fieldsToLoad{f}) = data(j);
                    end
            end
        end
        
        % -------------------------------------------------------------------------
        % Load and create maps
        % -------------------------------------------------------------------------
        if parameters.verbose
            display(['Loading ' dirPath 'data.bin']);
            tic;
        end
        
        % Open file
        fid = fopen([dirPath 'data.bin'], 'r');
        if fid < 1
            error('error opening file');
        end
        for i=1:length(numEntries)
            loadedData = fread(fid, [2, numEntries(i)], '*double');
            obj(i).data = parameters.mapType(loadedData);
        end
        % Close file
        fclose(fid);

        if parameters.verbose
            display(['.... finished in ' num2str(toc)]);
        end
        
        % -------------------------------------------------------------------------
        % Record map type and set hash base
        % -------------------------------------------------------------------------
        for i=1:length(obj)
            obj(i).mapType = parameters.mapType;
            obj(i).hashBase = fliplr(4.^[0:(obj(i).seedLength-1)]);
        end
    end
end

end


