 function MERFISHProbeDesign(varargin)



    %% ------------------------------------------------------------------------
    % Design of a MERFISH probe library
    % Jeffrey Moffitt
    % lmoffitt@mcb.harvard.edu
    % -------------------------------------------------------------------------
    % Copyright Presidents and Fellows of Harvard College, 2016.
    % Additional modifications by PRN, AIBS, 2019
    %% ------------------------------------------------------------------------
    % Setup the workspace
    %%-------------------------------------------------------------------------

    
    if nargin == 0
        % Call with no inputs to use this function as a script.
        % Input values in this block



        MERFISHAnalysisPath = 'C:\Users\ScanningLabAnalysis\Documents\MATLAB\MERFISH_analysis\';



        %% Set up paths
        % Raw data paths

            % MERFISHAnalysisPath is defined in the startup script for MERFISH_analysis;
            % MERFISH_Examples2 is a folder that contains several files required to
            % run this script.  These example files can be downloaded from http://zhuang.harvard.edu/merfish/MERFISHData/MERFISH_Examples2.zip


        % % Human (COS7 cells) MERFISH demo
        % basePath = [MERFISHAnalysisPath '\MERFISH_Examples2\']; % Base path where all required files can be found
        % analysisSavePath = SetFigureSavePath([basePath 'libraryDesign\'], 'makeDir', true);
        % rawTranscriptomeFasta = [basePath 'transcripts.fasta'];
        % fpkmPath = [basePath 'isoformsMod.fpkm_tracking'];
        % ncRNAPath = [basePath 'Homo_sapiens.GRCh38.ncrna.fa'];
        % readoutPath = [basePath 'readouts.fasta'];
        % codebookPath = [basePath 'codebook.csv'];

        % %-------------------------------------------------------------
        % % % Human smELT
        % basePath = 'D:\Data\MERFISH\Homosapiens\';
        % % basePath = [MERFISHAnalysisPath '\MERFISH_Examples2\']; % Base path where all required files can be found
        % % analysisSavePath = SetFigureSavePath([basePath 'libraryDesign\'], 'makeDir', true);
        % rawTranscriptomeFasta = [basePath 'Homo_sapiens_GRCh38_latest_rna.fna'];
        % fpkmPath = [basePath 'Homo_sapiens_proxyRandomFPKM.fpkm_tracking'];
        % ncRNAPath = [basePath 'Homo_sapiens.GRCh38.ncrna.fa'];
        % % readoutPath = [basePath 'readouts.fasta'];
        % % codebookPath = [basePath 'codebookHumansmELTv2.csv'];
        % 
        % readoutPath = 'D:\Data\MERFISH\Homosapiens\SMT-H-1002\SMT-H-1002_sequentialReadouts.fasta';
        % codebookPath = 'D:\Data\MERFISH\Homosapiens\SMT-H-1002\SMT-H-1002_sequential.csv';
        % 
        % analysisSavePath = SetFigureSavePath([basePath 'SMT-H-1002_isoSpec_70-100'], 'makeDir', true);
        % libraryName = ['SMT-H-1002_isoSpec_70-100'];

        % transcriptomeHeaderType = 'refseq'; % 'cufflinks', 'ensembl', or 'refseq'
        % transcriptomeIDType = 'NCBI'; % 'ENSEMBL' or 'NCBI'

        %-------------------------------------------------------------

        % % Mouse 
        species = 'Mus musculus';
        basePath = 'D:\Data\MERFISH\Musmusculus\'; % Base path where all required files can be found

        libraryName = ['SMT-M-1002_v2_noIsoSpecificity_70-100Specificity'];

        analysisSavePath = SetFigureSavePath([basePath, libraryName], 'makeDir', true);

        rawTranscriptomeFasta = [basePath 'Mus_musculus.GRCm38.cdna.all.fa'];
        fpkmPath = [basePath 'Mus_musculus_proxyRandomFPKM.fpkm_tracking'];
        ncRNAPath = [basePath 'Mus_musculus.GRCm38.ncrna.fa'];

        readoutPath = 'D:\Data\MERFISH\AIBSMERFISH_allReadouts.fasta';
        codebookPath = 'D:\Data\MERFISH\Musmusculus\SMT-M-1002_v2_noSpecificityv02\SMT-M-1002_v02_2hot.csv';

        transcriptomeHeaderType = 'ensembl'; % 'cufflinks', 'ensembl', or 'refseq'
        transcriptomeIDType = 'ENSEMBL'; % 'ENSEMBL' or 'NCBI'

        % -----------------------------------------------------------------------------

        % % % Human MERFISH + separate smELT, MTG
        % basePath = [MERFISHAnalysisPath 'MERFISH_Examples2\Homosapiens\']; % Base path where all required files can be found
        % % analysisSavePath = SetFigureSavePath([basePath 'libraryDesign\Musmusculus\'], 'makeDir', true);
        % 
        % 
        % 
        % rawTranscriptomeFasta = [basePath 'Homo_sapiens_GRCh38_latest_rna.fna'];
        % fpkmPath = [basePath 'Homo_sapiens_proxyRandomFPKM.fpkm_tracking'];
        % ncRNAPath = [basePath 'Homo_sapiens.GRCh38.ncrna.fa'];
        % % readoutPath = [basePath 'Human_MTG_Panel1_barcodedReadouts.fasta'];
        % % codebookPath = [basePath 'Human_MTG_Panel1_barcoded.csv'];
        % 
        % % readoutPath = [basePath 'Human_MTG_Panel1_sequentialReadouts.fasta'];
        % 
        % 
        % % codebookPath = [basePath 'Human_MTG_Panel1_sequential.csv'];



        %--------------------------------------
        % Input parameters:
        % Transcriptome:
        penaltyTableExactHomologyToExclude = 15;

        FPKMabundanceThreshold = 0;

        useUniformWeights = true;
        isoSpecificityTable_lengthOfExactHomology = 17;

        regionDesignParameters.regionLength = 30; % length of region (nt)
        regionDesignParameters.GC = [0 1]; % Fractional content
        regionDesignParameters.Tm = [70 100]; % deg C
        regionDesignParameters.isoSpecificity = [0 1]; % Within a gene
        regionDesignParameters.specificity = [.7 1]; % Across genes
        regionDesignParameters.monovalentSaltConcentration = 0.3; % mol/L
        regionDesignParameters.probeConcentration = 5e-9; % mol/L
        regionDesignParameters.probeSpacing = 3; % nt of gap between probes - set to negative to allow overlap

        numProbesPerGene = 48;


        primerDesignParameters.nPrimersToGenerate = 1e3;
        primerDesignParameters.primerLength = 20;
        primerDesignParameters.cutPrimers.Tm = [70 72];
        primerDesignParameters.cutPrimers.GC = [0.5 0.65];
        primerDesignParameters.cutPrimers.maxHomologySelf = 6;
        primerDesignParameters.cutPrimers.maxHomologyCross = 8;
        primerDesignParameters.primerMonovalentSaltConcentration = 0.3;
        primerDesignParameters.primerConcentration = 0.5e-6;
        

        targetRegionsFilter.geneIsoformListSource = 'allGenes';
        targetRegionsFilter.tRFilterMethod = 'default';
        targetRegionsFilter.tRFilterField = '';
        targetRegionsFilter.tRFilterParameters = [0, 1];
        
        % Indicate if this is allowed to ignore sequence version number in matching
        versionMatch = false;
        
        % For smELT probes where fewer than numProbesPerGene found, allow
        % adding second readout sequence to other end of probe 
        doubleHeadedsmELT = false;
        
        keepAllPossibleProbes = true;
        
        specifyReadouts = false;
        
        debugMode = false;
        
        readoutPermuteBySequence = true; % If true, use fileIO/permuteBySequence to define readout order on a sequence
                                         % If false, use default behavior
                                         % where randPerm defines readout
                                         % order.
        
    elseif (nargin == 1) && isa(varargin{1}, 'probeDesign')
        
        obj = varargin{1};
        % Input is probeDesign object
        
        % Call as method of that object
        % get all values from properties of that object
        
        MERFISHAnalysisPath = obj.MERFISHAnalysisPath;
        basePath = obj.basePath;

        libraryName = obj.libraryName;

        analysisSavePath = SetFigureSavePath([basePath, libraryName], 'makeDir', true);

        % Below block - if full path specified, use that; otherwise assume
        % relative to basePath
        if isempty(fileparts(obj.rawTranscriptomeFasta))
            rawTranscriptomeFasta = fullfile(basePath,  obj.rawTranscriptomeFasta);
        else
            rawTranscriptomeFasta = obj.rawTranscriptomeFasta;
        end
        
        if isempty(fileparts(obj.fpkmPath))
            fpkmPath = fullfile(basePath,  obj.fpkmPath);
        else
            fpkmPath = obj.fpkmPath;
        end
        
        if isempty(fileparts(obj.ncRNAPath))
            ncRNAPath = fullfile(basePath,  obj.ncRNAPath);
        else
            ncRNAPath = obj.ncRNAPath;
        end
        
        if isempty(fileparts(obj.readoutPath))
            readoutPath = fullfile(basePath,  obj.readoutPath);
        else
            readoutPath = obj.readoutPath;
        end
        
       if isempty(fileparts(obj.codebookPath))
            codebookPath = fullfile(basePath,  obj.codebookPath);
        else
            codebookPath = obj.codebookPath;
       end
        
       species = obj.species;
        
        transcriptomeHeaderType = obj.transcriptomeHeaderType;
        transcriptomeIDType = obj.transcriptomeIDType;
        
        penaltyTableExactHomologyToExclude = obj.penaltyTableExactHomologyToExclude;

        FPKMabundanceThreshold = obj.FPKMabundanceThreshold;

        useUniformWeights = obj.useUniformWeights;
        isoSpecificityTable_lengthOfExactHomology = obj.isoSpecificityTable_lengthOfExactHomology;

        regionDesignParameters.regionLength = obj.regionLength;
        regionDesignParameters.GC = obj.regionGC;
        regionDesignParameters.Tm = obj.regionTm;
        regionDesignParameters.isoSpecificity = obj.isoSpecificity;
        regionDesignParameters.specificity = obj.specificity;
        regionDesignParameters.monovalentSaltConcentration = obj.monovalentSaltConcentration;
        regionDesignParameters.probeConcentration = obj.probeConcentration;
        regionDesignParameters.probeSpacing = obj.probeSpacing;

        numProbesPerGene = obj.numProbesPerGene;


        primerDesignParameters.nPrimersToGenerate = obj.nPrimersToGenerate;
        primerDesignParameters.primerLength = obj.primerLength;
        primerDesignParameters.cutPrimers.Tm = obj.cutPrimersTm;
        primerDesignParameters.cutPrimers.GC = obj.cutPrimersGC;
        primerDesignParameters.cutPrimers.maxHomologySelf = obj.cutPrimersMaxHomologySelf;
        primerDesignParameters.cutPrimers.maxHomologyCross = obj.cutPrimersMaxHomologyCross;
        primerDesignParameters.primerMonovalentSaltConcentration = obj.primerMonovalentSaltConcentration;
        primerDesignParameters.primerConcentration = obj.primerConcentration;

        targetRegionsFilter.geneIsoformListSource = obj.geneIsoformListSource;
        targetRegionsFilter.tRFilterMethod = obj.tRFilterMethod;
        targetRegionsFilter.tRFilterField = obj.tRFilterField;
        targetRegionsFilter.tRFilterParameters = obj.tRFilterParameters;

        
        % Indicate if this is allowed to ignore sequence version number in matching
        versionMatch = obj.versionMatch;
        
        doubleHeadedsmELT = obj.doubleHeadedsmELT;
        
        keepAllPossibleProbes = obj.keepAllPossibleProbes;
        
        specifyReadouts = obj.specifyReadouts;
        
        debugMode = obj.debugMode;
        
        readoutPermuteBySequence = obj.readoutPermuteBySequence;
        
    else
        error('MERFISHProbeDesign inputs incorrect.');
        
    end

    %--------------------------------------

    % codebookPath = [basePath, 'M22E1_codebook.csv'];

    % transcriptomePath = 'C:\Users\ScanningLabAnalysis\Documents\MATLAB\MERFISH_analysis\MERFISH_Examples2\Musmusculus\150603_Transcriptome\transcriptomeC57Obj';


    % Paths at which to save created objects

    % analysisSavePath = SetFigureSavePath([basePath 'libraryDesign\Musmusculus'], 'makeDir', true);
    % % 
    % analysisSavePath = SetFigureSavePath('C:\Users\ScanningLabAnalysis\Documents\MATLAB\MERFISH_analysis\AIBSProbeDesign\HumanControlPanel1',...
    %                                         'makeDir', true);

    
    % The following objects can be re-used from a previous library, if ALL parameters for
    % a given object are identical in object and in provided path.
    % Most robust way would be to load file provided in path, double-check
    % that inputs in object are match.  If not, make new path since file
    % needs to be generated.
    % May be redundant with loading later if this is already loaded.
    % Add flags to skip loading below if already loaded here.
    
    
    
    specificityTablePathOK = false;
    
    trDesignerPathOK = false;
    
    %----------------------------------------------------------------------
    
    % Check if rRNAtRNA file provided can be reused
    rRNAtRNAPathOK = false;
    if (~isempty(obj.rRNAtRNAPath)) 
        if (exist(obj.rRNAtRNAPath, 'file') == 2) % File provided and path is file
    %         rRNAtRNA = fastaread(obj.rRNAtRNAPath);
            if true
                % Unfortunately ncRNA fasta file not saved for this part.  Also
                % this doesn't take long to generate.  So go ahead and re-make
                % if file not provided.
                rRNAtRNAPathOK = false;
                rRNAtRNAPath = [analysisSavePath 'rRNAtRNA.fa'];
            else 
                %
            end
        else % path is not to valid file
            rRNAtRNAPath = [analysisSavePath 'rRNAtRNA.fa']; % Regenerate path in this libraryName folder
            set(obj, 'rRNAtRNAPath', rRNAtRNAPath);
        end
    else % no file provided
        
        rRNAtRNAPath = [analysisSavePath 'rRNAtRNA.fa']; % Regenerate path in this libraryName folder
        set(obj, 'rRNAtRNAPath', rRNAtRNAPath);
    end
    
    
    %----------------------------------------------------------------------
    
    % Check if transcriptome provided can be reused
    transcriptomePathOK = false;
    if (~isempty(obj.transcriptomePath)) 
        if (exist(obj.transcriptomePath, 'dir') == 7)

            try 
                transcriptomePath = obj.transcriptomePath;
                transcriptome = Transcriptome.Load(transcriptomePath);

                checkVect = [strcmp(transcriptome.transPath,  obj.rawTranscriptomeFasta), ...     % raw transcriptome files match
                             strcmp(transcriptome.abundPath,  obj.fpkmPath), ...                  % FPKM files match
                             strcmp(transcriptome.headerType,  obj.transcriptomeHeaderType), ...  % Header type matches
                             strcmp(transcriptome.IDType,  obj.transcriptomeIDType)];             % ID types match
                if  all(checkVect)
                   % File and input parameters of object are same
                   % Use loaded file
                   transcriptomePathOK = true;
                   fprintf(1, 'Using existing transcriptome object.\n');
                else
                    fprintf(1, 'transcriptome object check failed.  Setting path to default\n');
                end
            catch
                fprintf(1, 'transcriptome object load error.  Setting path to default\n');
            end
        else
            fprintf(1, 'transcriptome object path not valid.  Setting path to default\n');
        end
    else
        fprintf(1, 'transcriptome object path not specified.  Setting path to default\n');
    end
    
    if ~transcriptomePathOK
        clear transcriptome
        transcriptomePath = [analysisSavePath 'transcriptomeObj'];
        set(obj, 'transcriptomePath', transcriptomePath);
    end
    
        
    
    %----------------------------------------------------------------------
    
    
    % Check if isospecificityTable provided can be reused
    isoSpecificityTablePathOK = false;
    if (~isempty(obj.isoSpecificityTablePath)) 
        if (exist(obj.isoSpecificityTablePath, 'dir') == 7)
        
            isoSpecificityTablePath = obj.isoSpecificityTablePath;
            try
                isoSpecificityTables = OTTable.Load(isoSpecificityTablePath, ...
                    'verbose', true, ...
                    'mapType', @OTMap2);

                tNames = transcriptome.GetNames();
                iNamesSame = false(length(isoSpecificityTables), 1);
                for iN = 1:length(isoSpecificityTables)
                    iNamesSame(iN) = strcmp(tNames{iN}, strtrim(isoSpecificityTables(iN).name));
                end
                checkVect = [length(tNames) == length(isoSpecificityTables), ...                                        %length matches
                             all(iNamesSame), ...                                                                       %all names match
                             isoSpecificityTables(1).seedLength == obj.isoSpecificityTable_lengthOfExactHomology, ...   %seed length matches
                             isoSpecificityTables(1).uniformWeight == obj.useUniformWeights];                           %use of uniformWeight matches

                if  all(checkVect)
                   % File and input parameters of object are same
                   % Use loaded file
                   isoSpecificityTablePathOK = true;
                   fprintf(1, 'Existing isoSpecificityTable to be used.\n')
                else
                    fprintf(1, 'isoSpecificityTable check failed. Setting path to default.\n')
                end
            catch
                fprintf(1, 'isoSpecificityTable load error.  Setting path to default\n');
            end
        else
            fprintf(1, 'isoSpecificityTable path not valid.  Setting path to default\n');
        end
    else
        fprintf(1, 'isoSpecificityTable path not specified.  Setting path to default\n');
    end
    
    if ~isoSpecificityTablePathOK
        clear isoSpecificityTable
        isoSpecificityTablePath = [analysisSavePath 'isoSpecificityTable'];
        set(obj, 'isoSpecificityTablePath', isoSpecificityTablePath);
    end
    
    %----------------------------------------------------------------------
    % Check if specificityTable provided can be reused
    specificityTablePathOK = false;
    if (~isempty(obj.specificityTablePath)) 
        if (exist(obj.specificityTablePath, 'dir') == 7) && isoSpecificityTablePathOK  
            
            specificityTablePath = obj.specificityTablePath;
            try
                specificityTables = OTTable.Load(specificityTablePath, ...
                    'verbose', true, ...
                    'mapType', @OTMap2);

                isoSpecCount = 0;
                for iN = 1:length(isoSpecificityTabales)
                    isoSpecCount = isoSpecCount + isoSpecificityTables(iN).numEntries;
                end
                checkVect = [isoSpecCount == specificityTables.numEntries];  % ???? Is this the right comparision to check if specPath = specObj?

                if  all(checkVect)
                   % File and input parameters of object are same
                   % Use loaded file
                   specificityTablePathOK = true;
                   fprintf(1, 'Using existing specificityTable.\n');
                else
                    fprintf(1, 'specificityTable check failed. Setting path to default.\n');
                end
            catch
                fprintf(1, 'specificityTable load error.  Setting path to default\n');
            end
        else
            fprintf(1, 'specificityTable not valid.  Setting path to default\n');
        end
    else
        fprintf(1, 'specificityTable path not specified.  Setting path to default\n');
    end
    
    if ~specificityTablePathOK
        clear specificityTables
        specificityTablePath = [analysisSavePath 'specificityTable'];
        set(obj, 'specificityTablePath', specificityTablePath);      
    end
    
    %----------------------------------------------------------------------
    % Check if trDesigner provided can be reused
    trDesignerPathOK = false;
    if (~isempty(obj.trDesignerPath)) 
        if (exist(obj.trDesignerPath, 'dir') == 7) && isoSpecificityTablePathOK   && specificityTablePathOK 
        
            trDesignerPath = obj.trDesignerPath;
            try
                trDesigner = TRDesigner.Load(trDesignerPath);

                checkVect = [false];  % This info not in object.  No way to validate.  

                if  all(checkVect)
                   % File and input parameters of object are same
                   % Use loaded file
                   trDesignerPathOK = true;
                   fprintf(1, 'Using existing trDesigner.\n');
                else
                    fprintf(1, 'trDesigner check failed. Using default path.\n');
                end
            catch
                fprintf(1, 'trDesigner load error.  Setting path to default\n');
            end
        else
            fprintf(1, 'trDesigner not valid.  Setting path to default\n');
        end
    else
        fprintf(1, 'trDesigner path not specified.  Setting path to default\n');

    end
    
    if ~trDesignerPathOK
        clear trDesigner
        trDesignerPath = [analysisSavePath 'trDesigner'];
        set(obj, 'trDesignerPath', trDesignerPath);
    end
    
    % trRegionsPath =  [analysisSavePath 'tr_GC_43_63_Tm_66_76_Len_30_30_IsoSpec_0.75_1_Spec_0.75_1'];


    trRegionsFileName = sprintf('tr_GC_%d_%d_Tm_%d_%d_Len_%d_%d_IsoSpec_%.2f_%.2f_Spec_%.2f_%.2f', regionDesignParameters.GC(1)*100, ...
                                                                                           regionDesignParameters.GC(2)*100, ...
                                                                                           regionDesignParameters.Tm(1), ...
                                                                                           regionDesignParameters.Tm(2), ...
                                                                                           regionDesignParameters.regionLength(1), ...
                                                                                           regionDesignParameters.regionLength(1), ...
                                                                                           regionDesignParameters.isoSpecificity(1), ...
                                                                                           regionDesignParameters.isoSpecificity(2), ...
                                                                                           regionDesignParameters.specificity(1), ...
                                                                                           regionDesignParameters.specificity(2));

    

    % Always rebuild these paths. Do not use from object.                                                                                   
    trRegionsPath =  [analysisSavePath trRegionsFileName];
    usedReadoutPath = [analysisSavePath libraryName '_used_readouts.fasta'];
    possibleOligosPath = [analysisSavePath libraryName '_possible_oligos.fasta'];
    oligosPath = [analysisSavePath libraryName '_oligos.fasta'];
    finalPrimersPath = [analysisSavePath libraryName '_primers.fasta'];
    primersPath = [analysisSavePath libraryName '_possible_primers.fasta'];
    allOligosPath = [analysisSavePath libraryName '_AllOligos.fasta'];
    
    % Check to see if any of the required or generated files already exist.
    % If any of these exist it would cause function to exit, so check
    % before code gets deep enough to be annoying when it stops.
    
    pathsToCheck = {rRNAtRNAPath, 'load'; ...
                    transcriptomePath, 'load'; ...
                    specificityTablePath, 'load'; ...
                    isoSpecificityTablePath, 'load'; ...
                    trDesignerPath, 'load'; ...
                    trRegionsPath, 'load'; ...
                    readoutPath, 'crashIfNone'; ...
                    codebookPath, 'crashIfNone'; ...
                    usedReadoutPath, 'delete'; ...
                    possibleOligosPath, 'crashIfPresent'; ...
                    oligosPath, 'crashIfPresent'; ...
                    primersPath, 'crashIfPresent'; ...
                    finalPrimersPath, 'crashIfPresent'; ...
                    allOligosPath, 'crashIfPresent'};

    % Check that paths generated exist or not as required by downstream
    % code
    badFilesFound = filePathCheck(pathsToCheck);
    
    if badFilesFound
        error('Address file path issues and restart processing.');
    end
                   
    % Start logging
    logFileName = fullfile(analysisSavePath, strcat(libraryName, '.log'));
    logFID = fopen(logFileName, 'w+');
    try

        fprintf(logFID, '%s log file\n', libraryName);
        fprintf(logFID, 'basePath : %s\n', basePath);
        fprintf(logFID, '-----------------------------\n');
        fprintf(logFID, 'species : %s\n', species);
        fprintf(logFID, 'rawTranscriptomeFasta : %s\n', rawTranscriptomeFasta);
        fprintf(logFID, 'fpkmPath : %s\n', fpkmPath);
        fprintf(logFID, 'ncRNAPath : %s\n', ncRNAPath);
        fprintf(logFID, 'readoutPath : %s\n', readoutPath);
        fprintf(logFID, 'codebookPath : %s\n', codebookPath);
        fprintf(logFID, 'fpkmPath : %s\n', fpkmPath);
        fprintf(logFID, 'rRNAtRNAPath : %s\n', rRNAtRNAPath);
        fprintf(logFID, 'transcriptomePath : %s\n', transcriptomePath);
        fprintf(logFID, 'specificityTablePath : %s\n', specificityTablePath);
        fprintf(logFID, 'isoSpecificityTablePath : %s\n', isoSpecificityTablePath);
        fprintf(logFID, 'trDesignerPath : %s\n', trDesignerPath);
        fprintf(logFID, 'trRegionsPath : %s\n', trRegionsPath);
        fprintf(logFID, 'usedReadoutPath : %s\n', usedReadoutPath);
        fprintf(logFID, 'possibleOligosPath : %s\n', possibleOligosPath);
        fprintf(logFID, 'oligosPath : %s\n', oligosPath);
        if keepAllPossibleProbes
            fprintf(logFID, 'allOligosPath : %s\n', allOligosPath);
        end
        fprintf(logFID, '-----------------------------\n');
        fprintf(logFID, 'transcriptomeHeaderType : %s\n', transcriptomeHeaderType);
        fprintf(logFID, 'transcriptomeIDType : %s\n', transcriptomeIDType);
        fprintf(logFID, 'transcriptome-codebook_version_match_enforced : %d\n', versionMatch);
        fprintf(logFID, '\n');
        fprintf(logFID, 'penaltyTableExactHomologyToExclude : %d\n', penaltyTableExactHomologyToExclude);
        fprintf(logFID, '\n');
        fprintf(logFID, 'FPKMabundanceThreshold : %d\n', FPKMabundanceThreshold);
        fprintf(logFID, '\n');
        fprintf(logFID, 'useUniformWeights : %d\n', useUniformWeights);
        fprintf(logFID, 'isoSpecificityTable_lengthOfExactHomology : %d\n', isoSpecificityTable_lengthOfExactHomology);
        fprintf(logFID, '\n');
        fprintf(logFID, 'regionDesignParameters\n');
        fprintf(logFID, 'regionLength : %d\n', regionDesignParameters.regionLength);
        fprintf(logFID, 'GC : [%d, %d]\n', regionDesignParameters.GC);
        fprintf(logFID, 'Tm : [%d, %d]\n', regionDesignParameters.Tm);
        fprintf(logFID, 'isoSpecificity : [%d, %d]\n', regionDesignParameters.isoSpecificity);
        fprintf(logFID, 'specificity : [%d, %d]\n', regionDesignParameters.specificity);
        fprintf(logFID, 'monovalentSaltConcentration : %d\n', regionDesignParameters.monovalentSaltConcentration);
        fprintf(logFID, 'probeConcentration : %d\n', regionDesignParameters.probeConcentration);
        fprintf(logFID, 'probeSpacing : %d\n', regionDesignParameters.probeSpacing);
        fprintf(logFID, '\n');
        fprintf(logFID, 'numProbesPerGene : %d\n', numProbesPerGene);
        fprintf(logFID, '\n');
        fprintf(logFID, 'primerDesignParameters\n');
        fprintf(logFID, 'nPrimersToGenerate : %d\n', primerDesignParameters.nPrimersToGenerate);
        fprintf(logFID, 'primerLength : %d\n', primerDesignParameters.primerLength);
        fprintf(logFID, 'GC : [%d, %d]\n', primerDesignParameters.cutPrimers.GC);
        fprintf(logFID, 'Tm : [%d, %d]\n', primerDesignParameters.cutPrimers.Tm);
        fprintf(logFID, 'primerMonovalentSaltConcentration : %d\n', primerDesignParameters.primerMonovalentSaltConcentration);
        fprintf(logFID, 'primerConcentration : %d\n', primerDesignParameters.primerConcentration);
        fprintf(logFID, 'maxHomologySelf : %d\n', primerDesignParameters.cutPrimers.maxHomologySelf);
        fprintf(logFID, 'maxHomologyCross : %d\n', primerDesignParameters.cutPrimers.maxHomologyCross);
        fprintf(logFID, 'doubleHeadedsmELT : %d\n', doubleHeadedsmELT);
        fprintf(logFID, 'keepAllPossibleProbes : %d\n', keepAllPossibleProbes);
        fprintf(logFID, 'specifyReadouts : %d\n', specifyReadouts);
        fprintf(logFID, 'readoutPermuteBySequence : %d\n', readoutPermuteBySequence);
        fprintf(logFID, 'geneIsoformListSource : %s\n', geneIsoformListSource);
        fprintf(logFID, 'tRFilterMethod : %s\n', tRFilterMethod);
        fprintf(logFID, 'tRFilterField : %s\n', tRFilterField);
        fprintf(logFID, 'tRFilterParameters : [%d, %d]\n', tRFilterParameters);
        
        
        fprintf(logFID, '-----------------------------\n');



        %% ------------------------------------------------------------------------
        % Step 1: Construct possible target regions 
        %  Below a set of transcripts, transcript abundances, and non-coding RNA 
        %  sequences will be used to design a set of possible target regions for all
        %  transcripts in the human transcriptome.
        %%-------------------------------------------------------------------------

        %% ------------------------------------------------------------------------
        % Load ncRNAs and cut to rRNA and tRNA: Needed to eliminate probes that 
        %   have homology to these abundant RNAs
        %%-------------------------------------------------------------------------
        %% Load and parse ncRNAs
        PageBreak();
        if ~exist(rRNAtRNAPath)
            display(['Loading: ' ncRNAPath]);
            fprintf(logFID, '%s - Loading %s\n', datestr(datetime), ncRNAPath);
            tic;
            ncRNAs = fastaread(ncRNAPath);
            display(['... completed in ' num2str(toc) ' s']);
            fprintf(logFID, '%s - ... completed in %s s\n', datestr(datetime), num2str(toc));
            display(['Found ' num2str(length(ncRNAs)) ' sequences']);
            fprintf(logFID, '%s - Found %d sequences\n', datestr(datetime), length(ncRNAs));

            % Parse out 'gene_biotype'
            biotypes = {};
            for i=1:length(ncRNAs)
                tempString= regexp(ncRNAs(i).Header, 'gene_biotype:\S+ ', 'match');
                strParts = strsplit(tempString{1}, {':', ' '});
                biotypes{i} = strParts{2};
            end

            % Identify features to keep
            biotypesToKeep = {'rRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA'};
            PageBreak();
            display(['Keeping the following types: ']);
            fprintf(logFID, '%s - Keeping the following types:\n', datestr(datetime));

            for i=1:length(biotypesToKeep)
                display(['   ' biotypesToKeep{i}]);
                fprintf(logFID, '%s\n', biotypesToKeep{i});
            end

            indsToKeep = ismember(biotypes, biotypesToKeep);

            rRNAtRNA = ncRNAs(indsToKeep);
            display(['Keeping ' num2str(length(rRNAtRNA)) ' ncRNAs']);
            fprintf(logFID, '%s - Keeping %d ncRNAs\n', datestr(datetime), length(rRNAtRNA));

            if exist(rRNAtRNAPath)
                delete(rRNAtRNAPath);
            end

            % Save ncRNA in a fasta file
            fastawrite(rRNAtRNAPath, rRNAtRNA);
            display(['Wrote: ' rRNAtRNAPath]);
            fprintf(logFID, '%s - Wrote %s\n', datestr(datetime), rRNAtRNAPath);
            
        elseif rRNAtRNAPathOK
            % Already loaded file above when doing path check
            display(['Previously loaded: ' rRNAtRNAPath]);
            fprintf(logFID, '%s - Previously loaded %s\n', datestr(datetime), rRNAtRNAPath);
            display(['Loaded ' num2str(length(rRNAtRNA)) ' sequences']);
            fprintf(logFID, '%s - Loaded %d sequences\n', datestr(datetime), length(rRNAtRNA));
        else
            % Load existing file if already created
            display(['Found and loading: ' rRNAtRNAPath]);
            fprintf(logFID, '%s - Found and loading %s\n', datestr(datetime), rRNAtRNAPath);
            tic;
            rRNAtRNA = fastaread(rRNAtRNAPath);
            display(['.... completed in ' num2str(toc) 's']);
            fprintf(logFID, '%s - ... completed in %s s\n', datestr(datetime), num2str(toc));
            display(['Loaded ' num2str(length(rRNAtRNA)) ' sequences']);
            fprintf(logFID, '%s - Loaded %d sequences\n', datestr(datetime), length(rRNAtRNA));
        end

        %% ------------------------------------------------------------------------
        % Build transcriptome object: This object collects information about
        %   all transcripts in the transcriptome so as to facilitate acess to 
        %   various properties.
        %%-------------------------------------------------------------------------
        %% Build transcriptome object


        if ~exist(transcriptomePath)

            fprintf(logFID, '%s - Building transcriptome object.\n', datestr(datetime));

            %Build transcriptome using existing abundance data
            transcriptome = Transcriptome(rawTranscriptomeFasta, ...
                'abundPath', fpkmPath, ...
                'verbose', true, ...
                'headerType', transcriptomeHeaderType, ...
                'IDType', transcriptomeIDType);
            %transcriptome = Transcriptome(rawTranscriptomeFasta, ...
             %   'abundPath', fpkmPath, ...
              %  'verbose', true);
            transcriptome.Save(transcriptomePath);
            fprintf(logFID, '%s - Transcriptome object saved to %s\n', datestr(datetime), transcriptomePath);
        
        elseif transcriptomePathOK % Loaded at path checking step
            fprintf(logFID, '%s - Prevously loaded transcriptome object from %s\n', datestr(datetime), transcriptomePath);
            fprintf(1, '%s - Prevously loaded transcriptome object from %s\n', datestr(datetime), transcriptomePath);
        
        else
            % Load transcriptome if it already exists
            fprintf(logFID, '%s - Loading transcriptome object from %s\n', datestr(datetime), transcriptomePath);
            transcriptome = Transcriptome.Load(transcriptomePath);
            fprintf(logFID, '%s - Transcriptome object loaded.\n', datestr(datetime));
        end

        %% ------------------------------------------------------------------------
        % Build specificity tables:  These tables are used to identify (and 
        %   penalize probes that contain) potential regions of cross-homology 
        %   between transcripts
        %%-------------------------------------------------------------------------
        %% Build isoform specificity table -- these tables, one per gene, calculate 
        %  the penalty associated with homology regions between isoforms of the 
        %  the same gene

        if ~exist(isoSpecificityTablePath)
            % Get isoform data
            names = transcriptome.GetNames();
            idsByName = transcriptome.GetIDsByName(names);
            tic;

            % Display progress header
            PageBreak();
            display('Starting construction of isoform specificity tables');
            fprintf(logFID, '%s - Constructing isoform specificity tables.\n', datestr(datetime));
            display(['Started on ' datestr(now)]);

            clear isoSpecificityTables

            for i=1:length(names) % Loop over all gene names -- RNAs that share the same gene name are considered isoforms
                % Generate a local transcriptome object that contains transcripts for a single gene
                localTranscriptome = transcriptome.Slice('geneID', idsByName{i});

                if mod(i, 100) == 0
                    fprintf(1, 'On gene %d of %d. First transcript %s\n', i, length(names), idsByName{i}{1});
                    fprintf(logFID, '%s - Gene %d of %d. First transcript %s\n', datestr(datetime), i, length(names), idsByName{i}{1});
                end

                if useUniformWeights

                    % Generate a OTTable for isoforms for the given gene
                    isoSpecificityTables(i) = OTTable(localTranscriptome, ...
                        isoSpecificityTable_lengthOfExactHomology, ...  % lengthOfExactHomology is the length of exact homology used to calculate penalties
                        'verbose', false, ...                        
                        'transferAbund', false);

                else

                    % Generate a OTTable for isoforms for the given gene
                    isoSpecificityTables(i) = OTTable(localTranscriptome, ...
                        isoSpecificityTable_lengthOfExactHomology, ...  % lengthOfExactHomology is the length of exact homology used to calculate penalties
                        'verbose', false, ...
                        'transferAbund', true);

                end

                % Name the table
                isoSpecificityTables(i).name = names{i}; 

                % Display progress
                if ~mod(i,500) 
                    display(['... completed ' num2str(i) ' of ' num2str(length(names)) ' genes']);
                    fprintf(logFID, '%s - Completed %d of %d genes.\n', datestr(datetime), i, length(names));
                end
            end
            display(['.... completed in ' num2str(toc) ' s']);
            fprintf(logFID, '%s - .\n', datestr(datetime), i, length(names));
            % Save tables
            isoSpecificityTables.Save(isoSpecificityTablePath);
            fprintf(logFID, '%s - isoSpecificity table saved to %s.\n', datestr(datetime), isoSpecificityTablePath);
        
        elseif isoSpecificityTablePathOK % Loaded at path checking step
            fprintf(logFID, '%s - Prevously loaded isoSpecificity tables from %s\n', datestr(datetime), isoSpecificityTablePath);            
            fprintf(1, '%s - Prevously loaded isoSpecificity tables from %s\n', datestr(datetime), isoSpecificityTablePath); 
        else
            fprintf(logFID, '%s - Loading isoSpecificity table from %s.\n', datestr(datetime), isoSpecificityTablePath);
            % Load tables if they already exist
            isoSpecificityTables = OTTable.Load(isoSpecificityTablePath, ...
                'verbose', true, ...
                'mapType', @OTMap2);
            fprintf(logFID, '%s - isoSpecificity table loaded!\n', datestr(datetime));
        end

        %% Build total specificity table --- this table contains a penalty associated
        %   with all possible sequences in the transcriptome
        if ~exist(specificityTablePath)

            fprintf(logFID, '%s - Generating specificty table\n', datestr(datetime));

            isoSpecificityTables(1).verbose = true; % Raise flag to display progress

            % Add isoform specificity tables to create total transcriptome
            % specificity table
            specificityTable = sum(isoSpecificityTables);  % Composed by summing the penalties calculated for all isoforms

            % Name the table
            specificityTable.name = 'Transcriptome Specificity';

            % Reset verbosity of isoform specificity table
            isoSpecificityTables(1).verbose = false;

            % Save table
            specificityTable.verbose = true;
            specificityTable.Save(specificityTablePath);
            fprintf(logFID, '%s - Specificty table saved to %s\n', datestr(datetime), specificityTablePath);
        
        elseif specificityTablePathOK
            fprintf(logFID, '%s - Specificty table previously loaded from %s\n', datestr(datetime), specificityTablePath);
            fprintf(1, '%s - Specificty table previously loaded from %s\n', datestr(datetime), specificityTablePath);
        else
            % Load table if it exists
            specificityTable = OTTable.Load(specificityTablePath, ...
                'verbose', true, ...
                'mapType', @OTMap2);
            fprintf(logFID, '%s - Specificty table loaded from %s\n', datestr(datetime), specificityTablePath);
        end

        %% ------------------------------------------------------------------------
        % Configure parallel.Pool
        %%-------------------------------------------------------------------------
        %% Create parallel pool... speeds up the construction of the TRDesigner and the construction of libraries
        if isempty(gcp('nocreate'))
            fprintf(logFID, '%s - Start parallel processing pool\n', datestr(datetime));
            p = parpool(4);  % Insert a number here appropriate to the used computational resources
        else
            p = gcp;
        end

        parHold = p;

        %% ------------------------------------------------------------------------
        % Build Penalty Tables
        %%-------------------------------------------------------------------------
        %% Build rRNA/tRNA penalty table

        fprintf(logFID, '%s - Generate rRNA+tRNA penalty table\n', datestr(datetime));

        OTrRNA15 = OTTable(fastaread(rRNAtRNAPath), penaltyTableExactHomologyToExclude, ...  % Any region of exact homology equal to or greater than 15 nt will be removed
            'verbose', true, 'parallel', p);

        %% ------------------------------------------------------------------------
        % Build TRDesigner
        %%-------------------------------------------------------------------------
        %% Slice the transcriptome to the desired expression range to lower 
        % computational complexity by not calculating target regions for transcripts
        % that are not expressed within the desired range.

        PageBreak();

        % abundanceThreshold = 1e-2;
        abundanceThreshold = FPKMabundanceThreshold;

        fprintf(1, 'Slicing transcriptome based on expression level: >= %.2f FPKM\n', abundanceThreshold);
        fprintf(logFID, 'Slicing transcriptome based on expression level: >= %.2f FPKM\n', abundanceThreshold);

        % Find ids with abund >= 1e-2
        ids = transcriptome.ids;
        abund = transcriptome.GetAbundanceByID(ids);
        goodIDs = ids(abund>=abundanceThreshold);

        % Slice transcriptome
        slicedTranscriptome = transcriptome.Slice('geneID', goodIDs);

		fprintf(1, '%s - Sliced transcriptome\n', datestr(datetime));
        fprintf(logFID, '%s - Sliced transcriptome\n', datestr(datetime));
		
		% Having memory errors after this part.  
		% What is needed here?
		% slicedTranscriptome
		% OTrRNA15
		% specificityTable
		% isoSpecificityTable
		%
		% Needed in memory to pass to TRDesigner
		
		
		% NOT NEEDED
		% transcriptome
		%
		%
		% Can clear transcriptome object from memory at this point
		clear transcriptome;
		
        
		
        %% Create Target Region Designer object
        if ~exist(trDesignerPath)

            fprintf(logFID, '%s - Generating trDesigner object\n', datestr(datetime));

            trDesigner = TRDesigner('transcriptome', slicedTranscriptome, ...
                'OTTables', [OTrRNA15], ...
                'OTTableNames', {'rRNA'}, ...
                'specificityTable', specificityTable, ...
                'isoSpecificityTables', isoSpecificityTables, ...
                'parallel', p);
            trDesigner.Save(trDesignerPath);

            fprintf(logFID, '%s - trDesigner object saved to %s\n', datestr(datetime), trDesignerPath);
            
        elseif trDesignerPathOK
            fprintf(logFID, '%s - Using previous loaded trDesigner %s\n', datestr(datetime), trDesignerPath);
            fprintf(1, '%s - Using previous loaded trDesigner %s\n', datestr(datetime), trDesignerPath);
        else
            trDesigner = TRDesigner.Load(trDesignerPath);
            fprintf(logFID, '%s - trDesigner object loaded from %s\n', datestr(datetime), trDesignerPath);
        end

        %% ------------------------------------------------------------------------
        % Create target regions for a specific set of probe properties
        %%-------------------------------------------------------------------------
        
        % Can clear large objects from memory at this point
        clear specificityTable
        clear isoSpecificityTable
        clear slicedTranscriptome
%         clear OTrRNA15   - used to calc penalties for ncRNA later! 
        drawnow;
        
        
        % Apparent that if trRegionsPath object exists and is correct, nothing upstream needs to be loaded.
		% Problem is that TRDesigner object does not carry relevant info to validate loaded vs to-build object.
        % right now trusting that if path exists w/ correct parameters,
        % then rest is OK.  
		% If this could be done, then trDesignerPathOK check should be done first, then cascade backwards to transcriptomePathOK check,
		% then loading only what is not complete. 
        % Internal hash to make? 
        
        
        if ~exist(trRegionsPath)
            % Design target regions
        % 	targetRegions = trDesigner.DesignTargetRegions(...
        % 		'regionLength', 30, ...
        % 		'GC', [43 63]/100, ...
        % 		'Tm', [66 76], ...
        % 		'isoSpecificity', [0.75 1], ...
        % 		'specificity', [0.75 1], ...
        % 		'OTTables', {'rRNA', [0, 0]});

        %     	targetRegions = trDesigner.DesignTargetRegions(...
        % 		'regionLength', 30, ...
        %         'GC', [1 100]/100, ...
        % 		'Tm', [1 100], ...
        % 		'isoSpecificity', [0.75 1], ...
        % 		'specificity', [0.75 1], ...
        % 		'OTTables', {'rRNA', [0, 0]});

            fprintf(logFID, '%s - Designing target regions\n', datestr(datetime));

                targetRegions = trDesigner.DesignTargetRegions(...
                'regionLength', regionDesignParameters.regionLength, ...
                'GC', regionDesignParameters.GC, ...
                'Tm', regionDesignParameters.Tm, ...
                'isoSpecificity', regionDesignParameters.isoSpecificity, ...
                'specificity', regionDesignParameters.specificity, ...
                'monovalentSalt', regionDesignParameters.monovalentSaltConcentration, ...
                'probeConc', regionDesignParameters.probeConcentration, ...
                'threePrimeSpace', regionDesignParameters.probeSpacing, ...
                'OTTables', {'rRNA', [0, 0]},...
                'debugMode', debugMode);

                    % NOTE: The ranges above were determined empirically to strike 
                    % the proper balance between stringency (narrow ranges) and 
                    % sufficient probe numbers to target the desired genes. We 
                    % recommend scanning many different ranges to identify the 
                    % optimal for each application. 

            display(['... completed: ' datestr(now)]);

            % Save target regions
            targetRegions.Save(trRegionsPath);

            fprintf(logFID, '%s - Target regions saved to %s\n', datestr(datetime), trRegionsPath);

        else
            display(['Found: ' trRegionsPath]);
            targetRegions = TargetRegions.Load(trRegionsPath);

            fprintf(logFID, '%s - Target regions loaded from %s\n', datestr(datetime), trRegionsPath);
        end

        %write the genename
        genenames={targetRegions.geneName}';
        geneids={targetRegions.id}';
        numRegions={targetRegions.numRegions}';
        writetable(table(genenames,geneids,numRegions),[analysisSavePath 'targetRegions.csv'],'WriteRowNames',false)
        %% ------------------------------------------------------------------------
        % Step 2: Compile the library 
        %  The target regions designed above will be compiled into template
        %  molecules that can be used to build the desired probe library
        %%-------------------------------------------------------------------------

        
        %% ------------------------------------------------------------------------
        % Filter targetRegions 
        %%-------------------------------------------------------------------------
        % If steps above were not sufficiently restrictive for
        % targetRegions design
        % -- OR -- 
        % you want to implement a new way of filtering targetRegions, do so here.
        

        % Can filter to only those genes in the codebook, 
        % or to top abundance isoform of all genes
        switch targetRegionsFilter.geneIsoformListSource

            case 'codebook'

                % Pull target isoforms for each gene
                [geneIsoformList, ~] = LoadCodebook(codebookPath, 'verbose', true);


            case 'allGenes'
                % Do a filtering for ALL TargetRegions in the slicedTranscriptome
                allGeneNames = slicedTranscriptome.GetNames();
                geneIsoformList(length(allGeneNames)) = struct('name', '', ...
                                                               'id', '');
                for k = 1:length(allGeneNames)
                    idsHere = slicedTranscriptome.GetIDsByName(allGeneNames{k});
                    abundsHere = slicedTranscriptome.GetAbundanceByID(idsHere{1});
                    topAbundanceIsoform = idsHere{1}{abundsHere == max(abundsHere)};

                    geneIsoformList(k).name = allGeneNames{k};
                    % Include only most abundant isoform
                    geneIsoformList(k).id = topAbundanceIsoform;
                end
                
            case 'default'
                % No filtering to do, so this can be default

            otherwise
                error('geneListDefinedBy must be "codebook" or "allGenes"\n');
        end
        
        
        switch targetRegionsFilter.tRFilterMethod
            case ('default')
                % Do nothing.  Trust upstream filtering approach
                
            case ('parameter')
                % Filter by defined parameter
                
                targetRegions = findFilteredTargetRegions(targetRegions, ...
                                                          geneIsoformList, ...
                                                          numProbesPerGene, ...
                                                          logFID, ...
                                                          targetRegionsFilter.tRFilterField, ...
                                                          targetRegionsFilter.tRFilterParameters);
                                                      
            case ('relaxIsospecificity')
                % Return regions of genes that have min isospecificity to
                % yield at least minNumberOfProbes on each target isoform
                
                targetRegions = findExpandedIsospecificityTargetRegions(targetRegions, ...
                                        geneIsoformList, ...
                                        numProbesPerGene, ... 
                                        logFID);
                
            case ('commonRegions')
                % Return regions of genes that are common across isoforms
                % and each isoform carries 0 or at least minNumberOfProbes
                % targetRegions.
                
                targetRegions = findCommonTargetRegions(targetRegions, ...
                                                        geneIsoformList, ...
                                                        numProbesPerGene, ... 
                                                        logFID);
                
                
            otherwise
                error('obj.targetRegionsFilter.method incorrectly specified');
        end
        
        
        

        
        
        %% ------------------------------------------------------------------------
        % Load readouts, target regions, codewords, and selected genes
        %%-------------------------------------------------------------------------
        %% Load readouts and strip out the 3-letter readouts
        PageBreak();
        tic;
        display(['Loading: ' readoutPath]);

        fprintf(logFID, '%s - Loading readouts from %s\n', datestr(datetime), readoutPath);

        readouts = fastaread(readoutPath);
        display(['Found ' num2str(length(readouts)) ' oligos in ' num2str(toc) ' s']);

        fprintf(logFID, '%s - Found %d readouts in %d s\n', datestr(datetime), length(readouts), toc);

        %% Load codebook (which defines the readouts to use, 
        % the isoforms to use, and the barcodes assigned to them)

        fprintf(logFID, '%s - Loading codebook from %s\n', datestr(datetime), codebookPath);

        [codebook, codebookHeader] = LoadCodebook(codebookPath, 'verbose', true);
           % NOTE: A codebook should be defined before the library is constructed. 
           % NOTE: See the code_construction example script for instructions on how
           % to generate barcodes for different encoding schemes
           
           
        
        if specifyReadouts
        % ------------------------------------------------------------------------
        % Select readouts from readouts.fasta file
        % Put in order as specified from codebook
        %------------------------------------------------------------------------- 
        
        fprintf(logFID, '%s - Trimming readouts to specified list from codebook.\n', datestr(datetime), codebookPath);
        
        % Make ordered list of readouts given readouts.fasta, codebook
        
        readoutOrder = zeros(length(codebookHeader.bit_names), 1);
        for k = 1:length(readoutOrder)
           
            readoutFound = find(~cellfun(@isempty, strfind(cellstr(vertcat(readouts.Header)), codebookHeader.bit_names{k})), 1);
            
            if isempty(readoutFound)
               fprintf(logFID, '%s - Readout %s is in codebook but not in readouts.fasta file.\n', datestr(datetime), codebookHeader.bit_names{k});
               error('Readout %s is in codebook but not in readouts.fasta file.', codebookHeader.bit_names{k}); 
            else
                readoutOrder(k) = readoutFound;
            end

        end
        
        % Trim readouts struct down to the ones we want to use, in order.
        readouts = readouts(readoutOrder);
        
        fprintf(logFID, '%s - Retaining readouts .\n', datestr(datetime), codebookPath);
        
        else
            fprintf(logFID, '%s - Using default readout ordering.\n', datestr(datetime), codebookPath);
            
            % Nothing else necessary to do here.  Ordering in
            % readouts.fasta file will be used regardless of codebook
            % bit_names specification. 
        end
        
  
        % ------------------------------------------------------------------------
        % Select isoforms
        %-------------------------------------------------------------------------
        % Identify the isoforms to keep from those requested in the codebook
        finalIds = {codebook.id}; % Extract isoform ids from codebook

        % For codebooks that use a transcript spanning multiple IDs, there is a
        % string miss-match.  Replacing '_' with ',' in original codebook solves
        % this issue.

        %if ~strcmp(transcriptomeIDType, 'NCBI')
         %   finalIds = strrep(finalIds, '_', ',');
        %end

        finalGenes = {codebook.name}; % Extract gene common names from codebook
        barcodes = char({codebook.barcode}) == '1'; % Extract string barcodes and convert to logical matrix

        %check if gene and ids match the reference
        genes_from_ids=targetRegions(ismember({targetRegions.id}, finalIds));
        genenames=targetRegions(ismember({targetRegions.geneName}, finalGenes));
        if length(setdiff({genenames.geneName},{genes_from_ids.geneName}))>0
        warning([char(setdiff({genenames.geneName},{genes_from_ids.geneName})),' ',char(setdiff({genes_from_ids.id},{genenames.id})),' is different from the reference']);
        fprintf(logFID, ['WARNING! ', char(setdiff({genenames.geneName},{genes_from_ids.geneName})),' ',char(setdiff({genes_from_ids.id},{genenames.id})),' is different from the reference\n']);
        end

        if versionMatch
            finalTargetRegions = targetRegions(ismember({targetRegions.id}, finalIds)); % Extract only the desired target regions
        else
            % Strip version number from finalIDs and do fuzzy match

            % Approach with ~isempty(fTr) may drop all cases with no match to
            % reference transcriptome.  This might be an issue later.


            finalIdsStripped = cell(1, length(finalIds));
            for k = 1:length(finalIds)
                if contains(finalIds{k}, '.')
                    finalIdsStripped{k} = finalIds{k}(1:strfind(finalIds{k}, '.'));
                else
                    finalIdsStripped{k} = finalIds{k};
                end
            end

            i = 1;
            for k = 1:length(finalIds)
                fTr = targetRegions(contains({targetRegions.id}, finalIdsStripped{k}));
                if ~(isempty(fTr)) && ~(isempty(finalIdsStripped{k}))
                    if i == 1
                        finalTargetRegions = fTr;
                    else
                        finalTargetRegions(i) = fTr;
                    end        
                    i = i + 1;
                end

            end


        end
        

        

        %% ------------------------------------------------------------------------
        % Construct the library
        %%-------------------------------------------------------------------------
        %% Define common properties

        PageBreak();
        display(['Designing oligos for ' libraryName]);
        fprintf(logFID, '%s - Desgining oligos for %s\n', datestr(datetime), libraryName);
        display(['... ' num2str(numProbesPerGene) ' probes per gene']);
        fprintf(logFID, '%s - Using %d probes per gene\n', datestr(datetime), numProbesPerGene);

        %% Record the used readout sequences
        
        if exist(usedReadoutPath)
            warning(['Found ' usedReadoutPath]);
            delete(usedReadoutPath);
        end
        fastawrite(usedReadoutPath, readouts);
        PageBreak();
        display(['Wrote ' num2str(length(readouts)) ' readouts to ' usedReadoutPath]);
        fprintf(logFID, '%s - Wrote %d readouts to %s\n', datestr(datetime), length(readouts), usedReadoutPath);

        %% Build possible probes -- more than are needed are constructed to allow those with homology to rRNA/tRNA to be removed

        fprintf(logFID, '%s - Building possible oligomers. \n', datestr(datetime));
        
        if keepAllPossibleProbes
            keepGoingFlag = exist(possibleOligosPath, 'file') ~= 2;
        else
            keepGoingFlag = (exist(possibleOligosPath, 'file') ~= 2) && (exist(allOligosPath, 'file') == 2);
        end
        
        
        if keepGoingFlag
            oligos = [];

            allOligos = [];
            lastGene = '';
            
            Genes={};
            ProbeNumbers={};  
            %if keepAllPossibleProbes
             %   allHeaders = cell(sum(vertcat(finalTargetRegions.numRegions)), 1);
              %  allSeqs = cell(sum(vertcat(finalTargetRegions.numRegions)), 1);
               % seqCount = 1;
            %end

            for i=1:length(finalIds)
                % Save local gene
                localGeneName = finalGenes{i};

                if ~strcmp(localGeneName, lastGene)

                    % Display progress
                    PageBreak();
                    fprintf(logFID, '--------------------------------------------\n');
                    display(['Designing probes for ' libraryName ': ' localGeneName]);
                    fprintf(logFID, '%s - Designing probes for %s : %s. \n', datestr(datetime), libraryName, localGeneName);

                    % Determine the bits to include for each word
                    possibleReadouts = readouts(barcodes(i,:)==1); 

                    % Determine targetRegion sequences
                    tRegionPull = finalTargetRegions(strcmp({finalTargetRegions.geneName}, localGeneName));


                    if ~isempty(tRegionPull) % Check to see if there are no target regions--only used for blanks

                        seqs = {};
                        headers = {};
                        p = 1;

                        indsToKeepForReal = [];

                        % Block from here down originally implicitly asserts single gene -> single tRegion
                        % If multiple transcript IDs given same geneName this would break.

                        % Added in loop over number of tRegions to allow multiple
                        % tRegions for single geneName.
                        for tR = 1:length(tRegionPull)
                            tRegion = tRegionPull(tR);

                            fprintf(logFID, '-------------\n');
                            fprintf(logFID, '%s - Transcript : %s. \n', datestr(datetime), tRegion.id);

                            % Check for overlaps
                            if any(diff(tRegion.startPos) < 30)
                                fprintf(1, 'Overlap found for %s!', localGeneName);
                                fprintf(logFID, '%s - Overlap found for %s!\n', datestr(datetime), localGeneName);
                            end

                            % Build all possible oligos
                            for nR=1:tRegion.numRegions
                                
                                
                                % smELT localReadouts enforcement
                                if sum(barcodes(i,:)) == 1
                                    localReadouts(1) = possibleReadouts;
                                    localReadouts(2).Header = '';
                                    localReadouts(2).Sequence = '';
                                    localReadouts(3).Header = '';
                                    localReadouts(3).Sequence = '';
                                    longSide = 0;

                                % 2-hot readout testing, smELT style
                                elseif (sum(barcodes(i, :)) == 2) 
                                    localReadouts(1) = possibleReadouts(1);
                                    localReadouts(2).Header = '';
                                    localReadouts(2).Sequence = ''; 
                                    localReadouts(3) = possibleReadouts(2);
                                    longSide = 0;
                                else
                                    
                                    if readoutPermuteBySequence
                                        % Use permuteBySequence to generate
                                        % pseudo-random readout order
                                        [readoutOrder, whichPermute] = permuteBySequence(uint8(tRegion.sequence{nR}), sum(barcodes(i, :)), 3);
                                        localReadouts = possibleReadouts(readoutOrder); 
                                        % Use whichPermute to decide which
                                        % side has 2 readouts vs 1 readout.
                                        longSide = mod(whichPermute, 2) == 0;
                                    else
                                        % Create random orientation and selection of readouts
                                        localReadouts = possibleReadouts(randperm(length(possibleReadouts), 3));
                                        longSide = rand(1);
                                    end
                                end



                                if longSide > 0.5
                                    % Create header 
                                    headers{p} = [libraryName ' ' ...
                                        localReadouts(1).Header ' ' ...
                                        tRegion.geneName '__' ...
                                        tRegion.id '__' ...
                                        num2str(tRegion.startPos(nR)) '__' ...
                                        num2str(length(tRegion.sequence{nR})) '__' ...
                                        num2str(tRegion.GC(nR)) '__' ...
                                        num2str(tRegion.Tm(nR)) '__' ...
                                        num2str(tRegion.specificity(nR)) ' ' ...
                                        localReadouts(2).Header ' ' ...
                                        localReadouts(3).Header];

                                    % Create sequence
                                    seqs{p} = ['A ' seqrcomplement(localReadouts(1).Sequence) ' '...
                                            seqrcomplement(tRegion.sequence{nR}) ' A ' ...
                                            seqrcomplement(localReadouts(2).Sequence) ' ' ...
                                            seqrcomplement(localReadouts(3).Sequence)];
                                else
                                    % Create header 
                                    headers{p} = [libraryName ' ' ...
                                        localReadouts(1).Header ' ' ...
                                        localReadouts(2).Header ' ' ...
                                        tRegion.geneName '__' ...
                                        tRegion.id '__' ...
                                        num2str(tRegion.startPos(nR)) '__' ...
                                        num2str(length(tRegion.sequence{nR})) '__' ...
                                        num2str(tRegion.GC(nR)) '__' ...
                                        num2str(tRegion.Tm(nR)) '__' ...
                                        num2str(tRegion.specificity(nR)) ' ' ...
                                        localReadouts(3).Header];

                                    % Create sequence
                                    seqs{p} = ['A ' seqrcomplement(localReadouts(1).Sequence) ' ' ...
                                        seqrcomplement(localReadouts(2).Sequence) ' ' ...
                                        seqrcomplement(tRegion.sequence{nR}) ' A ' ...
                                           seqrcomplement(localReadouts(3).Sequence)];
                                end


                                p = p + 1;

                            end


                            display(['... constructed ' num2str(length(seqs)) ' possible probes']);
                            fprintf(logFID, '%s - Constructed %d possible probes\n', datestr(datetime), length(seqs));

                            seqsWOSpace = cellfun(@(x) x(~isspace(x)), seqs, 'UniformOutput', false);

                            % Identify penalties
                            hasrRNAPenalty = cellfun(@(x) sum(OTrRNA15.CalculatePenalty(seqrcomplement(x)))>0, seqsWOSpace);

                            % Select probes
                            indsToKeep = find(~hasrRNAPenalty);
                            indsToRemove = setdiff(1:length(seqs), indsToKeep);
                            display(['... removing ' num2str(length(indsToRemove)) ' probes']);
                            fprintf(logFID, '%s - Removing %d probes\n', datestr(datetime), length(indsToRemove));
                            for r=1:length(indsToRemove)
                                display(['...     ' headers{indsToRemove(r)}]);
                            end


            %                 display(indsToKeep)

                            %indsToKeepForReal = [indsToKeepForReal, indsToKeep];
                            indsToKeepForReal =  indsToKeep;
                        end
                           indsToKeepForAll = indsToKeepForReal(1:length(indsToKeepForReal));
                        indsToKeepForReal = indsToKeepForReal(randperm(length(indsToKeepForReal), min([length(indsToKeepForReal) numProbesPerGene])));
                        display(['... keeping ' num2str(length(indsToKeepForReal)) ' probes']);
                        fprintf(logFID, '%s - Retaining %d probes\n', datestr(datetime), length(indsToKeepForReal));
                        %write the genename
                        Genes=[Genes,localGeneName];
                        ProbeNumbers=[ProbeNumbers,length(indsToKeepForReal)];       

                        % Check on number
                        if length(indsToKeepForReal) < numProbesPerGene
                            warning(' ');
                            display(['Not enough probes for ' num2str(i) ': ' tRegion.geneName]);
                            fprintf(logFID, '%s - Not enough probes for %s!\n', datestr(datetime), tRegion.geneName);
                            
                            % If doubleHeadedsmELT flag set and there
                            % aren't enough probes to satisfy
                            % numProbesPerGene for this transcript, 
                            % append second readout sequence of same
                            % type to all probes on this gene
                            if doubleHeadedsmELT 
                                fprintf(logFID, '%s - Adding second readout to each probe for gene %s!\n', datestr(datetime), tRegion.geneName);
                                
%                                 assignin('base', 'headers', headers);
%                                 assignin('base', 'seqs', seqs);
%                                 assignin('base', 'indsToKeepForReal', indsToKeepForReal);
%                                 
                                for revInds = indsToKeepForReal
                                    
                                    headerLengthCorrectForsmELT = false;
                                    % Update header
                                    currHeader = headers{revInds};
                                    currHeadSplit = strsplit(currHeader, ' ');
                                    currHeadSplit(cell2mat(cellfun(@isempty, currHeadSplit, 'UniformOutput', false))) = [];
                                    % Only apply doubleHeaded to smELT
                                    % probes.  These have header form
                                    % libraryName Readout Gene_ID_XX__XX__XX
                                    % Append another Readout to end
                                    if length(currHeadSplit) == 3 % is right length for smELT (MERFISH has 5 components)
                                        headerLengthCorrectForsmELT = true;
                                        currHeadSplit{end+1} = currHeadSplit{2};
                                    end
                                    headers{revInds} = strjoin(currHeadSplit, ' ');
                                    
                                    % Update sequence
                                    currSeq = seqs{revInds};
                                    % smELT probes have form
                                    % A RdoutSeq GeneComplement A
                                    % Append RdoutSeq to end again
                                    currSeqSplit = strsplit(currSeq, ' ');
                                    if headerLengthCorrectForsmELT
                                        currSeqSplit{end + 1} = currSeqSplit{2};
                                    end
                                    seqs{revInds} = strjoin(currSeqSplit, ' ');
                                    
                                end
                                
                                
                                possibleReadouts = [possibleReadouts, possibleReadouts];
                            end
                            
                            
                            
                        end




                        % Save new oligos in oligos struct
                        % In this block, only those probes with valid
                        % indsToKeepForReal are retained.  If we want *ALL*
                        % possible probes, then see next block.
                        for s=1:length(indsToKeepForReal)
                            oligos(end+1).Header = headers{indsToKeepForReal(s)};
                            oligos(end).Sequence = seqs{indsToKeepForReal(s)};
                        end
                        for s=1:length(indsToKeepForAll)
                            allOligos(end+1).Header = headers{indsToKeepForAll(s)};
                            allOligos(end).Sequence = seqs{indsToKeepForAll(s)};
                        end

                        %if keepAllPossibleProbes
                            % Append all headers + seqs to cell
                         %   allHeaders(seqCount:(seqCount + length(headers) - 1)) = headers;
                          %  allSeqs(seqCount:(seqCount + length(seqs) - 1)) = seqs;
                           % seqCount = seqCount + length(seqs);
                        %end


                    else
                       fprintf(1, 'Empty!\n');
                    end

                    fprintf(logFID, '-------------\n');
                else
                    % Do nothing - already did this gene!
                end

                lastGene = localGeneName;

            end
            PageBreak();
            display(['Writing: ' possibleOligosPath]);
            fprintf(logFID, '%s - Writing %s!\n', datestr(datetime), possibleOligosPath);
            writeTimer = tic;
            fastawrite(possibleOligosPath, oligos);
            display(['... completed in ' num2str(toc(writeTimer))]);
            fprintf(logFID, '%s - Completed in %d s\n', datestr(datetime), toc(writeTimer));
            %%%%%Write out how many probes per gene as a table
            Genes=Genes'
            ProbeNumbers=ProbeNumbers'
            writetable(table(Genes,ProbeNumbers),[analysisSavePath 'probesPerGeneMergedCode.csv'],'WriteRowNames',false)
             
            if keepAllPossibleProbes
                if obj.debugMode
                    assignin('base', 'allHeaders', allHeaders);
                    assignin('base', 'allSeqs', allSeqs);
                end
                
                %allSeqs(cellfun(@isempty, allHeaders)) = [];
                %allHeaders(cellfun(@isempty, allHeaders)) = [];
                
                %allOligos = cell2struct([allHeaders, allSeqs], {'Header', 'Sequence'}, 2);
                
                PageBreak();
                fprintf(1, 'Writing: %s\n', allOligosPath);
                fprintf(logFID, '%s - Writing %s!\n', datestr(datetime), allOligosPath);
                writeTimer = tic;
                fastawrite(allOligosPath, allOligos);
                display(['... completed in ' num2str(toc(writeTimer))]);
                fprintf(logFID, '%s - Completed in %d s\n', datestr(datetime), toc(writeTimer));
            end

        else % End oligo design
            fprintf(logFID, '%s - Possible oligos file found.  Exiting. \n', datestr(datetime), toc(writeTimer));
            error('Found existing possible oligos file!');

        end

        %% Design primers -- removing those that have homology to the probes designed above
        display('Designing primers!');
        fprintf(logFID, '%s - Designing primers for %s\n', datestr(datetime), libraryName);

        p = parHold;

        
        if ~exist(primersPath)

            % Display progress
            PageBreak();
            display(['Designing primers for ' libraryName]);

            fprintf(logFID, '%s - Generating primers for %s\n', datestr(datetime), libraryName);

            % Build Off-Target Table for existing sequences and their reverse
            % complements
            seqRcomplement = cellfun(@(x)seqrcomplement(x(~isspace(x))), {oligos.Sequence}, 'UniformOutput', false);
            allSeqs = cellfun(@(x) x(~isspace(x)), {oligos.Sequence}, 'UniformOutput', false);
            allSeqs((end+1):(end+length(seqRcomplement))) = seqRcomplement;

            encodingProbeOTTable = OTTable(struct('Sequence', allSeqs), penaltyTableExactHomologyToExclude, 'verbose', true, ...
                'parallel', p);

            % Build primer designer
            prDesigner = PrimerDesigner('numPrimersToGenerate', primerDesignParameters.nPrimersToGenerate, ...
                'primerLength', primerDesignParameters.primerLength, ...
                'OTTables', encodingProbeOTTable, ...
                'OTTableNames', {'encoding'}, ...
                'parallel', p, ...
                'monovalentSalt', primerDesignParameters.primerMonovalentSaltConcentration, ...
                'primerConc', primerDesignParameters.primerConcentration);

            % Cut primers
            prDesigner.CutPrimers('Tm', primerDesignParameters.cutPrimers.Tm, ...
                'GC', primerDesignParameters.cutPrimers.GC, ...
                'OTTables', {'encoding', [0,0]});
            prDesigner.RemoveForbiddenSeqs();
            prDesigner.RemoveSelfCompPrimers('homologyMax', primerDesignParameters.cutPrimers.maxHomologySelf);
            prDesigner.RemoveHomologousPrimers('homologyMax', primerDesignParameters.cutPrimers.maxHomologyCross);

            % Write fasta file
            prDesigner.WriteFasta(primersPath);

            fprintf(logFID, '%s - Saving primers to %s\n', datestr(datetime), primersPath);
        else
            fprintf(logFID, '%s - Found existing primers file. Exiting.\n', datestr(datetime));
            error('Found existing primers!');
        end

        %% Add primers to the possible encoding probes designed above to generate template molecules
        primers = fastaread(primersPath);
        usedPrimers = primers(1:2);  % Select the first two of the valid primers generated above

        % Add primers to encoding probes
        PageBreak();
        display('Adding primers');
        fprintf(logFID, '%s - Concatenating primers to probes.\n', datestr(datetime));
        
        if ~exist(finalPrimersPath)
            % Record the used primers
            fastawrite(finalPrimersPath, usedPrimers);
            display(['Wrote: ' finalPrimersPath]);

            fprintf(logFID, '%s - Wrote primers to %s.\n', datestr(datetime), finalPrimersPath);

            % Build the final oligos
            finalOligos = [];
            for i=1:length(oligos)
                stringParts = strsplit(oligos(i).Header, ' ');
                name1 = strsplit(usedPrimers(1).Header, ' ');
                name1 = name1{1};
                name2 = strsplit(usedPrimers(2).Header, ' ');
                name2 = name2{1};
                finalOligos(i).Header = [stringParts{1} ' ' ...
                name1 ' '];
                for j=2:length(stringParts)
                    finalOligos(i).Header = [finalOligos(i).Header ...
                        stringParts{j} ' '];
                end
                finalOligos(i).Header = [finalOligos(i).Header name2];

                finalOligos(i).Sequence = [usedPrimers(1).Sequence ' ' ...
                oligos(i).Sequence ' ' ...
                seqrcomplement(usedPrimers(2).Sequence)];
            end
        else
            fprintf(logFID, '%s - Found existing final primers file. Exiting.\n', datestr(datetime));
            error('Found existing final primers path!')
        end

        %% Select the final template molecules -- Remove any template molecules with homology to noncoding RNAs and select only the desired number to keep
        PageBreak();
        display('Running final cross checks and building final fasta file');
        fprintf(logFID, '%s - Running final cross checks and building final fasta file\n', datestr(datetime));
        % Write final fasta
        
        if ~exist(oligosPath)
            % Screen against the original tables
            tic;
            display(['Searching oligos for homology']);
            fprintf(logFID, '%s - Searching oligos for homology\n', datestr(datetime));
            hasrRNAPenalty = cellfun(@(x) sum(OTrRNA15.CalculatePenalty(seqrcomplement(x(~isspace(x)))))>0, {finalOligos.Sequence});
            display(['... completed in ' num2str(toc) ' s']);
            fprintf(logFID, '%s - Completed in %d s\n', datestr(datetime), toc);

            indsToKeep = ~hasrRNAPenalty;
            display(['... found ' num2str(sum(~indsToKeep)) ' oligos to remove ']);
            fprintf(logFID, '%s - Found %d oligos to remove\n', datestr(datetime), sum(~indsToKeep));
            indsToRemove = find(~indsToKeep);
            for r=1:length(indsToRemove)
                display(['...     ' finalOligos(indsToRemove(r)).Header]);
            end

            % Remove bad oligos
            finalOligos = finalOligos(indsToKeep);

            % Write final oligos
            fastawrite(oligosPath, finalOligos);
            display(['Wrote: ' oligosPath]);
            fprintf(logFID, '%s - Wrote oligos to %s\n', datestr(datetime), oligosPath);
            display('Wrote oligos file to log');

        else
            fprintf(logFID, '%s - Found existing oligos file!\n', datestr(datetime), oligosPath);
            warning('Found existing oligos');
        end

        %% ------------------------------------------------------------------------
        % Archival and logging
        %%-------------------------------------------------------------------------
        %% Archive script
%         copyfile( [mfilename('fullpath'),'.m'],[analysisSavePath,mfilename,'.m']);
%         display('------------------------------------------------------------------');
%         display(['Copied analysis script to ' analysisSavePath,mfilename,'.m']);
        display('Execution complete!');
        fprintf(logFID, '%s - Execution completed\n', datestr(datetime), oligosPath);
        fclose(logFID);

    catch er    
        fprintf(logFID, '%s - Error!\n', datestr(datetime));
        fprintf(logFID, '%s - %s\n', datestr(datetime), er.identifier);
        fprintf(logFID, '%s - %s\n', datestr(datetime), er.message);
        fprintf(logFID, '%s - Exiting.\n', datestr(datetime));

        fclose(logFID);
        rethrow(er);
        return;
    end
    
    if debugMode
        % Debug mode transfers all variables in this workspace into base
        % workspace
        myVarList=who;

        for indVar = 1:length(myVarList)
            assignin('base',myVarList{indVar},eval(myVarList{indVar}))
        end
    end
    
    display('Script completed!');
    % Return to calling function
    return; 
    %     error('Stop script here!')






end

% Sub-functions
function badFilesFound = filePathCheck(pathsToCheck)   
% Check if generated file paths exist
% If yes/no, give print statement depending on class of file dependency
% If 'crashIfNone' not found or 'crashIfPresent' is found, then
% return badFilesFound = true
    badFilesFound = false;
    for k = 1:size(pathsToCheck, 1)

        if exist(pathsToCheck{k, 1}, 'file') == 2
            % File found at path
            switch pathsToCheck{k, 2}

                case 'load'
                    fprintf(1, 'File %s found.\nWill use version on disk.\n', pathsToCheck{k, 1})
                case 'crashIfNone'
                    fprintf(1, 'File %s found.\nProceeding.\n', pathsToCheck{k, 1})
                case 'delete'
                    fprintf(1, 'File %s found.\nProcessing will delete existing copy.\n', pathsToCheck{k, 1})
                case 'crashIfPresent'
                    fprintf(1, 'Warning! File %s found.\nDelete old copy and re-run!\n', pathsToCheck{k, 1})
                    badFilesFound = true;
            end

        else 
            switch pathsToCheck{k, 2}
                % File NOT found at path
                case 'load'
                    fprintf(1, 'File not %s found.\nWill generate new file.\n', pathsToCheck{k, 1})
                case 'crashIfNone'
                    fprintf(1, 'Warning! File not %s found!\nCheck path and re-run!\n', pathsToCheck{k, 1})
                    badFilesFound = true;
                case 'delete'
                    fprintf(1, 'File not %s found.\nWill generate new file.\n', pathsToCheck{k, 1})
                case 'crashIfPresent'
                    fprintf(1, 'File not %s found.\nProceeding.\n', pathsToCheck{k, 1})
            end
        end
    end
end
