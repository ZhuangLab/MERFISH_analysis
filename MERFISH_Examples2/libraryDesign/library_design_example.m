%% ------------------------------------------------------------------------
% Demonstrate the design of a MERFISH probe library
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% -------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.

%% ------------------------------------------------------------------------
% Setup the workspace
%%-------------------------------------------------------------------------
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

% % Human smELT
% basePath = [MERFISHAnalysisPath '\MERFISH_Examples2\']; % Base path where all required files can be found
% analysisSavePath = SetFigureSavePath([basePath 'libraryDesign\'], 'makeDir', true);
% rawTranscriptomeFasta = [basePath 'transcripts.fasta'];
% fpkmPath = [basePath 'isoformsMod.fpkm_tracking'];
% ncRNAPath = [basePath 'Homo_sapiens.GRCh38.ncrna.fa'];
% readoutPath = [basePath 'readouts.fasta'];
% codebookPath = [basePath 'codebookHumansmELTv2.csv'];

% % Mouse MERFISH + smELT, hypothalamus
basePath = [MERFISHAnalysisPath 'MERFISH_Examples2\Musmusculus\']; % Base path where all required files can be found
% analysisSavePath = SetFigureSavePath([basePath 'libraryDesign\Musmusculus\'], 'makeDir', true);

analysisSavePath = SetFigureSavePath([basePath 'libraryDesign\Musmusculusv4\'], 'makeDir', true);

rawTranscriptomeFasta = [basePath 'Mus_musculus.GRCm38.cdna.all.fa'];
fpkmPath = [basePath 'Mus_musculus_proxyRandomFPKM.fpkm_tracking'];
ncRNAPath = [basePath 'Mus_musculus.GRCm38.ncrna.fa'];
readoutPath = [basePath 'readoutsHypothalamus.fasta'];
codebookPath = [basePath 'codebookMusmusculusHypothalamus_v01.csv'];

%--------------------------------------
% Input parameters:
% Transcriptome:
transcriptomeHeaderType = 'ensembl';
transcriptomeIDType = 'ENSEMBL';

penaltyTableExactHomologyToExclude = 15;

FPKMabundanceThreshold = 0;

useUniformWeights = true;
isoSpecificityTable_lengthOfExactHomology = 17;

regionDesignParameters.regionLength = 30; % length of region (nt)
regionDesignParameters.GC = [0 1]; % Fractional content
regionDesignParameters.Tm = [0 100]; % deg C
regionDesignParameters.isoSpecificity = [0.75 1];
regionDesignParameters.specificity = [0.75 1];

numProbesPerGene = 92;
libraryName = ['L1E6'];

primerDesignParameters.nPrimersToGenerate = 1e3;
primerDesignParameters.primerLength = 20;
primerDesignParameters.cutPrimers.Tm = [70 72];
primerDesignParameters.cutPrimers.GC = [0.5 0.65];
primerDesignParameters.cutPrimers.maxHomologySelf = 6;
primerDesignParameters.cutPrimers.maxHomologyCross = 8;

%--------------------------------------

% codebookPath = [basePath, 'M22E1_codebook.csv'];

% transcriptomePath = 'C:\Users\ScanningLabAnalysis\Documents\MATLAB\MERFISH_analysis\MERFISH_Examples2\Musmusculus\150603_Transcriptome\transcriptomeC57Obj';


% Paths at which to save created objects

% analysisSavePath = SetFigureSavePath([basePath 'libraryDesign\Musmusculus'], 'makeDir', true);
% % 
% analysisSavePath = SetFigureSavePath('C:\Users\ScanningLabAnalysis\Documents\MATLAB\MERFISH_analysis\AIBSProbeDesign\HumanControlPanel1',...
%                                         'makeDir', true);




rRNAtRNAPath = [analysisSavePath 'rRNAtRNA.fa'];
transcriptomePath = [analysisSavePath 'transcriptomeObj'];
specificityTablePath = [analysisSavePath 'specificityTable'];
isoSpecificityTablePath = [analysisSavePath 'isoSpecificityTables'];
trDesignerPath = [analysisSavePath 'trDesigner'];
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

trRegionsPath =  [analysisSavePath trRegionsFileName];

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
    tic;
    ncRNAs = fastaread(ncRNAPath);
    display(['... completed in ' num2str(toc) ' s']);
    display(['Found ' num2str(length(ncRNAs)) ' sequences']);

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
    for i=1:length(biotypesToKeep)
        display(['   ' biotypesToKeep{i}]);
    end
    indsToKeep = ismember(biotypes, biotypesToKeep);

    rRNAtRNA = ncRNAs(indsToKeep);
    display(['Keeping ' num2str(length(rRNAtRNA)) ' ncRNAs']);

    if exist(rRNAtRNAPath)
        delete(rRNAtRNAPath);
    end
    
    % Save ncRNA in a fasta file
    fastawrite(rRNAtRNAPath, rRNAtRNA);
    display(['Wrote: ' rRNAtRNAPath]);
else
    % Load existing file if already created
    display(['Found and loading: ' rRNAtRNAPath]);
    tic;
    rRNAtRNA = fastaread(rRNAtRNAPath);
    display(['.... completed in ' num2str(toc) 's']);
    display(['Loaded ' num2str(length(rRNAtRNA)) ' sequences']);
end

%% ------------------------------------------------------------------------
% Build transcriptome object: This object collects information about
%   all transcripts in the transcriptome so as to facilitate acess to 
%   various properties.
%%-------------------------------------------------------------------------
%% Build transcriptome object
if ~exist(transcriptomePath)
    % Build transcriptome using existing abundance data
    transcriptome = Transcriptome(rawTranscriptomeFasta, ...
        'abundPath', fpkmPath, ...
        'verbose', true, ...
        'headerType', transcriptomeHeaderType, ...
        'IDType', transcriptomeIDType);
    
    transcriptome.Save(transcriptomePath);
else
    % Load transcriptome if it already exists
    transcriptome = Transcriptome.Load(transcriptomePath);
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
    display(['Started on ' datestr(now)]);

    clear isoSpecificityTables
    
    for i=1:length(names) % Loop over all gene names -- RNAs that share the same gene name are considered isoforms
        % Generate a local transcriptome object that contains transcripts for a single gene
        localTranscriptome = transcriptome.Slice('geneID', idsByName{i});
        
        if mod(i, 100) == 0
            fprintf(1, 'On gene %d of %d. First transcript %s\n', i, length(names), idsByName{i}{1});
        end
        
        if useUniformWeights

            % Generate a OTTable for isoforms for the given gene
            isoSpecificityTables(i) = OTTable(localTranscriptome, ...
                isoSpecificityTable_lengthOfExactHomology, ...  % lengthOfExactHomology is the length of exact homology used to calculate penalties
                'verbose', false, ...
                'transferAbund', true, ...
                'weights', []);
        
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
        end
    end
    display(['.... completed in ' num2str(toc) ' s']);
    % Save tables
    isoSpecificityTables.Save(isoSpecificityTablePath);
else
    % Load tables if they already exist
    isoSpecificityTables = OTTable.Load(isoSpecificityTablePath, ...
        'verbose', true, ...
        'mapType', @OTMap2);
end

%% Build total specificity table --- this table contains a penalty associated
%   with all possible sequences in the transcriptome
if ~exist(specificityTablePath)
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
else
    % Load table if it exists
    specificityTable = OTTable.Load(specificityTablePath, ...
        'verbose', true, ...
        'mapType', @OTMap2);
end

%% ------------------------------------------------------------------------
% Configure parallel.Pool
%%-------------------------------------------------------------------------
%% Create parallel pool... speeds up the construction of the TRDesigner and the construction of libraries
if isempty(gcp('nocreate'))
    p = parpool(5);  % Insert a number here appropriate to the used computational resources
else
    p = gcp;
end

parHold = p;

%% ------------------------------------------------------------------------
% Build Penalty Tables
%%-------------------------------------------------------------------------
%% Build rRNA/tRNA penalty table
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



% Find ids with abund >= 1e-2
ids = transcriptome.ids;
abund = transcriptome.GetAbundanceByID(ids);
goodIDs = ids(abund>=abundanceThreshold);

% Slice transcriptome
slicedTranscriptome = transcriptome.Slice('geneID', goodIDs);

%% Create Target Region Designer object
if ~exist(trDesignerPath)
    trDesigner = TRDesigner('transcriptome', slicedTranscriptome, ...
        'OTTables', [OTrRNA15], ...
        'OTTableNames', {'rRNA'}, ...
        'specificityTable', specificityTable, ...
        'isoSpecificityTables', isoSpecificityTables, ...
        'parallel', p);
    trDesigner.Save(trDesignerPath);
else
    trDesigner = TRDesigner.Load(trDesignerPath);
end

%% ------------------------------------------------------------------------
% Create target regions for a specific set of probe properties
%%-------------------------------------------------------------------------
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
    
        targetRegions = trDesigner.DesignTargetRegions(...
		'regionLength', regionDesignParameters.regionLength, ...
        'GC', regionDesignParameters.GC, ...
		'Tm', regionDesignParameters.Tm, ...
		'isoSpecificity', regionDesignParameters.isoSpecificity, ...
		'specificity', regionDesignParameters.specificity, ...
		'OTTables', {'rRNA', [0, 0]});
    
            % NOTE: The ranges above were determined empirically to strike 
            % the proper balance between stringency (narrow ranges) and 
            % sufficient probe numbers to target the desired genes. We 
            % recommend scanning many different ranges to identify the 
            % optimal for each application. 
            
	display(['... completed: ' datestr(now)]);

	% Save target regions
	targetRegions.Save(trRegionsPath);

else
    display(['Found: ' trRegionsPath]);
    targetRegions = TargetRegions.Load(trRegionsPath);
end

%% ------------------------------------------------------------------------
% Step 2: Compile the library 
%  The target regions designed above will be compiled into template
%  molecules that can be used to build the desired probe library
%%-------------------------------------------------------------------------

%% ------------------------------------------------------------------------
% Load readouts, target regions, codewords, and selected genes
%%-------------------------------------------------------------------------
%% Load readouts and strip out the 3-letter readouts
PageBreak();
tic;
display(['Loading: ' readoutPath]);
readouts = fastaread(readoutPath);
display(['Found ' num2str(length(readouts)) ' oligos in ' num2str(toc) ' s']);

%% Load codebook (which defines the readouts to use, 
% the isoforms to use, and the barcodes assigned to them)
codebook = LoadCodebook(codebookPath);
   % NOTE: A codebook should be defined before the library is constructed. 
   % NOTE: See the code_construction example script for instructions on how
   % to generate barcodes for different encoding schemes

%% ------------------------------------------------------------------------
% Select isoforms
%%-------------------------------------------------------------------------
%% Identify the isoforms to keep from those requested in the codebook
finalIds = {codebook.id}; % Extract isoform ids from codebook

% For codebooks that use a transcript spanning multiple IDs, there is a
% string miss-match.  Replacing '_' with ',' in original codebook solves
% this issue.

finalIds = strrep(finalIds, '_', ',');

finalGenes = {codebook.name}; % Extract gene common names from codebook
barcodes = char({codebook.barcode}) == '1'; % Extract string barcodes and convert to logical matrix

% Indicate if this is allowed to ignore sequence version number in matching
versionMatch = false;
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
display(['... ' num2str(numProbesPerGene) ' probes per gene']);

%% Record the used readout sequences
usedReadoutPath = [analysisSavePath libraryName '_used_readouts.fasta'];
if exist(usedReadoutPath)
    warning(['Found ' usedReadoutPath]);
    delete(usedReadoutPath);
end
fastawrite(usedReadoutPath, readouts);
PageBreak();
display(['Wrote ' num2str(length(readouts)) ' readouts to ' usedReadoutPath]);

%% Build possible probes -- more than are needed are constructed to allow those with homology to rRNA/tRNA to be removed
oligosPath = [analysisSavePath libraryName '_possible_oligos.fasta'];
if ~exist(oligosPath)
    oligos = [];
    
    lastGene = '';
    
    for i=1:length(finalIds)
        % Save local gene
        localGeneName = finalGenes{i};
        
        if ~strcmp(localGeneName, lastGene)

            % Display progress
            PageBreak();
            display(['Designing probes for ' libraryName ': ' localGeneName]);

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

                    % Check for overlaps
                    if any(diff(tRegion.startPos) < 30)
                        fprintf(1, 'Overlap found for %s!', localGeneName);
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
                        else
                            % Create random orientation and selection of readouts
                            localReadouts = possibleReadouts(randperm(length(possibleReadouts), 3));
                        end

                        if rand(1) > 0.5
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

                    seqsWOSpace = cellfun(@(x) x(~isspace(x)), seqs, 'UniformOutput', false);

                    % Identify penalties
                    hasrRNAPenalty = cellfun(@(x) sum(OTrRNA15.CalculatePenalty(seqrcomplement(x)))>0, seqsWOSpace);

                    % Select probes
                    indsToKeep = find(~hasrRNAPenalty);
                    indsToRemove = setdiff(1:length(seqs), indsToKeep);
                    display(['... removing ' num2str(length(indsToRemove)) ' probes']);
                    for r=1:length(indsToRemove)
                        display(['...     ' headers{indsToRemove(r)}]);
                    end

    %                 display(indsToKeep)

                    indsToKeepForReal = [indsToKeepForReal, indsToKeep];

                end

                indsToKeepForReal = indsToKeepForReal(randperm(length(indsToKeepForReal), min([length(indsToKeepForReal) numProbesPerGene])));
                display(['... keeping ' num2str(length(indsToKeepForReal)) ' probes']);

                % Check on number
                if length(indsToKeep) < numProbesPerGene
                    warning(' ');
                    display(['Not enough probes for ' num2str(i) ': ' tRegion.geneName]);
                end




                % Save new oligos in oligos struct
                for s=1:length(indsToKeepForReal)
                    oligos(end+1).Header = headers{indsToKeepForReal(s)};
                    oligos(end).Sequence = seqs{indsToKeepForReal(s)};
                end


            else
               fprintf(1, 'Empty!\n');
            end
            
        else
            % Do nothing - already did this gene!
        end
        
        lastGene = localGeneName;
        
    end
    PageBreak();
    display(['Writing: ' oligosPath]);
    writeTimer = tic;
    fastawrite(oligosPath, oligos);
    display(['... completed in ' num2str(toc(writeTimer))]);

else % End oligo design
    error('Found existing possible oligos file!');
end

%% Design primers -- removing those that have homology to the probes designed above
display('Designing primers!');

p = parHold;

primersPath = [analysisSavePath libraryName '_possible_primers.fasta'];
if ~exist(primersPath)

    % Display progress
    PageBreak();
    display(['Designing primers for ' libraryName]);

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
        'parallel', p);

    % Cut primers
    prDesigner.CutPrimers('Tm', primerDesignParameters.cutPrimers.Tm, ...
        'GC', primerDesignParameters.cutPrimers.GC, ...
        'OTTables', {'encoding', [0,0]});
    prDesigner.RemoveForbiddenSeqs();
    prDesigner.RemoveSelfCompPrimers('homologyMax', primerDesignParameters.cutPrimers.maxHomologySelf);
    prDesigner.RemoveHomologousPrimers('homologyMax', primerDesignParameters.cutPrimers.maxHomologyCross);

    % Write fasta file
    prDesigner.WriteFasta(primersPath);
else
    error('Found existing primers!');
end

%% Add primers to the possible encoding probes designed above to generate template molecules
primers = fastaread(primersPath);
usedPrimers = primers(1:2);  % Select the first two of the valid primers generated above

% Add primers to encoding probes
PageBreak();
display('Adding primers');
finalPrimersPath = [analysisSavePath libraryName '_primers.fasta'];
if ~exist(finalPrimersPath)
    % Record the used primers
    fastawrite(finalPrimersPath, usedPrimers);
    display(['Wrote: ' finalPrimersPath]);

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
    error('Found existing final primers path!')
end

%% Select the final template molecules -- Remove any template molecules with homology to noncoding RNAs and select only the desired number to keep
PageBreak();
display('Running final cross checks and building final fasta file');
% Write final fasta
oligosPath = [analysisSavePath libraryName '_oligos.fasta'];
if ~exist(oligosPath)
    % Screen against the original tables
    tic;
    display(['Searching oligos for homology']);
    hasrRNAPenalty = cellfun(@(x) sum(OTrRNA15.CalculatePenalty(seqrcomplement(x(~isspace(x)))))>0, {finalOligos.Sequence});
    display(['... completed in ' num2str(toc) ' s']);

    indsToKeep = ~hasrRNAPenalty;
    display(['... found ' num2str(sum(~indsToKeep)) ' oligos to remove ']);
    indsToRemove = find(~indsToKeep);
    for r=1:length(indsToRemove)
        display(['...     ' finalOligos(indsToRemove(r)).Header]);
    end

    % Remove bad oligos
    finalOligos = finalOligos(indsToKeep);

    % Write final oligos
    fastawrite(oligosPath, finalOligos);
    display(['Wrote: ' oligosPath]);

else
    warning('Found existing oligos');
end

%% ------------------------------------------------------------------------
% Archival and logging
%%-------------------------------------------------------------------------
%% Archive script
copyfile( [mfilename('fullpath'),'.m'],[analysisSavePath,mfilename,'.m']);
display('------------------------------------------------------------------');
display(['Copied analysis script to ' analysisSavePath,mfilename,'.m']);

%% Beep
error('Stop script here!')
