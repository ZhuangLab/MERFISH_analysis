% probeDesign class for MERFISHProbeDesign.m input

classdef probeDesign < matlab.mixin.SetGet
    
    properties
        % Set parameter validation and defaults
        MERFISHAnalysisPath char = 'C:\Users\ScanningLabAnalysis\Documents\MATLAB\MERFISH_analysis\';
        basePath char = 'D:\Data\MERFISH\Musmusculus\'; % Base path where all required files can be found

        libraryName char = '';
        
        species char {mustBeMember(species, {'Mus musculus', 'Homo sapiens', 'human', 'mouse', ''})} = 'Mus musculus';

        rawTranscriptomeFasta char = 'Mus_musculus.GRCm38.cdna.all.fa';
        fpkmPath char = 'Mus_musculus_proxyRandomFPKM.fpkm_tracking';
        ncRNAPath char = 'Mus_musculus.GRCm38.ncrna.fa';

        readoutPath char = 'D:\Data\MERFISH\AIBSMERFISH_allReadouts.fasta';
        codebookPath char = '';

        transcriptomeHeaderType char {mustBeMember(transcriptomeHeaderType,{'ensembl','refseq', 'cufflinks', ''})} = 'ensembl'; % 'cufflinks', 'ensembl', or 'refseq'
        transcriptomeIDType char {mustBeMember(transcriptomeIDType,{'ENSEMBL','NCBI', ''})} = 'ENSEMBL'; % 'ENSEMBL' or 'NCBI'
        %--------------------------------------
        % Input parameters:
        % Transcriptome:
        penaltyTableExactHomologyToExclude {mustBeNumeric, mustBeNonnegative} = 15;

        FPKMabundanceThreshold double {mustBeNumeric, mustBeNonnegative} = 0;

        useUniformWeights {mustBeNumericOrLogical} = true;
        isoSpecificityTable_lengthOfExactHomology {mustBeNumeric, mustBeNonnegative} = 17;

        regionLength {mustBeNumeric} = 30; % length of region (nt)
        regionGC = [0, 1]; % Fractional content
        regionTm = [70, 100]; % deg C
        isoSpecificity = [0, 1]; % Within a gene
        specificity = [0, 1]; % Across genes
        monovalentSaltConcentration {mustBeNumeric, mustBeNonnegative} = 0.3; % mol/L
        probeConcentration {mustBeNumeric, mustBeNonnegative} = 5e-9; % mol/L
        probeSpacing {mustBeNumeric} = 3; % nt of gap between probes - set to negative to allow overlap

        numProbesPerGene {mustBeNumeric, mustBeNonnegative} = 48;


        nPrimersToGenerate {mustBeNumeric, mustBeNonnegative} = 1e3;
        primerLength {mustBeNumeric, mustBeNonnegative} = 20;
        primerMonovalentSaltConcentration{mustBeNumeric, mustBeNonnegative} = 0.3;
        primerConcentration {mustBeNumeric, mustBeNonnegative} = 0.5e-6;
        cutPrimersTm = [70, 72];
        cutPrimersGC = [0.5, 0.65];
        cutPrimersMaxHomologySelf {mustBeNumeric, mustBeNonnegative} = 6;
        cutPrimersMaxHomologyCross {mustBeNumeric, mustBeNonnegative} = 8;

        % Indicate if this is allowed to ignore sequence version number in matching
        versionMatch {mustBeNumericOrLogical} = false;
        
        doubleHeadedsmELT {mustBeNumericOrLogical} = false;
        
        keepAllPossibleProbes {mustBeNumericOrLogical} = true;
        
        %---------------------------------------------------------------
        % Generated paths that could be reused if existing, things match
        % 
        rRNAtRNAPath char = ''; % Depends on ncRNA FASTA file
        transcriptomePath char = ''; % raw transcriptome and fpkm file
        specificityTablePath char = ''; % specificityTablePath
        isoSpecificityTablePath char = ''; % transcriptome, isoSpecificityTable_lengthOfExactHomology
        trDesignerPath char = ''; % rRNAtRNAPath, penaltyTableExactHomologyToExclude, FPKMabundanceThreshold, transcriptome, specificityTable, isoSpecificityTables

    end
    
    methods
        % Constructor
        % When called with 3 arguments, initialize species paths from
        % defaults
        % Otherwise when no arguments, initialize w/o path defaults
        function pd = probeDesign(libraryName, species, codebookPath)
            if nargin > 0
                pd.libraryName = libraryName;
                pd.species = species;
                pd.codebookPath = codebookPath;
                
                if (strcmpi((pd.species), 'mus musculus')) || strcmpi((pd.species), 'mouse')
                    
                    pd.basePath = 'D:\Data\MERFISH\Musmusculus\';
                    pd.rawTranscriptomeFasta = 'Mus_musculus.GRCm38.cdna.all.fa';
                    pd.fpkmPath = 'Mus_musculus_proxyRandomFPKM.fpkm_tracking';
                    pd.ncRNAPath = 'Mus_musculus.GRCm38.ncrna.fa';
                    
                    pd.transcriptomeHeaderType = 'ensembl';
                    pd.transcriptomeIDType = 'ENSEMBL';
                    
                elseif (strcmpi((pd.species), 'homo sapiens')) || strcmpi((pd.species), 'human')
                    
                    pd.basePath = 'D:\Data\MERFISH\Homosapiens\';
                    pd.rawTranscriptomeFasta = 'Homo_sapiens_GRCh38_latest_rna.fna';
                    pd.fpkmPath = 'Homo_sapiens_proxyRandomFPKM.fpkm_tracking';
                    pd.ncRNAPath = 'Homo_sapiens.GRCh38.ncrna.fa';
                    
                    pd.transcriptomeHeaderType = 'refseq';
                    pd.transcriptomeIDType = 'NCBI';
                    
                else
                    
                    error('Species %s not supported', pd.species);
                    
                end
              
            else
                pd.libraryName = '';
                pd.species = '';
                pd.codebookPath = '';
                
                pd.rawTranscriptomeFasta = '';
                pd.fpkmPath = '';
                pd.ncRNAPath = '';

                pd.transcriptomeHeaderType = '';
                pd.transcriptomeIDType = '';
                
            end
        end
            
        % Make parameters match those from log file 
        function matchLogFile(obj, filePath)
            
            % Read in file header
            % All in property : value pairs
            % Header ends at third '--------...' line
            
            allowedSplits = 3;
            
            fID = fopen(filePath, 'r');
            libName = strsplit(fgetl(fID), ' ');
            obj.libraryName = libName{1};
            
            splitsFound = 0;
            primerBlockFound = false;
            while (splitsFound < allowedSplits)
                lineHere = fgetl(fID);
                
                if (lineHere > 0)
                    ret = parseLogLine(lineHere);
                else
                    continue;
                end
                
%                 display(ret{1})
                    
                switch ret{1}
                    
                    
                    case ''
                        % Pass
                        
                    case '-'

                        splitsFound = splitsFound + 1;
%                         disp('Found split')

                    case {'basePath', ...
                          'rawTranscriptomeFasta', ...
                          'fpkmPath', ...
                          'ncRNAPath',...
                          'readoutPath', ...
                          'codebookPath', ...
                          'transcriptomeHeaderType', ...
                          'transcriptomeIDType', ...
                          'penaltyTableExactHomologyToExclude', ...
                          'FPKMabundanceThreshold', ...
                          'useUniformWeights', ...
                          'isoSpecificityTable_lengthOfExactHomology', ...
                          'regionLength', ...
                          'monovalentSaltConcentration', ...
                          'probeConcentration', ...
                          'probeSpacing',...
                          'numProbesPerGene', ...
                          'nPrimersToGenerate', ...
                          'primerLength', ...
                          'rRNAtRNAPath', ...
                          'transcriptomePath', ...
                          'specificityTablePath', ...
                          'isoSpecificityTablePath', ...
                          'trDesignerPath'}
                        
                        if all(isstrprop(ret{2}, 'digit') | (ret{2} == '.') | (ret{2} == 'e') | (ret{2} == '-'))
                            obj.(ret{1}) = str2double(ret{2});
                        else
                            obj.(ret{1}) = ret{2};
                        end

                    case {'GC', 'Tm', 'isoSpecificity', 'specificity'}
                        
                        % To do - convert from strings to pairs of doubles
                        
                        
                        % could be readout or primer
                        if primerBlockFound
                            if strcmp(ret{1}, 'GC')
                                obj.cutPrimersGC = str2double(strsplit(ret{2}(2:(end-1)), ','));
                            elseif strcmp(ret{1}, 'Tm')
                                obj.cutPrimersTm = str2double(strsplit(ret{2}(2:(end-1)), ','));
                            else
                                error('Bad value at primerBlockFound as true in switch-case');
                            end
                        else 
                            if strcmp(ret{1}, 'GC')
                                obj.regionGC = str2double(strsplit(ret{2}(2:(end-1)), ','));
                            elseif strcmp(ret{1}, 'Tm')
                                obj.regionTm = str2double(strsplit(ret{2}(2:(end-1)), ','));
                            elseif strcmp(ret{1}, 'isoSpecificity')
                                obj.isoSpecificity = str2double(strsplit(ret{2}(2:(end-1)), ','));
                            elseif strcmp(ret{1}, 'specificity')
                                obj.specificity = str2double(strsplit(ret{2}(2:(end-1)), ','));
                            else
                                error('Bad value at primerBlockFound as false in switch-case');
                            end
                        end
                        

                    case 'primerDesignParameters'
                        % Flag for GC, Tm going to region or primer parameters
                        primerBlockFound = true;

                    case 'maxHomologySelf'
                        obj.cutPrimersMaxHomologySelf = str2double(ret{2});

                    case 'maxHomologyCross'
                        obj.cutPrimersMaxHomologyCross = str2double(ret{2});

                    case 'transcriptome-codebook_version_match_enforced'
                        obj.versionMatch = (ret{2} == '1');

                    case 'doubleHeadedsmELT'
                        assignin('base', 'line2', ret)
                        obj.doubleHeadedsmELT = (ret{2} == '1');
                        
                end
            end
            fclose(fID);
            
            
            function ret = parseLogLine(nextLine)
                splt = strsplit(nextLine, ' : ');
                if length(splt) == 2
                    % property : value pair
                    ret  = splt;
                else 
                    if nextLine(1) == '-'
                        ret = {'-'};
                    elseif strcmp(nextLine, 'primerDesignParameters')
                        ret = {'primerDesignParameters'};
                    else
                        ret = {''};
                    end
                end
            end
            
            
        end

        
%     end
%     
%     methods (Access = private)
        % Check that class has required components, then build library
        % using MERFISHProbeDesign.m function
        
        function buildLibrary(obj, varargin)
            
            % Check that required properties are all set
            
            props = properties(obj);
            validProp = zeros(length(props), 1);
            for iprop = 1:length(props)
              thisprop = props{iprop};
              
              propValue = obj.(thisprop);
              
              switch thisprop
                  case {'rRNAtRNAPath', ...
                         'transcriptomePath', ...
                         'specificityTablePath', ...
                         'isoSpecificityTablePath', ...
                         'trDesignerPath'}
                     % These properties can be empty at start.
                     
                     validProp(iprop) = 0;
                     
                  otherwise
                      
                    validProp(iprop) = isempty(propValue);
              end
              
            end
            
            if all(validProp == 0)
                fprintf(1, 'Building library\n');
                if (nargin == 2) && strcmp(varargin{1}, 'test')
                    returnTest(obj);
                else
                     MERFISHProbeDesign(obj);
                end
            else
                fprintf(1, 'Empty parameters in probeDesign object.\n');
                fprintf(1, 'Correct and restart buildLibrary method.\n');
            end
        end % buildLibrary
    end % methods
end % classDef