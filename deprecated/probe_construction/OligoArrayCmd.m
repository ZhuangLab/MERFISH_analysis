function [oligoArrayCommand, savePath] = OligoArrayCmd(varargin)
%--------------------------------------------------------------------------
% oligoArrayCommand = OligoArrayCmd('minTm',80,'probeLength',42,...
%   'crosshybeT',72,'sectstructT',72,'fastaName','myFasta');
%  Returns string formatted to send command to oligoarray.
%  The input file is expected to be in savePath
% 
%--------------------------------------------------------------------------
% Optional Parameters
%--------------------------------------------------------------------------
% GbMem = 4; % max memory in Gb an instance of oligoArray may use. 
% maxFragment = 1E3; 
% crosshybeT = 72;
% secstructT = 72;
% minTm = 80; 
% maxTm = 100;
% minPercentGC = 10;
% maxPercentGC = 90; 
% probeLength = 42;
% maskedSeq = '"GGGGGG;CCCC;TTTTTT;AAAA"'; 
% oligoArrayPath = 'C:\Users\Alistair\Documents\Research\Projects\Genomics\OligoArray\';
% oligoArrayExe = [oligoArrayPath '\OligoArray2.jar']; % Local oligarray;
% blastLib  = 'D:\Data\Genomics\DmelGenome\legacyBLASTlib\Dmel_Genome.fasta';
% savePath = ''; 
% saveDir = '';
% fastaName = '';
% inputPath = '';
% numParallel = 1;
% verbose = true;
%
%--------------------------------------------------------------------------
% Alistair Boettiger
% boettiger.alistair@gmail.com
% June 07 2014
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.

%% Default Parameters
GbMem = 3; % max memory in Gb an instance of oligoArray may use. 
maxFragment = 1E3; 
crosshybeT = 72;
secstructT = 72;
minTm = 80; 
maxTm = 100;
minPercentGC = 10;
maxPercentGC = 90; 
probeLength = 42;
probeLengthMax = [];
minProbeDist = [];
maskedSeq = '"GGGG;CCCC"'; % ;TTTT;AAAA
oligoArrayPath = 'C:\Users\Alistair\Documents\Research\Software\OligoArray\';
oligoArrayExe = [oligoArrayPath 'OligoArray2.jar']; % Local oligarray;
blastLib  = 'D:\Data\Genomics\DmelGenome\legacyBLASTlib\Dmel_Genome.fasta';
savePath = ''; 
saveDir = '';
fastaName = '';
inputPath = '';
numParallel = 1;
verbose = true;

%% Parse Optional Parameters

%--------------------------------------------------------------------------
% Parse variable input
%--------------------------------------------------------------------------
if nargin > 0
    if (mod(length(varargin), 2) ~= 0 ),
        error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
    end
    parameterCount = length(varargin)/2;
    for parameterIndex = 1:parameterCount,
        parameterName = varargin{parameterIndex*2 - 1};
        parameterValue = varargin{parameterIndex*2};
        switch parameterName
            case 'blastLib'
                blastLib = CheckParameter(parameterValue,'string','blastLib');
            case 'savePath'
                savePath = CheckParameter(parameterValue,'string','savePath');
            case 'saveDir'
                saveDir = CheckParameter(parameterValue,'string','saveDir');
            case 'inputPath'
                inputPath = CheckParameter(parameterValue,'string','inputPath');
            case 'oligoArrayExe'
                oligoArrayExe = CheckParameter(parameterValue,'string','oligoArrayExe');
            case 'probeLength'
                probeLength = CheckParameter(parameterValue,'positive','probeLength');
            case 'probeLengthMax'
                probeLengthMax = CheckParameter(parameterValue,'positive','probeLengthMax');
            case 'minTm'
                minTm = CheckParameter(parameterValue,'positive','minTm');
            case 'crosshybeT'
                crosshybeT = CheckParameter(parameterValue,'positive','crosshybeT');
            case 'secstructT'
                secstructT = CheckParameter(parameterValue,'positive','secstructT');
            case 'maxFragment'
                maxFragment = CheckParameter(parameterValue,'positive','maxFragment');
            case 'GbMem'
                GbMem = CheckParameter(parameterValue,'positive','GbMem');
            case 'maskedSeq'
                maskedSeq = CheckParameter(parameterValue,'string','maskedSeq');
            case 'minPercentGC'
                minPercentGC = CheckParameter(parameterValue,'positive','minPercentGC');
            case 'maxPercentGC'
                maxPercentGC = CheckParameter(parameterValue,'positive','maxPercentGC');
            case 'maxTm'
                maxTm = CheckParameter(parameterValue,'positive','maxTm');
            case 'minProbeDist'
                minProbeDist = CheckParameter(parameterValue,'positive','minProbeDist');
            case 'fastaName'
                fastaName = CheckParameter(parameterValue,'string','fastaName'); 
            case 'numParallel'
                numParallel= CheckParameter(parameterValue,'positive','numParallel'); 
            case 'verbose'
                verbose = CheckParameter(parameterValue,'boolean','verbose');
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

if isempty(probeLengthMax)
    probeLengthMax = probeLength;
end

if isempty(saveDir) && isempty(savePath)
    error(['Please specify a `saveDir` in which OligoArrayCmd can generate new folders ',...
        'or a `savePath` in which OligoArrayCmd can save data']);
end


if ~isempty(saveDir)
    saveFolder = strcat(datestr(date,'yy-mm-dd'),...
        '_Probes_',num2str(probeLength),'mers_t',...
        num2str(minTm),'Cx',num2str(crosshybeT),'Cs',num2str(secstructT),'C/');
    savePath = [saveDir,saveFolder]; % Location to save gene fasta files
    if ~exist(savePath,'dir')
        mkdir(savePath);
        if verbose
            disp(['Creating Folder ', savePath]); 
        end
    end
end
if isempty(inputPath)
    inputPath = savePath;
end

if isempty(minProbeDist)
    minProbeDist = probeLengthMax+1;
end

%% OligoArray2.1 Options

% note 'fastaInput' and 'genename' will get replaced by the appropriate
% gene and fasta file for which OligoArray2 is called prior to submitting
% the bsub command.  

oligoArrayCommand = ['java -Xmx',ceil(num2str(GbMem)),'g -jar ' oligoArrayExe ...    % Xmx6g is 6 Gb memory allowance. 
    ' -i ' inputPath 'genename' '.fasta' ...
    ' -d ' blastLib ...
    ' -o ' savePath 'genename' '_oligos','.txt' ... % found probes
    ' -r ' savePath 'genename' '_failed','.txt' ... % failed sequences
    ' -R ' savePath 'genename' '_log','.txt' ...    % log file
    ' -n ' num2str(round(maxFragment/probeLength)) ...              % The maximum number of oligos per entry   (30)
    ' -l ' num2str(probeLength) ...                          % Minimum oligo length      (32)
    ' -L ' num2str(probeLengthMax) ...                          % Maximum length            (32) 
    ' -D ' num2str(maxFragment) ...                        % Maximum distance between 5' of oligo and 3' of query
    ' -t ' num2str(minTm) ...                        % Minimum Tm
    ' -T ' num2str(maxTm) ...                         % Maximum Tm
    ' -s ' num2str(secstructT) ...                  % Temperature for secondary structure prediction
    ' -x ' num2str(crosshybeT) ...                  % Threshold for cross-hybridization
    ' -p ' num2str(minPercentGC) ...                          % Minimum GC content
    ' -P ' num2str(maxPercentGC) ...                          % Maximum GC content
    ' -m ',maskedSeq ...                            % Masked sequences   
    ' -g ' num2str(minProbeDist) ...           % Minimum distance between 5' end of oligos  %  (40) 
    ' -N ' num2str(numParallel)];                             % Number of parallel processors

if ~isempty(fastaName)
    fastaName = regexprep(fastaName,'.fasta','');
    oligoArrayCommand = regexprep(oligoArrayCommand,'genename',fastaName);
end

