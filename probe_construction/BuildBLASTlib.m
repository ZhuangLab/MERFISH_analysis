function dataBaseOut = BuildBLASTlib(fastaFile,varargin)
% ------------------------------------------------------------------------
% dataBaseOut = BuildBLASTlib(fastaFile,varargin)
% This function builds a blast library using the specified fasta file
%--------------------------------------------------------------------------
% Necessary Inputs
% fastaFile: String to a fasta file
%--------------------------------------------------------------------------
% Outputs
% dataBaseOut: Path to the created dataBase 
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Alistair Boettiger
% boettiger.alistair@gmail.com
% June 07 2014
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.



% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);

% BLAST version and path options
defaults(end+1,:) = {'blastPath', 'string', ''}; % Base tag for all images
defaults(end+1,:) = {'legacy', 'boolean', false};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1 || ~exist(fastaFile)
    error('matlabFunctions:invalidArguments', 'A valid fasta path is required.');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

%-------------------------------------------------------------------------
% Handle default blast paths
%-------------------------------------------------------------------------
if isempty(parameters.blastPath)
    if parameters.legacy
        parameters.blastPath = 'C:\"Program Files"\NCBI\LegacyBLAST\bin\';
    else
        parameters.blastPath = 'C:\"Program Files"\NCBI\blast-2.2.27+\bin\';
    end
end

%--------------------------------------------------------------------------
% Test Installation
%--------------------------------------------------------------------------
% Test BLAST installations
if ~parameters.legacy
    if ~exist([regexprep(parameters.blastPath,'"','') 'makeblastdb.exe'],'file')
        error(['Could not find NCBI blast+.  Please update blastPlusPath in ',mfilename]);
    end
else
    if ~exist([regexprep(parameters.blastPath,'"','') 'formatdb.exe'],'file')
        error(['Could not find NCBI blast.  Please update legacyBlastPath in ',mfilename]);
    end
end

%--------------------------------------------------------------------------
% Build database
%--------------------------------------------------------------------------
if parameters.legacy
    % Generate BLASTdatabase 
    command = [parameters.blastPath 'formatdb.exe' ...
        ' -i ' fastaFile ...
        ' -o T' ...  % -o T sets parse SeqID to true
        ' -p F']; %  -p F sets protein to false
    display('-----------------------------------------------------------------');
    display('Issuing:')
    display(['     ' command]);
    display('-----------------------------------------------------------------');
    system(command);
    dataBaseOut = regexprep(fastaFile,'.fasta','');
else
    dataBaseOut = regexprep(fastaFile,'.fasta','');

    % Build BLAST+ library
    command = [parameters.blastPath 'makeblastdb.exe' ...
        ' -dbtype "nucl"' ...
        ' -in ' fastaFile ...
        ' -parse_seqids' ...
        ' -out ' dataBaseOut]; % -o T sets parse SeqID to true -p F sets protein to false
    display('-----------------------------------------------------------------');
    display('Issuing:')
    display(['     ' command]);
    display('-----------------------------------------------------------------');
    system(command);
end