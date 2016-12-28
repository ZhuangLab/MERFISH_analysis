function [oligoSeq, oligoName, codebookTable] = AssembleProbes(ProbeData,varargin)
% [oligoSeq, oligoName, codebookTable] = AssembleProbes(ProbeData,varargin)
%--------------------------------------------------------------------------
% Outputs
%--------------------------------------------------------------------------
% AssembleProbes(ProbeData)
% 1. E1_codebook  - fasta file recording the codeword assigned to each
%                   gene, and the names of oligos assigned to each bit.
% 2. E1_oligos   - fasta file containing all the oligos in the library. the
%                   header records the experiment number, the and lists the
%                   binding regions included in the probe as a
%                   space-delimited entry.  The sequence is also space
%                   delimited.  This fasta file may be concatinated with
%                   other fasta files which employ orthoganol index
%                   primers.
%    Example probe output:
%    > E1 primer_002 B02 ABCA2____107 B03 primer_003 probe005
%    CGCGGGCTATATGCGAACCG  GCTATCGTTCGTTCGAGGCCAGAGCATTCG TGGCCCAAGAAGGAGACCACCTGGTTGCGG 
%    CCCATGATCGTCCGATCTGGTCGGATTTGT GCGTTGTATGCCCTCCACGC
%    is the 5th probe in the library, this one binds starting at bp 107 of the gene
%    ABCA2.  This probe will be bound by readout-sequence 2 and readout-
%    sequence 3.  It is part of experiment 1, which can be amplified using 
%    primer sequences 2 and 3 (the 20 mers at each end).  
%
% Optionally, return cell arrays of the oligo sequences, the oligo names,
% and a table recording the codebook used: 
% [oligoSeq, oligoName, codebookTable] = AssembleProbes(ProbeData);
% 
%--------------------------------------------------------------------------
% Required Inputs
%--------------------------------------------------------------------------
% ProbeData - a matlab structure, which contains the name of the gene and 
%          all the unique probes found by OligoArray. 
% Either a data folder containing the following fasta files or the paths to
% the following fasta files must also be specified.  See Optional inputs
% for more information:
% 
% 1. TargetGeneSeqs.fasta - FASTA file containing all genes of interest.  
%    This may be the whole transcriptome or just a subset of genes. 
% 2. PrimerSeqs.fasta - FASTA file containing orthoganol sequences to use
%    and index primers. These will be used in the inital PCR reaction in
%    the probe synthesis protocol and allow multiple libraries to be
%    constructed from the same oligo-pool.
% 3. ReadoutSeqs.fasta - FASTA file containing 16 orthoganol sequence to
%    use in the readout hybridizations.
% 4. BLASTedRandomSeqs.fasta - FASTA file containing sequences which have
%    little or no homology to the human genome. Used to construct negative
%    controls for the library.
% 
%--------------------------------------------------------------------------
% Optional Inputs
%--------------------------------------------------------------------------
%  parameter flag, default value, description 
%    'numCntrls',5,... % number of random sequences to design for negative controls 
%    'numBlanks',5,...  % minimum number of codewords left blank as controls
%    'numOligos',200,... % total oligos per probe (will be divided across the readout sequences).  
%    'subLibNum',1,... % starting sublibrary (increase to avoid using the same index primers as for a different experiment); 
%    'numTotalBits',16,...   % Total number of bits to read out (must have at least this number of secondary sequences 
%    'onBits',4,...  % Number of bits per gene which will be "ON". try 6 with 16 total bits for 448 genes  
%    'bitsPerProbe',2,... % bits per targeting sequence.  Increase to 3 or 4 if you have only a few oligos per gene.
%    'universal',false,... % toggle true to use the first primer sequence as a univeral primer for all sublibraries   
%    'primerFasta','', Path to primer sequences, for example: [dataFolder,'PrimerSeqs.fasta'],...
%    'readoutFasta','', Path to readout sequences, for example: [dataFolder,'ReadoutSeqs.fasta'],...
%    'cntrlFasta','',Path to random sequences for negative controls, for example: [dataFolder,'BLASTedRandomSeqs.fasta'],...
%    'libSaveFolder','',Save path [saveFolder,'AssembledProbes\']);
%    'dataFolder','',Path containing fasta files: BLASTedRandomSeqs.fasta, ReadoutSeqs.fasta, and PrimerSeqs.fasta.    
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
% Sequences to Assemble
defaults(end+1,:) = {'primerFasta', 'string', ''};
defaults(end+1,:) = {'readoutFasta', 'string', ''};
defaults(end+1,:) = {'cntrlFasta', 'string', ''};
defaults(end+1,:) = {'dataFolder', 'string', ''};
defaults(end+1,:) = {'libSaveFolder', 'string', ''};
% Library properties
defaults(end+1,:) = {'numCntrls', 'integer', 0};  % random sequences for negative controls 
defaults(end+1,:) = {'numBlanks', 'integer', 0}; % codewords left blank as controls
defaults(end+1,:) = {'numOligos', 'integer', 200}; % total oligos per probe
defaults(end+1,:) = {'subLibNum', 'integer', 1}; % starting sublibrary (increase to avoid using the same index primers as for a different experiment); 
defaults(end+1,:) = {'universal', 'boolean', true};

% Codebook Properties
defaults(end+1,:) = {'numTotalBits', 'integer', 16};  % Total number of bits to read out (must have at least this number of secondary sequences 
defaults(end+1,:) = {'numDataBits', 'integer', 11}; % 
defaults(end+1,:) = {'onBits', 'integer', 4}; % Number of bits per gene which will be "ON"
defaults(end+1,:) = {'bitsPerProbe', 'integer', 2}; % bits per targeting sequence

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabSTORM:invalidArguments', 'required: oligo_folder');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

if ~isempty(parameters.dataFolder)
    parameters.primerFasta = [parameters.dataFolder,'PrimerSeqs.fasta'];
    parameters.readoutFasta = [parameters.dataFolder,'ReadoutSeqs.fasta'];
    parameters.cntrlFasta = [parameters.dataFolder,'BLASTedRandomSeqs.fasta'];   
elseif isempty(parameters.readoutFasta)
    error('either "dataFolder" or "readoutFasta" is required');
else
    parameters.dataFolder = [fileparts(parameters.readoutFasta),filesep];
end

if isempty(parameters.libSaveFolder)
    parameters.libSaveFolder = [dataFolder,'AssembledProbes\'];  %  save folder
end

if ~exist(parameters.libSaveFolder,'dir')
    mkdir(parameters.libSaveFolder);
end

% some book properties
numGenes = length(ProbeData.GeneName);
numWords = numGenes + parameters.numCntrls + parameters.numBlanks; % total number of genes and controls in library


%% Main Function


% Load primers and chose forward and reverse primers randomly
indexPrimers = fastaread(parameters.primerFasta);

% Chose a universal and remove from list
if parameters.universal
    universal = indexPrimers(1).Sequence;
else
    universal = '';
end

% Load random non-homology sequences
cntrl = fastaread(parameters.cntrlFasta);

% Load secondary sequences
readoutSeqs = fastaread(parameters.readoutFasta); %These are the secondary sequences

%% Experiment 1

% A little rounding so the tile size distribution works out even

parameters.numOligos = floor(parameters.numOligos/nchoosek(parameters.onBits,parameters.bitsPerProbe))*nchoosek(parameters.onBits,parameters.bitsPerProbe);


% Generate SECDED codebook with indicated # of ON bits
codebookTable = GenSECDED(parameters.numTotalBits,parameters.numDataBits,parameters.onBits);
[numPosGenes,numHybes] = size(codebookTable);
if numPosGenes < numWords
    error('codebook is not large enough to encode all the indicated genes and controls');
elseif numPosGenes > numWords  
    parameters.numBlanks = numPosGenes - parameters.numCntrls - numGenes;
end
disp(['constructed codebook for up to ',num2str(numPosGenes),' genes in ',num2str(numHybes),' total bits']);



% Add control and Blank sequences
pickedGenes = ProbeData;
pickedGenes = AddCntrlSeqs(pickedGenes,cntrl,...
    parameters.numOligos,parameters.numCntrls,parameters.numBlanks);

[oligoSeq, oligoName] = GenProbe(codebookTable,pickedGenes,readoutSeqs,indexPrimers,...
    parameters.numOligos,numGenes+parameters.numCntrls,parameters.libSaveFolder,...
    'bitsPerProbe',parameters.bitsPerProbe,...
    'subLibNum',parameters.subLibNum,...
    'universal',universal);

%
disp('Assembly complete');
disp(['Library data written to ',parameters.libSaveFolder]);

%%

     
