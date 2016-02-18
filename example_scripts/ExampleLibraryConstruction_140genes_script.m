%% Example Designing Probes to Target 140 different genes
% -------------------------------------------------------------------------
% Alistair Boettiger
% February 18, 2016
% boettiger.alistair@gmail.com
% -------------------------------------------------------------------------
% Purpose: To illustrate the construction of a MERFISH probe library.
% -------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.
%
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
% Required Files
% -------------------------------------------------------------------------
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
% 5. OligoArray 2.1 installation (see http://berry.engin.umich.edu/oligoarray2_1/) 
% 6. Legacy BLAST installation (also required for OligoArray 2.1).
% 
% Important: the probe
% designer will search for unique regions of each gene in this fasta file
% from which to build probes.  Including multiple isoforms of the same gene
% in this fasta file may result in very few or no unique sequences for that
% gene.  Additional genes included in this fasta which are not
% selected for the library will still be used in the unique region
% screening. 
% 
% -------------------------------------------------------------------------
% Output Files: 
% -------------------------------------------------------------------------
% 1. E1_codebook  - fasta file recording the codeword assigned to each
%                   gene, and the names of oligos assigned to each bit.
% 2. E1_oligos   - fasta file containing all the oligos in the library. the
%                   header records the experiment number, the and lists the
%                   binding regions included in the probe as a
%                   space-delimited entry.  The sequence is also space
%                   delimited.  This fasta file may be concatinated with
%                   other fasta files which employ orthoganol index
%                   primers.
%    Example probe:
%    > E1 primer_002 B02 ABCA2____107 B03 primer_003 probe005
%    CGCGGGCTATATGCGAACCG  GCTATCGTTCGTTCGAGGCCAGAGCATTCG TGGCCCAAGAAGGAGACCACCTGGTTGCGG 
%    CCCATGATCGTCCGATCTGGTCGGATTTGT GCGTTGTATGCCCTCCACGC
%    is the 5th probe in the library, this one binds starting at bp 107 of the gene
%    ABCA2.  This probe will be bound by readout-sequence 2 and readout-
%    sequence 3.  It is part of experiment 1, which can be amplified using 
%    primer sequences 2 and 3 (the 20 mers at each end).  
% -------------------------------------------------------------------------


% Add all MERFISH functions to your filepath (if not already active):
addpath(genpath('C:\Users\Alistair\Documents\Research\Projects\MERFISH-Release\MERFISH-public\'))

% Tell matlab where to find OliogArray
oligoArrayExe = 'C:\Users\Alistair\Documents\Research\Software\OligoArray\OligoArray2.jar'; % UPDATE THIS

% Specify the location of the target sequence data
exampleData = 'C:\Users\Alistair\Documents\Research\Projects\MERFISH-Release\MERFISH-data\MERFISH_Examples\';  % UPDATE THIS 
dataFolder = [exampleData,'probe_construction_data\'];
geneFasta = [dataFolder,'TargetGeneSeqs.fasta'];  % Download our demo sequences or choose your own.  
blastLib =  [dataFolder,'TargetGeneSeqs.fasta']; % Download our demo sequences or choose your own.  
% If there are additional sequences you want to ensure do not capture your probes, you can add these 
% sequences to the fasta file used for the blastLib.  Do not add them to the geneFasta unless you want
% to design probes for them.

% Specify a folder to save the ouptut data in
saveFolder = [exampleData,'probe_construction_output\'];

timeRun = tic;
% Build a BLAST database for OligoArray to use based on the fasta file 
BuildBLASTlib(blastLib,'legacy',true);

[oligoArrayCommand, savePath] = OligoArrayCmd('savePath',saveFolder,...   'saveDir',saveFolder,... 
    'minTm',70,...                      % Minimum Tm of probe for its target
    'probeLength',30,...                % Desired length of targeting region 
    'crosshybeT',72,...                 % Cross-hybridization with a Tm of this value or higher will be considered "off-target" and the probe region ignored
    'secstructT',76,...                 % Minimum Tm of secondary structure allowed in viable probes
    'blastLib',blastLib,...            % BLAST library to use for screening against off-target genes and library cross-talk
    'maskedSeq','"GGGG;CCCC;TTTTTT;AAAA"',... % forbidden sequences (will not be used for probes) 
    'numParallel',5,...                 %  number of cores to use (max number of fragments per gene to use at once);
    'GbMem',1,...                       %  max memory allowed per core 
    'oligoArrayExe',oligoArrayExe);     % complete filepath to OligoArray2.jar

% Read in the fasta file containing the genes for which we want to design probes. 
geneList = fastaread(geneFasta);

BatchLaunchOligoArray(oligoArrayCommand,geneList,...
    'batchsize',20,... % max number of genes to launch at once (If you run into CPU problems, set this to less than the number of cores you have).
    'maxTime',30,...  % in minutes, if an individual gene takes longer than this, probe construction for the gene will be aborted.
    'savePath',savePath,... % place to save the txt files of oligos generated by OligoArray.   
    'maxCPU',95,... % this is the preferred way to throttle CPU, note there is a delay on how frequently CPU checks are run.
    'runExternal',true); % set false to print oligoArray output to the command window. This is useful for troubleshooting. 


%% Parse OligoArray results
% The output structure will be saved as ProbeData.mat inside the savePath.
CompileOligoArrayOutput(savePath,... 
    'blastLib',blastLib,... % blastlib
    'zeroOffTarget',true); % zero-tolerance for off targets - any sequence with an off target hit will be rejected.  
                           % This runs much faster if you set this to false.  

%% Assemble target sequences and readout sequences into library
load([saveFolder,'ProbeData.mat'],'ProbeData'); 

AssembleProbes(ProbeData,...
    'numCntrls',5,... % random sequences for negative controls 
    'numBlanks',5,...  % min number of codewords left blank as controls
    'numOligos',200,... % total oligos per probe
    'subLibNum',1,... % starting sublibrary (increase to avoid using the same index primers as for a different experiment); 
    'numTotalBits',16,...   % Total number of bits to read out (must have at least this number of secondary sequences 
    'onBits',4,...  % Number of bits per gene which will be "ON". try 6 with 16 total bits for 448 genes  
    'bitsPerProbe',2,... % bits per targeting sequence.  Increase to 3 or 4 if you have only a few oligos per gene.
    'universal',false,...
    'primerFasta',[dataFolder,'PrimerSeqs.fasta'],...
    'readoutFasta',[dataFolder,'ReadoutSeqs.fasta'],...
    'cntrlFasta',[dataFolder,'BLASTedRandomSeqs.fasta'],...
    'libSaveFolder',[saveFolder,'AssembledProbes\']);

totTime = toc(timeRun)/60;
disp(['Probe library assembly completed in ',num2str(totTime,5),' min.']);
