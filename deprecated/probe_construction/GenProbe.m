function [OligoSeq, OligoName, primerNum, parameters] = GenProbe(libCodebook,pickedGenes,secondaries,indexPrimers,numOligos,numGenesAndCntrls,saveFolder,varargin)
%--------------------------------------------------------------------------
% Generate Probes based on libray codebook for the indicated genes using
% the indicated secondaries and indicated index primers.
%
%--------------------------------------------------------------------------
% Required Inputs 
%
% libCodebook - matrix codebook: a logical, numGenes x numHybes
% pickedGenes - structure containing the genes to make probes for
% indexPrimers - fasta file structure of index primers to use
% secondaries - fasta file structure of secondary sequences to use
% numOligos - number of probes per gene
% numGenesAndCntrls - total number of sequences to build probes for
% saveFolder - a place to export the data as a fasta files
% 
%--------------------------------------------------------------------------
% Optional Inputs 
% subLibNum / integer / 1
%                      - index of sub library.  used to keep track of which
%                      fwd and reverse primers to use.  
% universal / string / '' 
%                      - sequence to
% shuffleCodewords
% shufflePriSec
% shufflePri
% changeSecondaries
%

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'bitsPerProbe', 'integer',2 };
defaults(end+1,:) = {'subLibNum', 'integer',1 };
defaults(end+1,:) = {'primerNum', 'integer',1 };
defaults(end+1,:) = {'universal', 'string','' };
defaults(end+1,:) = {'distanceSort','boolean',true};
defaults(end+1,:) = {'shuffleCodewords', 'freeType',1 };
defaults(end+1,:) = {'shufflePriSec', 'freeType',1 };
defaults(end+1,:) = {'shufflePri', 'freeType',1 };
defaults(end+1,:) = {'changeSecondaries', 'freeType',0};
defaults(end+1,:) = {'saveIsoName','boolean',true};
% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 7
    error('matlabSTORM:invalidArguments', 'required: libCodebook,pickedGenes,secondaries,indexPrimers,numOligos,numGenesAndCntrls,saveFolder');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

saveIsoName = parameters.saveIsoName;
changeSecondaries = parameters.changeSecondaries;
universal = parameters.universal;
subLibNum = parameters.subLibNum;
primerNum = parameters.primerNum;
distanceSort =parameters.distanceSort;
shuffleCodewords =parameters.shuffleCodewords;
shufflePriSec = parameters.shufflePriSec;
shufflePri = parameters.shufflePri;
bitsPerProbe =  parameters.bitsPerProbe;


%% Main Function
%-------------------------------------------------------------------------
[numWords,numLetters] = size(libCodebook);
numShuffles = max([length(shuffleCodewords),length(shufflePriSec),length(shufflePri),length(changeSecondaries)]);

if ~isempty(universal)
    commonPrimerFlag  = ['universal', ' '];
else
    commonPrimerFlag  = '';
end

% save the random draw events 
origSecondaries = secondaries; 
priSecShuffle = cell(numShuffles,numWords);
priDrawShuffle = cell(numShuffles,numWords); 
codewordShuffle = cell(numShuffles); 

OligoSeq=cell(numGenesAndCntrls*numOligos,numShuffles);
OligoName=cell(numGenesAndCntrls*numOligos,numShuffles);

for r = 1:numShuffles;
    
disp(['Processing sublibrary: ' num2str(subLibNum)]);
    % Select Index Primers
    primerNum = primerNum + 1;
    fwdPrimer = indexPrimers( primerNum ).Sequence;
    fwdPrimerName = indexPrimers( primerNum ).Header;
    primerNum = primerNum + 1;
    revPrimer = seqrcomplement( indexPrimers( primerNum ).Sequence );
    revPrimerName = indexPrimers( primerNum ).Header;

    
    % Select Secondaries
    if changeSecondaries(min(length(changeSecondaries),r))
        secondaries = origSecondaries((r-1)*numLetters+1:r*numLetters); % draw from end
    else
        secondaries = origSecondaries(1:numLetters); 
    end
    
    % shuffle codewords and write codebook file
    if shuffleCodewords(min(length(shuffleCodewords),r))
        libCodebook = libCodebook(randperm(length(libCodebook)),:);  % randomize codeword assignments 
        if distanceSort
            % assign a codeword at random to the first gene.
            % codewords for the remaining genes are listed in order of
            % increasing distance from the first gene.  Therefore if genes
            % are passed in order of descending FPKMs, high expressed genes
            % will be closer to other highly expressed genes and low
            % expressed genes closer to other low expressed genes.  
            wordDistances = squareform(pdist(double(libCodebook),'hamming'));
            % figure(1); clf; imagesc(wordDistances); colormap jet; colorbar;
            distSorted = sortrows([wordDistances(1,:)',(1:numWords)']);
            libCodebook = libCodebook(distSorted(:,2),:);
        end
     
        codewordShuffle{r} = libCodebook;
    else
        libCodebook = codewordShuffle{1};
    end
     if saveIsoName
         geneNames = cellfun(@(x,y) [x,'    ',y],pickedGenes.CommonName,pickedGenes.IsoformName,'UniformOutput',false);
     else
         geneNames = pickedGenes.CommonName;
     end
     codebook.Sequence = strcat(geneNames,['    ',CSL2str({secondaries(1:numLetters).Header})]) ; 
     codebook.Header = mat2cell( num2str(libCodebook),ones(1,size(libCodebook,1)),3*numLetters-2);
     WriteFasta([saveFolder, 'E',num2str(subLibNum),'_','codebook.fasta'],codebook.Header,codebook.Sequence,'Append',false); 
     oligoFile = [saveFolder, 'E',num2str(subLibNum),'_','oligos.fasta'];

    if exist(oligoFile,'file') ~= 0
        warning(['overwriting existing file: ',oligoFile]);
        delete(oligoFile);
    end

    % Select primaries, group secondaries, match primaries and secondaries.
    j=1;
    for n = 1:numWords
        
        % don't try to build sequences from blanks, we don't have any primaries to build on.   
        if isempty( strfind(pickedGenes.CommonName{n},'blank') )

            geneOnBits = find(libCodebook(n,:));  % indices of the 1 bits (out of 16 letters)
            secCombos = nchoosek(geneOnBits,bitsPerProbe);              % all possible pairs combinations (since each tile has 1 pair)
            numSecCombos = size(secCombos,1); 

            % shuffle which secondaries are attached to which targeting sequences (of the numOligo used targeting sequences)
            if shufflePriSec(min(length(shufflePriSec), r))
                onBits = length(geneOnBits); 
                numOligosUsed = min(length([pickedGenes.Sequence{n}]),numOligos);
                numOligosUsed = floor(numOligosUsed/nchoosek(onBits,bitsPerProbe))*nchoosek(onBits,bitsPerProbe);
                randInds = randperm(numOligosUsed);
                priIdxPerSecCombo = reshape(randInds, numSecCombos,numOligosUsed/numSecCombos );  
                priSecShuffle{r,n} = priIdxPerSecCombo; % save to reuse in other shuffles
            else
               priIdxPerSecCombo = priSecShuffle{1,n};
            end

            % randomly draw numOligo subset of sequences to use as targeting sequences. 
            if shufflePri(min(length(shufflePri),r))
                randSeq = randperm(length([pickedGenes.Sequence{n}]));  
                priDrawShuffle{r,n} = randSeq; 
            else
                randSeq = priDrawShuffle{1,n} ;
            end

            seq = [pickedGenes.Sequence{n}]; 
            numOligosUsed = min(length(seq),numOligosUsed);
            seq = seq(randSeq(1:numOligosUsed));
            fiveprime = [pickedGenes.FivePrimeEnd{n}];
            % fiveprime = [pickedGenes.MeltingTemp{n}];
            try
            fiveprime = fiveprime(randSeq(1:numOligosUsed));
            catch
                fiveprime = zeros(1,max(priIdxPerSecCombo(:)));
            end

            for k = 1:numSecCombos
                for i = 1:numOligosUsed/numSecCombos
                    oligosIndex = priIdxPerSecCombo(k,:);
                    secIndex = secCombos(k,:);
                    
                    readoutSeq1 = [seqrcomplement(secondaries(secIndex(1)).Sequence), ' '];
                    readoutSeq2 =  [seqrcomplement(secondaries(secIndex(2)).Sequence), ' '];
                    readoutName1 = [secondaries(secIndex(1)).Header ' '];
                    readoutName2 = [secondaries(secIndex(2)).Header ' '];
                    if length(secIndex) > 2
                        readoutSeq3 = [seqrcomplement(secondaries(secIndex(3)).Sequence),' '];
                        readoutName3 = [secondaries(secIndex(3)).Header ' '];
                    else
                        readoutSeq3 = '';
                        readoutName3 = '';
                    end
                    if length(secIndex) > 3
                        readoutSeq4 = [seqrcomplement(secondaries(secIndex(4)).Sequence) ' '];        
                        readoutName4 = [secondaries(secIndex(4)).Header ' '];
                    else
                        readoutSeq4 = '';
                        readoutName4 = ''; 
                    end
                    OligoSeq(j,r) = ...
                        {[fwdPrimer, ' ',...
                          universal, ' ',...
                          readoutSeq3,...
                          readoutSeq1,...
                          seqrcomplement(char(seq(oligosIndex(i)))), ' ',...
                          readoutSeq2,...
                          readoutSeq4,...
                          revPrimer]};
                    
                    OligoName(j,r) = ...
                          {['E',num2str(subLibNum), ' ',...
                            fwdPrimerName,' ',...
                            commonPrimerFlag,...
                            readoutName3,...
                            readoutName1,...
                            char(pickedGenes.CommonName{n}) '__',...
                            char(pickedGenes.IsoformName{n}) '__',...
                            num2str(fiveprime(oligosIndex(i))),' ',... % ,'i',num2str(1E12*rand,12), 
                            readoutName2,...
                            readoutName4,...
                            revPrimerName,...
                            ' probe',num2str(j,'%.3d')]};
                    WriteFasta([saveFolder, 'E',num2str(subLibNum),'_','oligos.fasta'],OligoName{j,r},OligoSeq{j,r},'Append',true,'Warnings',false); 
                    j = j + 1;

                end    
            end
        end 
    end
    subLibNum = subLibNum+1;
end
