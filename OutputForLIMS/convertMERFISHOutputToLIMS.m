% Native MERFISH probe design output into LIMS-compatible CSVs

%-----------------------------------------------------------%
% Inputs
% Change for each file
oligosFASTA = 'C:\Users\ScanningLabAnalysis\Documents\MATLAB\MERFISH_analysis\MERFISH_Examples2\Musmusculus\libraryDesign\Musmusculusv3\L1E5_oligos.fasta';
fpkmFile = 'C:\Users\ScanningLabAnalysis\Documents\MATLAB\MERFISH_analysis\MERFISH_Examples2\Musmusculus\Mus_musculus_proxy.fpkm_tracking';
codebookCSV = 'C:\Users\ScanningLabAnalysis\Documents\MATLAB\MERFISH_analysis\MERFISH_Examples2\Musmusculus\codebookMusmusculusHypothalamus_v01.csv';

% oligosFASTA = 'C:\Users\ScanningLabAnalysis\Documents\MATLAB\MERFISH_analysis\MERFISH_Examples2\libraryDesign\Human_smELT_v3_oligos.fasta';
% codebookCSV = 'C:\Users\ScanningLabAnalysis\Documents\MATLAB\MERFISH_analysis\MERFISH_Examples2\codebookHumansmELTv2.csv';

outputFolder = fileparts(codebookCSV);

%-----------------------------------------------------------%
% Config data
% Change only if necessary

nameConversionFile = 'C:\Users\ScanningLabAnalysis\Documents\MATLAB\MERFISH_analysis\OutputForLIMS\mFISH Probes_Allen Institute - Readout + Ear Sequences.csv';

spacerNucleotide = 'A';

%%
%-----------------------------------------------------------%
% Processing code

%-----------------------------------------------------------%
% Name conversion file

fID = fopen(nameConversionFile, 'r');
nameConversion = textscan(fID, '%s%s%s%s%q%s', 'delimiter', ',', 'emptyvalue',0);
fclose(fID);

%%
%-----------------------------------------------------------%
% Codebook processing
fID = fopen(codebookCSV, 'r');
codebook = textscan(fID, '%s', 'Delimiter', '\n');
fclose(fID);

% Parse codebook file

% Library name at top of codebook file
libraryName = codebook{1}{~cell2mat(cellfun(@isempty, strfind(codebook{1}, 'codebook_name'), 'UniformOutput', 0))};
libraryName = libraryName(strfind(libraryName, ', ')+2:end);

% Readouts on 'bit_names' line
readoutsUsed = codebook{1}{~cell2mat(cellfun(@isempty, strfind(codebook{1}, 'bit_names'), 'UniformOutput', 0))};
readoutsUsed = strtrim(strsplit(readoutsUsed, ',')); % Split comma-separated list into cell array
readoutsUsed(1) = []; % First one is header ID string, so delete

% Convert all readouts into internal and external names
% Internal names always start with AIBS0...
readoutsExternal = cell(length(readoutsUsed), 1);
for k = 1:length(readoutsUsed)
    externalCheck = ~cellfun(@isempty, strfind(nameConversion{1}, readoutsUsed{k}));
    internalCheck = ~cellfun(@isempty, strfind(nameConversion{2}, readoutsUsed{k}));
    
    if any(externalCheck)
        readoutsExternal{k} = readoutsUsed{k}; % this readout name is the external one
        readoutsUsed{k} = nameConversion{2}{externalCheck}; % convert used name to the internal name
    else
        readoutsExternal{k} = nameConversion{1}{internalCheck}; % This readout is an internal name, so fill in the external name (if exists)
    end
end

% Genes and code after 'name, id, barcode' line
entries = codebook{1}((find(~cellfun(@isempty, strfind(codebook{1}, 'name, id, barcode')))+1):length(codebook{1}));
genes = cell(length(entries), 1);
code = cell(length(entries), 1);

lastBitUsed = 1;

for k = 1:length(entries)
   
    splitString = strsplit(entries{k}, ',');
    genes{k} = strtrim(splitString{1});
    code{k} = strtrim(splitString{3});
    
    lastBitUsed = max([lastBitUsed, max(strfind(code{k}, '1'))]);
    
end

%%
%-----------------------------------------------------------%
% Oligomers file to LIMS format
oligos = fastaread(oligosFASTA);

% Read fpkm_tracking file into array
fpkmData = readtable(fpkmFile, 'Delimiter', '\t', 'HeaderLines', 1, 'FileType', 'text');
geneIDassoc = fpkmData{:,[1,4]};

%%
processedOligos = struct('gene', [], ...
                         'geneID', [], ...
                         'seqContents', [], ...
                         'sequenceID', [], ...
                         'targetCount', [], ...
                         'fullSequence', [], ...
                         'hybSequence', [], ...
                         'readoutsHere', [], ...
                         'readoutSequences', []);

targetCount = 0; 
prevGene = '';
maxReadouts = 1;

for k = 1:length(oligos)
    
    contents = strsplit(oligos(k).Header, ' ');
    sequence = strsplit(oligos(k).Sequence, ' ');
    
    contents(1) = []; % First entry in header is always library name
                      % Current issue in code has this not match that in
                      % input codebook in frustrating way. Best to ignore
                      % what is here and stick with codebook line.
    
    % Give some ID values to designate sequence components
    %   Readouts = 1
    %   Target sequence = 2
    %   Primers = 3
    %   Spacers = 4
    sequenceID = zeros(length(sequence), 1);
    
    % Spacers are those in sequence that match spacerNucleotide string
    areSpacers = strcmp(sequence, spacerNucleotide);
    sequenceID(areSpacers) = 4;
    
    
    % Find and cut out primers and spacers from sequence
    % Cut out these fields in header as well
    spacersInSeq = find(areSpacers);
    
    if length(spacersInSeq) ~= 2
        error('Number of spacers in sequence found ~= 2');
    end
    
    sequenceID((spacersInSeq(2)+1):end) = 3;
    sequenceID(1:(spacersInSeq(1)-1)) = 3;
    
    seqContents = cell(length(sequenceID), 1);
    contentsCounter = 1;
    for s = 1:length(seqContents)
       
        if sequenceID(s) == 4
            seqContents{s} = 'spacer';
        else
            seqContents{s} = contents{contentsCounter};
            contentsCounter = contentsCounter + 1;
        end
        
    end
                    
    % Readout sequences are those that match to 'Readouts' string in
    % codebookCSV header
    areReadouts = (ismember(seqContents, readoutsUsed) | ismember(seqContents, readoutsExternal));
    sequenceID(areReadouts) = 1;
    
    % Confirm that each potential gene in sequence is actually on the gene 
    % list from the codebook
    
    potentialGenes = find(sequenceID == 0);
        
    for pG = potentialGenes'
       
        geneString = seqContents{pG};
        transcriptString = geneString((strfind(geneString, '__')+2):end);
        transcriptString = transcriptString(1:(strfind(transcriptString, '__')-1));
        geneString = geneString(1:(strfind(geneString, '__')-1));
        
        if ~ismember(geneString, genes)
            error('Potential gene %s not found in gene list', geneString);
        else
            seqContents{pG} = geneString;
            sequenceID(pG) = 2;
        end  
        
    end
    
    if sum(sequenceID == 2) > 1
        error('More than 1 gene found in oligomer');
    end
    
    processedOligos(k).sequenceID = sequenceID;
    processedOligos(k).seqContents = seqContents;
    
    % Increment counter for which oligomer this is on the target
    % transcript.  Assuming that these are in FASTA file in order.
    if strcmp(prevGene, geneString)
        targetCount = targetCount + 1;
        processedOligos(k).geneID = processedOligos(k-1).geneID;
    else
        prevGene = geneString;
        targetCount = 0;
        processedOligos(k).geneID = geneIDassoc{(strcmp(geneIDassoc(:,1), transcriptString)), 2};
    end
    
    processedOligos(k).gene = geneString;
    processedOligos(k).targetCount = targetCount;
    
    processedOligos(k).fullSequence = strrep(oligos(k).Sequence, ' ', ''); % As ordered; target as in gene
    processedOligos(k).hybSequence = seqrcomplement(sequence{sequenceID == 2}); % In probe, not gene
    
    % Pull out readout names
    % Convert to internal names if needed
    readoutNamesHere = seqContents(sequenceID  == 1);
    readoutsInternalHere = cell(length(readoutNamesHere), 1);
    for rOh = 1:length(readoutNamesHere)
        if (ismember(readoutNamesHere{rOh}, readoutsExternal))
            readoutsInternalHere{rOh} = readoutsUsed{strcmp(readoutNamesHere{rOh}, readoutsExternal)};
        else
            readoutsInternalHere{rOh} = readoutNamesHere{rOh};
        end
    end
    processedOligos(k).readoutsHere = readoutsInternalHere;
    
    processedOligos(k).readoutSequences = cell(sum(sequenceID == 1), 1); % As readout probe
    readoutCounter = 1;
    for rOsQ = find(sequenceID == 1)'
        processedOligos(k).readoutSequences{readoutCounter} = seqrcomplement(sequence{rOsQ});
        readoutCounter = readoutCounter + 1;
    end
    
    maxReadouts = max([maxReadouts, readoutCounter-1]);
    
end

%% 
% Write parsed codebook to LIMS format
%  Header lines : gene panel name, readout strings in order
%  Contents : gene names, binary code

fileName = fullfile(outputFolder, strcat(libraryName, '_LIMScodebook.csv'));
fID = fopen(fileName, 'w+');

fprintf(fID, '#CodebookFile\n');
fprintf(fID, '#version\t:\t1.0\n');
fprintf(fID, '#GenePanelOligoPoolName\t:\t%s\n', libraryName);
readoutFormatString = repmat('%s,', 1, lastBitUsed-1);
readoutFormatString = strcat('#ReadoutsUsed\t:\t', readoutFormatString, '%s\n');
fprintf(fID, readoutFormatString, readoutsUsed{1:lastBitUsed});

roundFormatString = repmat('%d,', 1, lastBitUsed-1);
roundFormatString = strcat('#Round\t:\t', roundFormatString, '%d\n');
fprintf(fID, roundFormatString, floor(0:0.5:((lastBitUsed-1)*0.5)));

channelString = strrep(num2str(2*(mod(((0:0.5:((lastBitUsed-1)*0.5))), 1))), '  ', ',');
channelString = strcat('#Channel\t:\t', channelString, '\n');
fprintf(fID, channelString);

wavelengthString = strrep(channelString, 'Channel', 'ExWavelength');
wavelengthString = strrep(wavelengthString, '0', '640');
wavelengthString = strrep(wavelengthString, '1', '750');
fprintf(fID, wavelengthString);

fprintf(fID, 'GeneTarget,Code\n');

for k = 1:length(genes)
    
    fprintf(fID, '%s,%s\n', genes{k}, strtrim(code{k}(1:lastBitUsed)));
end

fclose(fID);


%% 
% Write parsed oligmers file to LIMS format
%  Header lines : gene panel name
%  Contents : gene name, targetCount, full sequence, hyb sequence, [readout
%  name EXT, readout name AIBS, readout sequence] * N

fileName = fullfile(outputFolder, strcat(libraryName, '_LIMSGenePanelOligoPool.csv'));
fID = fopen(fileName, 'w+');

fprintf(fID, '#GenePanelOligoPoolFile\n');
fprintf(fID, '#GenePanelOligoPool.name\t:\t%s\n', libraryName);
fprintf(fID, '#ReferenceGenome.name\t:\t%s\n', 'REFERENCEGENOME');
fprintf(fID, '#Species.name\t:\t%s\n', 'Mus musculus');

readoutHeaderString = '';
for k = 1:maxReadouts
    readoutHeaderString = strcat(readoutHeaderString, ...
                                sprintf(',ReadoutName_%d,ReadoutExternalName_%d,ReadoutSequence_%d', k, k, k));
end

fprintf(fID, 'GeneTarget,TargetCount,GeneID,FullSequence,HybSequence%s\n', ...
        readoutHeaderString);



for k = 1:length(processedOligos)

    % Internal name from GenePanelOligoPool generation output file
    readoutsHere = processedOligos(k).readoutsHere;
    
    readoutsExternalHere = cell(length(readoutsHere), 1);
    for r = 1:length(readoutsExternalHere)
        readoutsExternalHere{r} = readoutsExternal{strcmp(readoutsUsed, readoutsHere{r})};
    end
    
    readoutCell = cell(3, length(readoutsExternalHere));
    for r = 1:length(readoutsExternalHere)
       
        if r <= length(readoutsHere)
            readoutCell{1, r} = readoutsHere{r};
            readoutCell{2, r} = readoutsExternalHere{r};
            readoutCell{3, r} = processedOligos(k).readoutSequences{r};
        else
            readoutCell{1, r} = '';
            readoutCell{2, r} = '';
            readoutCell{3, r} = '';
        end
        
    end
    
    lineFormatString = repmat('%s,%s,%s,', 1, length(readoutsExternalHere)-1);
    lineFormatString = strcat('%s,%03d,%s,%s,%s,', lineFormatString, '%s,%s,%s\n');
    
    fprintf(fID, lineFormatString, ...
        processedOligos(k).gene, ...
        processedOligos(k).targetCount, ...
        processedOligos(k).geneID, ...
        processedOligos(k).fullSequence, ...
        processedOligos(k).hybSequence, ...
        readoutCell{:});
end

fclose(fID);





