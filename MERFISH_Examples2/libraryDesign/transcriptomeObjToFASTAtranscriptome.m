% Convert transcriptome object on disk back to a FASTA file and abundance
% file

% Path supplied by Moffitt:
transcriptomePath = 'C:\Users\ScanningLabAnalysis\Documents\MATLAB\MERFISH_analysis\MERFISH_Examples2\Musmusculus\150603_Transcriptome\transcriptomeC57Obj';

% OutPath and file name for FASTA file
outPath = 'C:\Users\ScanningLabAnalysis\Documents\MATLAB\MERFISH_analysis\MERFISH_Examples2\Musmusculus';
outName = 'Mus_musculus_regenTranscriptome.fasta';

% Out name for isoforms abundance file
outputFile = 'Mus_musculus_proxyByWayOfMoffitt.fpkm_tracking';

%% Load transcriptome on disk
transcriptome = Transcriptome.Load(transcriptomePath);

%% Formulate FASTA object and along the way write to isoforms file
transOut = struct('Header', '', ...
                  'Sequence', '');

nameList = transcriptome.GetNames();
onGene = 1;

% isoforms file
fID = fopen(fullfile(outPath, outputFile), 'w+');
fmt = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n';
fprintf(fID, 'tracking_id\tclass_code\tnearest_ref_id\tgene_id\tgene_short_name\ttss_id\tlocus\tlength\tcoverage\tFPKM\tFPKM_conf_lo\tFPKM_conf_hi\tFPKM_status\n');
              
for k = 1:transcriptome.numGenes
    
    geneHere = nameList{k};
    IDsHere = transcriptome.GetIDsByName(nameList{k});
    
    for i = 1:length(IDsHere{1})
        
        IDnow = IDsHere{1}{i};
        
        
        % Formulate header w/ ID number and gene name
        headerString = sprintf('%s gene=%s', IDnow, geneHere);
        
        seq = transcriptome.GetSequenceByID(IDnow);
        
        transOut(onGene).Sequence = seq;
        transOut(onGene).Header = headerString;
        
        onGene = onGene + 1;
        
        
        % Isoforms file
        % tracking_id - transcript ID, 'ENSMUST00000140889.7'
        % class_code - '-'
        % nearest_ref_id - '-'
        % gene_id - gene ID, 'ENSMUSG00000044700.15'.  Need to BS value here.
        % gene_short_name - gene_symbol in header, 'Tmem201'
        % tss_id - unclear, unimportant? 'TSS0000'
        % locus - chromosome and bp location, 'chr4:149728790-149737998'. (Certainly incorrect)
        % length - length of transcript, get from length of Sequence unit 
        % coverage - ? '0'
        % FPKM - '0'
        % FPKM_conf_lo - '0'
        % FPKM_conf_hi - '0'
        % FPKM_status - 'OK'    
        
        abund = transcriptome.GetAbundanceByID({IDnow});
        
        fprintf(fID, fmt, ...
            IDnow, ...
            '-', ...
            '-', ...
            'ENSMUSG00000000000.0', ...
            geneHere, ...
            'TSS0000', ...
            sprintf('chr%s:%s-%s', '0', '0', num2str(numel(seq))), ...
            num2str(numel(seq)), ...
            '1', ...
            sprintf('%f', abund), ...
            '0', ...
            '0', ...
            'OK');
        
        
        
    end
end

fclose(fID);



%% Write to FASTA file

fastawrite(fullfile(outPath, outName), transOut);