% Fake a bulk sequencing file from a cDNA transcriptome file
% Script name is misnomer - all genes should be included in output fpkm file. 
%							At minimum all isoforms of gene appearing in codebook must appear in fpkm_tracking file.

baseFolder = 'C:\Users\ScanningLabAnalysis\Documents\MATLAB\MERFISH_analysis\MERFISH_Examples2\Musmusculus';

inputTranscriptome = 'Mus_musculus.GRCm38.cdna.all.fa';

outputFile = 'Mus_musculus_proxyRandomFPKM.fpkm_tracking';

%% Read in FASTA cDNA file
cDNAAll = fastaread(fullfile(baseFolder, inputTranscriptome));

%% Output to fake fpkm_tracking file

% Each line in header needs to decompose to:
% tracking_id - transcript ID, 'ENSMUST00000140889.7'
% class_code - '-'
% nearest_ref_id - '-'
% gene_id - gene ID, 'ENSMUSG00000044700.15'
% gene_short_name - gene_symbol in header, 'Tmem201'
% tss_id - unclear, unimportant? 'TSS0000'
% locus - chromosome and bp location, 'chr4:149728790-149737998'. (Certainly incorrect)
% length - length of transcript, get from length of Sequence unit 
% coverage - ? '0'
% FPKM - '0'
% FPKM_conf_lo - '0'
% FPKM_conf_hi - '0'
% FPKM_status - 'OK'

fID = fopen(fullfile(baseFolder, outputFile), 'w+');

fmt = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n';

fprintf(fID, 'tracking_id\tclass_code\tnearest_ref_id\tgene_id\tgene_short_name\ttss_id\tlocus\tlength\tcoverage\tFPKM\tFPKM_conf_lo\tFPKM_conf_hi\tFPKM_status\n');
for k = 1:length(cDNAAll)
    
    headerSplit = strsplit(cDNAAll(k).Header);
    
    
    chromSplit = strsplit(headerSplit{3}, ':');
    
    fprintf(fID, fmt, ...
        headerSplit{1}, ...
        '-', ...
        '-', ...
        headerSplit{4}(strfind(headerSplit{4}, ':') + 1:end), ...
        headerSplit{7}(strfind(headerSplit{7}, ':') + 1:end), ...
        'TSS0000', ...
        sprintf('chr%s:%s-%s', chromSplit{3}, chromSplit{4}, chromSplit{5}), ...
        num2str(numel(cDNAAll(k).Sequence)), ...
        '1', ...
        sprintf('%f', 10^(randn(1))), ...
        '0', ...
        '0', ...
        'OK');
    
end
        

