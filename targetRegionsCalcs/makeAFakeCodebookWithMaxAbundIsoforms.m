% Make a fake codebook from provided codebook, FPKM_tracking file
% Want 1 line per gene, with isoform ID = max abundance isoform
baseFolder = 'C:\Users\ScanningLabAnalysis\Documents\MATLAB\MERFISHProbeDesign\targetRegionsCalcs';

codebookPath = fullfile(baseFolder, 'Human_MTG_Panel1_corrected_3_barcoded_truncated.csv');
[codebook, codebookHeader] = LoadCodebook(codebookPath, 'verbose', true);


fprintf(1, 'Loading FPKM file %s\n', ...
    'D:\Data\MERFISH\Homosapiens\Homo_sapiens_proxyRandomFPKM.fpkm_tracking');
% [isoformID, -, -, geneID, geneName, -, -, -, -, FPKM, -, -, -]
fpkm = readtable(fullfile(baseFolder, fpkmFile), 'headerLines', 1, 'filetype', 'text');

%% Isoform selection method
% 'mostAbundant' or 'longest'

isoformPickMethod = 'mostAbundant';
versionMatch = true;


%%
codebookNames = unique({codebook.name});

cdbk = struct('name', '', 'id', '', 'barcode', '');

for t = 1:length(codebookNames)
    cdbk(t).name = codebookNames{t};
        
    cdbk(t).barcode = codebook(strcmp({codebook.name}, codebookNames{t})).barcode;

    switch isoformPickMethod
        case 'mostAbundant'
            % Pick the isoform with the highest value in the fpkm table
            % Best when you have real fpkm data to use
            
            isoforms = fpkm{strcmp(fpkm{:,5}, cdbk(t).name), 1};
            abunds = zeros(length(isoforms), 1);
            for k = 1:length(isoforms)

               abunds(k) =  fpkm{strcmp(fpkm{:,1}, isoforms{k}), 10};

            end
            
            fpkmIsoform = isoforms{find(abunds == max(abunds), 1, 'first')};
            
            if versionMatch
                cdbk(t).id = fpkmIsoform;
                
            else
                % Find the right isoform, which may not version match
                
                tRTrunc = tR(strcmp({tR.geneName}, cdbk(t).name));
                findCloseIsoformID = ~cellfun(@isempty, strfind({tRTrunc.id}, fpkmIsoform(1:strfind(fpkmIsoform, '.'))));
                if any(findCloseIsoformID)
                    findAnIsoform = tRTrunc(findCloseIsoformID).id;
                    cdbk(t).id = findAnIsoform;
                else
                    cdbk(t).name = '';
                end
            end
            
        case 'longest'
            % Choose longest isoform as target
            % Best when you have no fkpm data or useUniformWeights = TRUE
            tRTrunc = tR(strcmp({tR.geneName}, codebook(t).name));
            isoforms = {tRTrunc.id};
            
            % numRegions is proxy for length
            numRegs = vertcat(tRTrunc.numRegions);
            
            if ~isempty(tRTrunc)
                cdbk(t).id = isoforms{find(numRegs == max(numRegs), 1, 'first')};
            else
                cdbk(t).name = '';
            end
            
            
    end
    
end

%%
% Output to file
fID = fopen(fullfile(baseFolder, 'fakeMTGcodebook.csv'), 'w+');

fprintf(fID, 'version, %s\n', codebookHeader.version);
fprintf(fID, 'codebook_name, fakeMTGcodebook.csv - %s\n', isoformPickMethod);
fprintf(fID, ['bit_names, ', ...
    repmat('%s, ', 1, length(codebookHeader.bit_names) - 1),...
    '%s\n'], ... 
    codebookHeader.bit_names{:});
fprintf(fID, 'name, id, barcode\n');
for k = 1:length(cdbk)
    if ~(strcmp(cdbk(k).name, ''))
        fprintf(fID, '%s, %s, %s\n', cdbk(k).name, cdbk(k).id, cdbk(k).barcode);    
    end
end
fclose(fID);