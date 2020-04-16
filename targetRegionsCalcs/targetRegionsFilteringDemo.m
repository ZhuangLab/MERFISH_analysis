trPath = 'D:\Data\MERFISH\Homosapiens\HumanMTG_TargetRegionsTester\tr_GC_43_63_Tm_66_76_Len_30_30_IsoSpec_0.00_1.00_Spec_0.75_1.00';

minNumberOfProbes = 67;
printTarget = 1;

geneListDefinedBy = 'allGenes'; % 'codebook' or 'allGenes'

%%
% Load TargetRegions object
tR = TargetRegions.Load(trPath);

%% Output starting TargetRegions to file
baseFolder = 'C:\Users\ScanningLabAnalysis\Documents\MATLAB\MERFISHProbeDesign\targetRegionsCalcs';
outputTargetRegionsTable(fullfile(baseFolder, 'preFilterTargetRegions.csv'), tR, 'FPKMonly', 0);
%%

clear geneIsoformList

switch geneListDefinedBy
    
    case 'codebook'
        
        % Pull target isoforms for each gene
        baseFolder = 'C:\Users\ScanningLabAnalysis\Documents\MATLAB\MERFISHProbeDesign\targetRegionsCalcs';
        cdbkName = 'fakeMTGcodebook.csv';
        codebookPath = fullfile(baseFolder, cdbkName);

        [geneIsoformList, codebookHeader] = LoadCodebook(codebookPath, 'verbose', true);


    case 'allGenes'
        % Do a filtering for ALL TargetRegions in the slicedTranscriptome
        allGeneNames = slicedTranscriptome.GetNames();
        geneIsoformList(length(allGeneNames)) = struct('name', '', ...
                                                       'id', '');
        for k = 1:length(allGeneNames)
            idsHere = slicedTranscriptome.GetIDsByName(allGeneNames{k});
            abundsHere = slicedTranscriptome.GetAbundanceByID(idsHere{1});
            topAbundanceIsoform = idsHere{1}{abundsHere == max(abundsHere)};

            geneIsoformList(k).name = allGeneNames{k};
            % Include only most abundant isoform
            geneIsoformList(k).id = topAbundanceIsoform;
        end
        
    otherwise
        error('geneListDefinedBy must be "codebook" or "allGenes"\n');
end



%% Return regions within given parameter range
filterField = 'isoSpecificity';
parameterRange = [0.75, 1.0];

filteredTargetRegions = findFilteredTargetRegions(tR, geneIsoformList, minNumberOfProbes, printTarget, filterField, parameterRange);
outputTargetRegionsTable(fullfile(baseFolder, 'isospecThresholdedTargetRegions.csv'), ...
    filteredTargetRegions, sprintf('isospecificityThreshold=[%.2f, %.2f]', min(parameterRange), max(parameterRange)), minNumberOfProbes);

%% Return common regions 
commonTargetRegions = findCommonTargetRegions(tR, geneIsoformList, minNumberOfProbes, printTarget);
outputTargetRegionsTable(fullfile(baseFolder, 'commonTargetRegions.csv'), ...
    commonTargetRegions, 'commonRegions', minNumberOfProbes);

%% Return regions by relaxing isospecificity
expandedTargetRegions = findExpandedIsospecificityTargetRegions(tR, geneIsoformList, minNumberOfProbes, printTarget);
outputTargetRegionsTable(fullfile(baseFolder, 'relaxIsospecTargetRegions.csv'), ...
    expandedTargetRegions, 'relaxIsospec', minNumberOfProbes);








