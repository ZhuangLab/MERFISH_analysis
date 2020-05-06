% Harvard way of returning sufficient probes when isospecificity is <
% threshold
% Instead set threshold to maximum value which still returns sufficient probes
% 
% They implement interatively.  Of course this is going to be slow.
% 
% Instead we can look at targetRegions output and tabulate N targetRegions with 
% isospecificity <= unique value.  Then take max that yields num probes above 
% desired cut-off.
% 
% This biases towards probes with higher isospecificity.  The highest, single-isoform
% probes will be included first, then those that overlap with lowest-expressed
% other isoform, and so on until sufficient number is found.  If this is never found,
% then that target is not good for MERFISH experiment. 

baseFolder = 'C:\Users\Rusty Nicovich\Documents\MATLAB\targetRegionsCalcs';

targetRegionsFile = 'target_regions_CUX2_LHX6_GAD1_iso_0_1.csv';
fpkmFile = 'MTG.cleaned.isoforms.fpkm_tracking';

%%
% Load files

% Table with colums
% [geneName, isoformID, isospecificity, specificity, sequence]
fprintf(1, 'Loading target regions file %s\n', fullfile(baseFolder, targetRegionsFile));
tR = readtable(fullfile(baseFolder, targetRegionsFile), 'headerLines', 1);

% Table with colums
% [isoformID, -, -, geneID, geneName, -, -, -, -, FPKM, -, -, -]
% '-' either dummy data or not needed here
fprintf(1, 'Loading FPKM file %s\n', fullfile(baseFolder, targetRegionsFile));
fpkm = readtable(fullfile(baseFolder, fpkmFile), 'headerLines', 1, 'filetype', 'text');


% Pull target isoforms for each gene
% Stand-in for codebook information on targets
% Highest-expressed isoform chosen for each (remember proxy FPKM files are
% BS)
targets = struct('geneName', '', ...
                 'isoformID', '');

targets(1).geneName = 'GAD1'; 
targets(1).isoformID = 'NM_000817.2';

targets(2).geneName = 'LHX6'; 
targets(2).isoformID = 'NM_001242334.1';

targets(3).geneName = 'CUX2'; 
targets(3).isoformID = 'XM_011538071.1';

% Minimum number of probes that must be present to be above 'background'
minNumberOfProbes = 67; 

%%
% Loop over target genes
for t = 1:length(targets)
    goodProbes = false;
    fprintf(1, '-----------------------------------------\n');
    fprintf(1, 'Searching for targets in gene %s\n', targets(t).geneName);

    rowsHere = ~cellfun(@isempty, strfind(tR{:,1}, targets(t).geneName)) & ...
               ~cellfun(@isempty, strfind(tR{:,2}, targets(t).isoformID));

    % For target isoform, descending sort of available isospecificity
    % scores
    isospecsHere = sort(unique(tR{rowsHere, 3}), 'ascend');
 
    % Count occurences of each of these isospecificity scores
    isospecCounts = histc(tR{rowsHere, 3}, (isospecsHere));
    
    % Find threshold isospecificity to cross minNumberOfProbes value
    isospecsDescend = flipud(isospecsHere);
    cumCounts = cumsum(flipud(isospecCounts));
    whichValHasEnoughCounts = find(cumCounts > minNumberOfProbes, 1, 'first');
    maxIsospecWithEnoughProbes = isospecsDescend(whichValHasEnoughCounts);    
    
    if isempty(maxIsospecWithEnoughProbes)
        % No isospecificty value returned enough probes
        fprintf(1, 'No isospecificity value found to yield enough probes for gene %s and isoform %s\n',...
            targets(t).geneName, targets(t).isoformID);
        
    else
        % We found a good isospecificity value!
        fprintf(1, 'Gene %s gives %.d probes with isospecificty value %f\n', ...
            targets(t).geneName, cumCounts(whichValHasEnoughCounts), maxIsospecWithEnoughProbes);
    end
end
    
    
    
    
    
    