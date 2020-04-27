% ------------------------------------------------------------------------
% filteredTargetRegions = findCommonTargetRegions(targetRegions, ...
%                                                 geneIsoformList, ...
%                                                 minNumberOfProbes, ...
%                                                 printTarget, ...
%                                                 filterField, ...
%                                                 parameterRange)
% 
% Filtering function for targetRegions object. Returned target regions in 
% filteredTargetRegions have following criteria enforced:
% - 1 isoform per gene (which isoform defined by geneIsoformList)
% - This isoform appears in targetRegions
% - The selected isoform has at least minNumberOfProbes target regions
% - Regions returned have value of filterField in range parameterRange
%   
% Filtering is agnostic to other isoforms of the same gene.
% 
% Inputs : targetRegions - TargetRegions object to filter
%          geneIsoformList - struct; fields .name (gene name) and .id (isoform ID)
%          minNumberOfProbes - numeral; threshold for minimum number of probes on an isoform
%          printTarget - '1' or file identifier; where to print output messages. '1' to terminal.
%          filterField - string; match to property in TargetRegions object
%          parameterRange - 1 x 2 vector; pass range for values in filterField
% 
% Returns : filteredTargetRegions - TargetRegions object filtered as described above
%          
% PRN AIBS 2020

function filteredTargetRegions = findFilteredTargetRegions(tR, geneIsoformList, ...
                                                           minNumberOfProbes, printTarget, ...
                                                           filterField, parameterRange)

regsFound = cell(length(geneIsoformList), 4);

% Loop over target genes
    for t = 1:length(geneIsoformList)
        goodProbes = false;
        maxFoundProbes = 0;

        fprintf(printTarget, '-----------------------------------------\n');
        fprintf(printTarget, 'Searching for targets in gene %s\n', geneIsoformList(t).name);

        rowsHere = strcmp({tR.geneName}, geneIsoformList(t).name) & ...
            strcmp({tR.id}, geneIsoformList(t).id);

        if ~any(rowsHere)
            % Isoform not found in targetRegions object
            fprintf(printTarget, 'Gene %s and isoform %s not found in targetRegions object\n', ...
                geneIsoformList(t).name, geneIsoformList(t).id);
            regsFound(t,:) = {geneIsoformList(t).name, maxFoundProbes, -1, []};
            continue;
        end
        
        
        probesToRetain = tR(rowsHere).(filterField) > min(parameterRange) & ...
                         tR(rowsHere).(filterField) <= max(parameterRange);
        
       if sum(probesToRetain) <= minNumberOfProbes
            % No isospecificty value returned enough probes
            fprintf(printTarget, 'Insufficient probes in %s range [%.2f, %.2f] for gene %s and isoform %s\n',...
                filterField, parameterRange, geneIsoformList(t).name, geneIsoformList(t).id);
            regsFound(t,:) = {geneIsoformList(t).name, sum(probesToRetain), -2, []};
            continue;
       end
        
        % We found an isoform with enough probes!
        fprintf(printTarget, 'Gene %s gives %.d probes with %s range [%.2f, %.2f]\n', ...
            geneIsoformList(t).name, sum(probesToRetain), filterField, parameterRange);

        regsFound(t,:) = {geneIsoformList(t).name, sum(probesToRetain), 1, probesToRetain};
       
        
        
    end
    
    % Display how many genes have insufficient probes
    statusTally = vertcat(regsFound{:,3}) < 0;
    fprintf(printTarget, '%d of %d targets fail common regions selection\n', ...
        sum(statusTally), numel(statusTally));

    % Filter out targetRegions given results here
    % Want only genes + isoforms that appear in geneIsoformList
    % Want only those targetRegions that are common across all isoforms

    allIsoformIds = {tR.id}';
    passedIsoformIds = {geneIsoformList(~statusTally).id};

    toKeepTRs = tR(ismember(allIsoformIds, passedIsoformIds));
    
    if ~isempty(toKeepTRs)
        filteredTargetRegions(length(toKeepTRs)) = TargetRegions();

        for k = 1:length(toKeepTRs)
            % Keep only the targetRegions selected to retain in above algorithm
            goodTargetRegions = regsFound{strcmp(regsFound(:,1), toKeepTRs(k).geneName), 4};

            filteredTargetRegions(k) = TargetRegions('geneName', toKeepTRs(k).geneName, ...
                'id', toKeepTRs(k).id, ...
                'geneSequence', toKeepTRs(k).sequence(goodTargetRegions), ...
                'startPos', toKeepTRs(k).startPos(goodTargetRegions), ...
                'regionLength', toKeepTRs(k).regionLength(goodTargetRegions), ...
                'GC', toKeepTRs(k).GC(goodTargetRegions), ...
                'Tm', toKeepTRs(k).Tm(goodTargetRegions), ...
                'specificity', toKeepTRs(k).specificity(goodTargetRegions), ...
                'isoSpecificity', toKeepTRs(k).isoSpecificity(goodTargetRegions));

        end
    else
        % Output empty TargetRegions object
        filteredTargetRegions = TargetRegions();
    end
end