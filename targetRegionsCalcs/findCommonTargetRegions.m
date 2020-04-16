% ------------------------------------------------------------------------
% filteredTargetRegions = findCommonTargetRegions(targetRegions, ...
%                                                 geneIsoformList, ...
%                                                 minNumberOfProbes, ...
%                                                 printTarget)
% 
% Filtering function for targetRegions object. Returned target regions in 
% filteredTargetRegions have following criteria enforced:
% - 1 isoform per gene (which isoform defined by geneIsoformList
% - This isoform appears in targetRegions
% - The selected isoform has at least minNumberOfProbes target regions
% - Regions returned have same isospecificity
% - There are at least minNumberOfProbes at this isospecificity on target isoform
% - If these probes appear on any other isoform of this gene, that isoform carries at least
%   minNumberOfProbes probes
%   
% The goal is to return a set of 'common' sequences in all isoforms of this gene.
% At the very least, no isoform should carry fewer than minNumberOfProbes probes.
% 
% Inputs : targetRegions - TargetRegions object to filter
%          geneIsoformList - struct; fields .name (gene name) and .id (isoform ID)
%          minNumberOfProbes - numeral; threshold for minimum number of probes on an isoform
%          printTarget - '1' or file identifier; where to print output messages. '1' to terminal.
% 
% Returns : filteredTargetRegions - TargetRegions object filtered as described above
%          
% PRN AIBS 2020


function filteredTargetRegions = findCommonTargetRegions(tR, geneIsoformList, minNumberOfProbes, printTarget)

% Loop over genes in geneIsoformList

regsFound = cell(length(geneIsoformList), 4);

for t = 1:length(geneIsoformList)
    
    % Initialize variables for this gene
    goodProbes = false;
    maxFoundProbes = 0;
    %regionsAlignmentStatus = 0;
    
    
    fprintf(printTarget, '-----------------------------------------\n');
    fprintf(printTarget, 'Searching for common targets in gene %s\n', geneIsoformList(t).name);
    tRTrunc = tR(strcmp(geneIsoformList(t).name, {tR.geneName}));
    whichEntry = strcmp(geneIsoformList(t).id, {tRTrunc.id}) & strcmp(geneIsoformList(t).name, {tRTrunc.geneName});
    
    if ~any(whichEntry)
        % Gene + isoform does not appear in targetRegions object
        % Likely is filtered out at FPKM threshold step
        fprintf(printTarget, 'No target regions found for gene %s and isoform %s\n', ...
            geneIsoformList(t).name, geneIsoformList(t).id);
        regsFound(t,:) = {geneIsoformList(t).name, maxFoundProbes, -1, []};
        continue;
    end
    
    
    % Check that a solution is even possible, given the numRegions
    % available for this gene
    numRegs = horzcat(tRTrunc.numRegions);
    maxLength = max(numRegs);
    fprintf(printTarget, 'Longest isoform of %d has %d target regions\n', length(numRegs), maxLength);
    
    if maxLength <= minNumberOfProbes
        % All isoforms of gene are too short.  There can never be sufficient probes for
        % this isoform
        regsFound(t,:) = {geneIsoformList(t).name, maxFoundProbes, -2, []};
        continue;
        
    end
    
    
    % Can attempt to find sufficient common probes for this isoform.
    
    fprintf(printTarget, 'Target isoform has %d target regions\n', tRTrunc(whichEntry).numRegions);
    
    if tRTrunc(whichEntry).numRegions <= minNumberOfProbes
        % Selected isoform has too few targetRegions.
        % No solution exists for this isoform, but might for another
        % isoform.
        fprintf(printTarget, 'Target isoform %s has too few target regions.\n', geneIsoformList(t).id);
        regsFound(t,:) = {geneIsoformList(t).name, maxFoundProbes, -3, []};
        continue;
    end
    
    
    % For target isoform, ascending sort of isospecificity score
    isospecsHere = unique(tRTrunc(whichEntry).isoSpecificity);
    
    onIsospec = 1;
    
    % Loop over isospec
    % Pick lowest available isospecificity score
    % If none available, then selection fails
    %     for onIsospec = 1:length(isospecsHere)
    while (onIsospec < (length(isospecsHere) + 1)) && ~goodProbes
        
        isospecRows = tRTrunc(whichEntry).isoSpecificity == isospecsHere(onIsospec);
        % Need count of number of isoforms for this gene
        
        isoformsOfThisGene = {tRTrunc(strcmp(geneIsoformList(t).name, {tRTrunc.geneName})).id};
        % nb - if there are only 1 or 2 isoforms then the above test catches issues
        
        % Are there enough probes on target isoform with this score?
        if sum(isospecRows) <= minNumberOfProbes
            % If no, move to next-lowest isospecificity score
            
            fprintf(printTarget, 'Only %.d probes in isoform %s at isospecificity %f\n', ...
                sum(isospecRows), geneIsoformList(t).id, isospecsHere(onIsospec));
            onIsospec = onIsospec + 1;

            regsFound(t,:) = {geneIsoformList(t).name, maxFoundProbes, -4, []};
            maxFoundProbes = max([maxFoundProbes, sum(isospecRows)]);
            continue;
        end
        
        
        % If yes, proceed with cross-isoform checks
        
        
        % Find sequences with this score
        probeSequences = tRTrunc(whichEntry).sequence(isospecRows);
        
        % For each, see how many show up in each of the other isoforms
        nIsoformsHere = zeros(sum(isospecRows), length(isoformsOfThisGene));
        
        for k = 1:length(probeSequences)
            for i = 1:length(isoformsOfThisGene)
                
                % Find entry for gene name, isoform i
                % Isoform ID is unique, so no need to search for geneName
                
                %                   probeFindEntry = tRTrunc(strcmp(isoformsOfThisGene{i}, {tRTrunc.id})).sequence;
                %                   % True if probe sequence k in isoform i
                %                   nIsoformsHere(k, i) = any(~cellfun(@isempty, strfind(probeFindCat, probeSequences{k})));
                %
                probeFindCat = strcat(tRTrunc(strcmp(isoformsOfThisGene{i}, {tRTrunc.id})).sequence{:});
                nIsoformsHere(k, i) = ~isempty(strfind(probeFindCat, probeSequences{k}));
                
            end
        end
        
        % Total these overlap occurances
        nGoodCommonIsoforms = sum(nIsoformsHere);
        
        % Special case - if two isoforms with same fpkm value (or if
        % useUniformWeights = TRUE, two isoforms) have abutting
        % sequences, then you can have same isospecificity on target
        % isoform, but different number of probes shared between
        % different isoforms.  Something like this:
        % Seq A :         ---------------------------
        % Seq B :     ------------
        % Seq C :                 -----------------------
        % If seqs B and C have same abundance weight then A gets same
        % isospecificity for region overlapped by either.
        % We don't want probes at overlap with A + B (will give
        % background from probe B) so should pick those only for A + C.
        % To fix, remove columns for those 'not overlapping enough'
        % isoforms, and then any rows with one or more '0's should be
        % all 0s.
        
        if any(nGoodCommonIsoforms <= minNumberOfProbes)
            % Which isoforms to exclude by this test
            %excludeIsoforms =  nGoodCommonIsoforms < minNumberOfProbes;
            % Print them here
            
            nIsoformsHere(:, nGoodCommonIsoforms <= minNumberOfProbes) = [];
            
            nIsoformsHere(any(nIsoformsHere == 0, 2), :) = 0;
        end
        
        % Keep probes that appear in all isoforms
        probesToRetain = all(nIsoformsHere, 2);
        
        
        % If target isoform has insufficient probes, move on to next
        % isospecificity value
        if sum(probesToRetain) <= minNumberOfProbes
            
            fprintf(printTarget, 'Only %.d probes from isoform %s at isospecificity %f common across all isoforms\n', ...
                sum(probesToRetain), geneIsoformList(t).id, isospecsHere(onIsospec));
            
            onIsospec = onIsospec + 1;
            regsFound(t,:) = {geneIsoformList(t).name, maxFoundProbes, -5, []};
            continue;
        end
        
        % If sufficient probes on all overlapping isoforms,
        goodProbes = true;

        maxFoundProbes = sum(probesToRetain);
        
        regsFound(t,:) = {geneIsoformList(t).name, maxFoundProbes, 1, probesToRetain};        
        
    end
    
    % Keep probes that appear in all isoforms
    
    if goodProbes
        % Finding good probes worked
        % return approved targetRegions, best isospecificity score for this isoform
        fprintf(printTarget, 'Found %.d common probes for gene %s\n', sum(probesToRetain), geneIsoformList(t).name);
        fprintf(printTarget, 'Isospecificity score = %f for isoform %s\n', isospecsHere(onIsospec), geneIsoformList(t).id);
        
    else
        % Set of probes not found
        
        fprintf(printTarget, 'Insufficient common sequences found for gene %s\n', geneIsoformList(t).name);
        
    end
    
end

% Display how many genes have insufficient probes
statusTally = vertcat(regsFound{:,3}) < 0;
fprintf(printTarget, '%d of %d targets fail common regions selection\n', ...
    sum(statusTally), numel(statusTally));
searchStatusReport(regsFound, printTarget);

% Filter out targetRegions given results here
% Want only genes + isoforms that appear in geneIsoformList
% Want only those targetRegions that are common across all isoforms

allIsoformIds = {tR.id}';
passedIsoformIds = {geneIsoformList(~statusTally).id};

toKeepTRs = tR(ismember(allIsoformIds, passedIsoformIds));
filteredTargetRegions(length(toKeepTRs)) = TargetRegions();

for k = 1:length(toKeepTRs)
    % Keep only the targetRegions selected to retain in above algorithm
    goodTargetRegions = regsFound{strcmp(regsFound(:,1), toKeepTRs(k).geneName), 4};
    %
    %         assignin('base', 'goodTargetRegions', goodTargetRegions);
    %         assignin('base', 'filteredTargetRegions', filteredTargetRegions);
    %
    
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

end


function searchStatusReport(regsFound, varargin)
% Decode status codes in search status table
% Print groups to terminal

printTarget = varargin{1}; % = 1 is to terminal
% = fID to file

% Gene + isoform not found in targetRegions object
% Likely FPKM is below threshold
% -1
failsHere = vertcat(regsFound{:,3}) == -1;
if any(failsHere)
    fprintf(printTarget, 'Gene + isoform does not appear in targetRegions object:\n');
end
listGeneNames(printTarget, regsFound, failsHere);

% Genes that are too short
% -2
failsHere = vertcat(regsFound{:,3}) == -2;
if any(failsHere)
    fprintf(printTarget, 'All isoforms of gene are too short:\n');
end
listGeneNames(printTarget, regsFound, failsHere);

% Gene is possibly long enough, but isoform is not
% -3
failsHere = vertcat(regsFound{:,3}) == -3;
if any(failsHere)
    fprintf(printTarget, 'Target isoform is too short:\n');
end
listGeneNames(printTarget, regsFound, failsHere);

% Insufficient probes at selected isospecificty
% -4
failsHere = vertcat(regsFound{:,3}) == -4;
if any(failsHere)
    fprintf(printTarget, 'Too few probes at single isospecificity:\n');
end
listGeneNames(printTarget, regsFound, failsHere);

% Insufficient probes are common across all isoforms
% -5
failsHere = vertcat(regsFound{:,3}) == -5;
if any(failsHere)
    printf(printTarget, 'Too few probes in common across all isoforms:\n');
end
listGeneNames(printTarget, regsFound, failsHere);

% Passing genes
% 1

    function listGeneNames(printTarget, regsFound, mask)
        for k = 1:length(mask)
            if mask(k)
                fprintf(printTarget, '%s\n', regsFound{k, 1});
            end
        end
        
    end

end


