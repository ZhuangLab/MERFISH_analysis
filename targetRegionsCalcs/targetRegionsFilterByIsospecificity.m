% Given set of targetRegions, figure out which ones to keep 
% The critical filter is based on isospecificity score.
% Score is (effectively) defined as: 
% 
% (FPKM of currentIsoform) / ? (FPKM of isoforms with same sequence)
% 
% There is a slight modification due to an overlap score btw sequences, 
% but this affects on those sequences on 'edges' of common and non-common 
% sequence regions.
% 
% This metric effectively measures the probability that a probe of complementary 
% sequence will bind to the isoform we want and not to the ones we don't. 
% For isoform-specific experiments, we want that probability to be very high.  The 
% current dogma of isospecificity = [0.75, 1] is appropriate.
% Issue is that for 'isoform agnostic' design, there is almost always overlap
% between isoforms.  In that case it is desired to find common regions within
% the provided isoforms and retain those target regions. We want to ignore isoform-
% specific target regions in favor of these common regions. 
% 
% What we *really* want is the most number of transcripts each carrying the most
% number of probes each.  In that scenario, no isoform should have fewer than a 
% threshold number of probes.  Otherwise that isoform would contribute to background, 
% not to signal.  This would be the case if an isoform overlaps some with the 
% 'common' region, but not enough to contribute sufficient probes. In that scenario
% we would want not the most common of common regions, but the common regions that 
% ignore the zone with this partially-overlapping isoform.
% 
% Finding common regions is straightforward.  Given the FPKM table we can find which
% regions overlap with N other isoforms (nb - this would be easier if 
% 'useUniformWeights' flag was TRUE).  We can then filter out all other non-common 
% regions.  Would still design probes against a single chosen isoform, but 
% retained target regions would target other isoforms at a similar rate of 
% probes/transcript. 
% 
% Before this is final, a check is needed that all isoforms would have sufficient
% probes/transcript by targeting the common regions.  If this is not true, then 
% regions that target the short transcript have to be removed.  'Common' regions
% becomes 'common amongst sufficiently long transcripts' regions. 
% 
% Most of the time we expect overlap to look something like this:
% Seq A :    ----------------------------------------------------------
% Seq B :          ------------         ------------------------
% Seq C :              ----------------------------------------
% with some truncation and intron removal, but 'common' regions are OK to
% probe.
%
% The 'common of long transcripts' situation would be something like this:
% Seq A :          --------------------------------------
% Seq B :    ------------------------------------------------------
% Seq C :                                            --------------
% There we'd want to keep common regions between A and B but not C because 
% overlap of A and B and C does not contain sufficient number of probes.
%
% Bad would be this:
% Seq A :          --------------------------------------
% Seq B :    ----------------------
% Seq C :                                          -----------
% where nothing showing up as common would be sufficiently long.  
% In such a case it might be best to design against unique regions in A.
% 
% This implies that 'right' way to choose best-overlapping isoforms is to
% start with lowest isospecificity score in a target transcript. If enough of those
% sequences in the target isoform AND enough targetRegions with same sequence 
% in all other isoforms, then those targetRegions are OK. 
% If not, then try next-most common isospecificity score in target isoform. Do check again.
% Keep working up list of isoform score until test is passed. 
% 
% nb - target isoform can remain as 'most abundant' by FPKM and be fine.
%
% Real way to do this looks at OTTable object in isospecificity calc.  Hack  
% is to operate on targetRegions object between targetRegions and probe generation
% steps in workflow. True hack is to work on CSV output table as in here.


baseFolder = 'C:\Users\ScanningLabAnalysis\Documents\MATLAB\MERFISHProbeDesign\targetRegionsCalcs';

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
    fprintf(1, 'Searching for common targets in gene %s\n', targets(t).geneName);

    rowsHere = ~cellfun(@isempty, strfind(tR{:,1}, targets(t).geneName)) & ...
               ~cellfun(@isempty, strfind(tR{:,2}, targets(t).isoformID));

    % For target isoform, ascending sort of isospecificity score
    isospecsHere = unique(tR{rowsHere, 3});
    onIsospec = 1;

    % While loop on isospec
    % Pick lowest available isospecificity score
    % If none available, then selection fails
    while onIsospec < (length(isospecsHere) + 1)

          isospecRows = ~cellfun(@isempty, strfind(tR{:,1}, targets(t).geneName)) & ...
                        ~cellfun(@isempty, strfind(tR{:,2}, targets(t).isoformID)) & ...
                        tR{:,3} == isospecsHere(onIsospec);

        % Are there enough probes on target isoform with this score?
          if sum(isospecRows) < minNumberOfProbes
              % If no, move to next-lowest isospecificity score
              
              fprintf(1, 'Only %.d probes in isoform %s at isospecificity %f\n', ...
                  sum(isospecRows), targets(t).isoformID, isospecsHere(onIsospec));
              onIsospec = onIsospec + 1;
              continue;
          else
             % If yes, proceed with cross-isoform checks
          end

        % Find sequences with this score
        probeSequences = tR{isospecRows, 5};

        % For each, see how many show up in each of the other isoforms
        isoformsOfThisGene = unique(tR{~cellfun(@isempty, strfind(tR{:,1}, targets(t).geneName)), 2});
        % nb - if there are only 2 isoforms then the above test catches issues
        nIsoformsHere = zeros(sum(isospecRows), length(isoformsOfThisGene));

        for k = 1:length(probeSequences)
           for i = 1:length(isoformsOfThisGene)

              probeFindRows = ~cellfun(@isempty, strfind(tR{:,1}, targets(t).geneName)) & ...
                              ~cellfun(@isempty, strfind(tR{:,2}, isoformsOfThisGene{i}));

              nIsoformsHere(k, i) = any(~cellfun(@isempty, strfind(tR{probeFindRows, 5}, probeSequences{k})));

           end
        end

        % Total these overlap occurances
        nGoodCommonIsoforms = sum(nIsoformsHere);


        % If each overlapped isoform has less than min number of probes, then
        % move to next-lowest isospecificity score for target isoform.
        if any(nGoodCommonIsoforms < minNumberOfProbes)
           onIsospec = onIsospec + 1; 
           fprintf(1, 'Not enough shared probes for %s at isospecificity %f\n', ...
                targets(t).isoformID, isospecsHere(onIsospec));
           
           continue;
        else
           % If sufficient probes on all overlapping isoforms, 
           goodProbes = true;

           % Keep probes that appear in all isoforms 
           probesToRetain = all(nIsoformsHere, 2);

           break;
        end

    end

    if goodProbes
       % Finding good probes worked 
       % return approved targetRegions, best isospecificity score for this isoform 
       fprintf(1, 'Found %.d common probes for gene %s\n', sum(probesToRetain), targets(t).geneName);
       fprintf(1, 'Isospecificity score = %f for isoform %s\n', isospecsHere(onIsospec), targets(t).isoformID);

    else
       % Set of probes not found

       fprintf(1, 'Insufficient common sequences found for gene %s\n', targets(t).geneName);

    end
end



















































