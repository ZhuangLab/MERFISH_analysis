geneName = 'LHX6';

fpkmIDs = fpkm{~cellfun(@isempty, strfind(fpkm{:,5}, geneName)), 1};
fpkmVals = fpkm{~cellfun(@isempty, strfind(fpkm{:,5}, geneName)), 10};
isoformsOfThisGene = unique(tR{~cellfun(@isempty, strfind(tR{:,1}, geneName)), 2});

tt = zeros(length(fpkmIDs), 1);
for k = 1:length(fpkmIDs)
tt(k) = any(~cellfun(@isempty, strfind(isoformsOfThisGene, fpkmIDs{k})));
end

excludedIsoformFPKMs = fpkmVals(~tt)