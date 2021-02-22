
gene = 'PPIB';

tt = strfind({targetRegions.geneName}, gene);
ff = [];
for k = 1:length(tt)
    if ~isempty(tt{k})
        ff = [ff; k];
    end
end

lg = cell(length(ff), 2);
for k = 1:numel(ff)
    lg{k,1} = targetRegions(ff(k)).id;
    lg{k, 2} = length(targetRegions(ff(k)).sequence);
end

lg