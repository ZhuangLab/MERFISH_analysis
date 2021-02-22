

tt = strfind({targetRegions.geneName}, 'UBC');
ff = [];
for k = 1:length(tt)
    if ~isempty(tt{k})
        ff = [ff; k];
    end
end