check = zeros(length(NicovichfinalTargetRegions), 1, 'int8');
for k = 1:length(check)
    if fTR(k).numRegions > 0
        check(k) = int8(strcmp(MoffittfinalTargetRegions(k).sequence{1}, NicovichfinalTargetRegions(k).sequence{1}));
    else
        check(k) = -1;
    end
end


%%
check = zeros(length(fTR), 4, 'int8');
for k = 1:size(check, 1)
    
    % Newest headers are equivalent between JM, PRN code
    mPOHeader = strsplit(moffittPossibleOligos_again.headers{k}, ' ');
    mLong = cell2mat(cellfun(@length, mPOHeader, 'UniformOutput', false));
    mwhichLong = find(mLong == max(mLong(:)));
    tPOHeader = strsplit(testPossibleOligos_again.headers{k}, ' ');
    tLong = cell2mat(cellfun(@length, tPOHeader, 'UniformOutput', false));
    twhichLong = find(tLong == max(tLong(:)));
    
    check(k, 1) = int8(strcmp(mPOHeader{mwhichLong}, tPOHeader{twhichLong}));
    
    % Headers are equivalent between consecutive runs of JM code
    mPOHeader = strsplit(moffittPossibleOligos.headers{k}, ' ');
    mLong = cell2mat(cellfun(@length, mPOHeader, 'UniformOutput', false));
    mwhichLong = find(mLong == max(mLong(:)));
    tPOHeader = strsplit(moffittPossibleOligos_again.headers{k}, ' ');
    tLong = cell2mat(cellfun(@length, tPOHeader, 'UniformOutput', false));
    twhichLong = find(tLong == max(tLong(:)));
    
    check(k, 2) = int8(strcmp(mPOHeader{mwhichLong}, tPOHeader{twhichLong}));
    
    % Headers are equivalent between consecutive runs of PRN code
    mPOSeq = strsplit(testPossibleOligos.seqs{k}, ' ');
    mLong = cell2mat(cellfun(@length, mPOSeq, 'UniformOutput', false));
    mwhichLong = find(mLong == max(mLong(:)));
    tPOSeq = strsplit(testPossibleOligos_again.seqs{k}, ' ');
    tLong = cell2mat(cellfun(@length, tPOSeq, 'UniformOutput', false));
    twhichLong = find(tLong == max(tLong(:)));
    
    check(k, 3) = int8(strcmp(mPOSeq{mwhichLong}, tPOSeq{twhichLong}));
    
    % Newest sequences are equivalent between JM, PRN code
    mPOSeq = strsplit(moffittPossibleOligos_again.seqs{k}, ' ');
    mLong = cell2mat(cellfun(@length, mPOSeq, 'UniformOutput', false));
    mwhichLong = find(mLong == max(mLong(:)));
    tPOSeq = strsplit(testPossibleOligos_again.seqs{k}, ' ');
    tLong = cell2mat(cellfun(@length, tPOSeq, 'UniformOutput', false));
    twhichLong = find(tLong == max(tLong(:)));
    
    check(k, 4) = int8(strcmp(mPOSeq{mwhichLong}, tPOSeq{twhichLong}));
end

%%
i = 140;
localGeneName = 'NHSL2';
all(tRegionPull.startPos == tRegionPull_Nicovich.startPos)
all(tRegionPull.startPos == tRegionPull_Moffitt.startPos)


























