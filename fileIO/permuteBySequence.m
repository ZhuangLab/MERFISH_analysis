% permuteBySequence - generate pseudorandom permutation of nBins of
% 1:nOutputsPossible outputs.
% Selection of possible permutation is defined by 
%       mod(sum(sequence), nPossiblePermutations)
% where 
%       nPossiblePermutations = nOutputsPossible!
% 
% sequence should be vector of m integers, where m > nPossiblePermutations
% and m(i) > nPossiblePermutations.  Sum(sequence) should be much greater
% than nPossiblePermutations.
%
% Ex : 
% seq = 'GATTAGGGGTGATAACCAATTGGTACTTTT'; % 30-mer sequence 
% intSeq = uint8(seq);
% returnPermutation = permuteBySequence(intSeq, 4, 3);
%
% returns 
% 
%   [4 1 2]
%


function [returnPermutation, whichPerm] = permuteBySequence(sequence, nOutputsPossible, nBins)

    readoutPerms = perms(1:nOutputsPossible);
    whichPerm = mod(11*sum(sequence), size(readoutPerms, 1)) + 1;
    
    returnPermutation = readoutPerms(whichPerm, 1:nBins);
end