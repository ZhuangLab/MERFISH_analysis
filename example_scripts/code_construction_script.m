% Example Code Construction 
% Jeffrey Moffitt
% January 30, 2016
% lmoffitt@mcb.harvard.edu
% -------------------------------------------------------------------------
% Purpose: To illustrate construction of a MHD4 code and a MHD2 code.
% -------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.

%% Generate MHD4 Code
% Generate the Extended Hamming Code
numDataBits = 11;
EHwords = GenerateExtendedHammingWords(numDataBits);

% Find Hamming Weight
hammingWeight = sum(EHwords,2);

% Cut words
MHD4words = EHwords(hammingWeight==4,:);

% Display properties of code
display('-----------------------------------------------------------------');
display(['Constructed ' num2str(size(MHD4words,1)) ' barcodes/words']);
display(['Found the following hamming weights']);
display([num2str(unique(sum(MHD4words,2)))]);

% Check HD
hammingDistance = @(x,y)sum(abs(x-y));
measuredDistances = Inf(size(MHD4words,1), size(MHD4words,1));
for i=1:size(MHD4words,1)
    for j=1:size(MHD4words,1)
        if i==j
            continue;
        else
            measuredDistances(i,j) = hammingDistance(MHD4words(i,:), ...
                MHD4words(j,:));
        end
    end
end
minDist = min(measuredDistances);
display('Found the following minimum HD');
display(num2str(unique(minDist)));
        

%% Generate MHD2 Code
numBits = 14;
onBitInds = nchoosek(1:numBits, 4);
MHD2words = zeros(size(onBitInds,1), numBits);
for i=1:size(MHD2words, 1)
    MHD2words(i,onBitInds(i,:)) = 1;
end

% Display Properties
display('-----------------------------------------------------------------');
display(['Constructed ' num2str(size(MHD2words,1)) ' barcodes/words']);
display(['Found the following hamming weights']);
display([num2str(unique(sum(MHD2words,2)))]);

% Check HD
hammingDistance = @(x,y)sum(abs(x-y));
measuredDistances = Inf(size(MHD2words,1), size(MHD2words,1));
for i=1:size(MHD2words,1)
    for j=1:size(MHD2words,1)
        if i==j
            continue;
        else
            measuredDistances(i,j) = hammingDistance(MHD2words(i,:), ...
                MHD2words(j,:));
        end
    end
end
minDist = min(measuredDistances);
display('Found the following minimum HD');
display(num2str(unique(minDist)));
        

