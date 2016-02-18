function usedSECDEDcodewords = GenSECDED(numLetters,numDataBits,onBits,varargin)
%--------------------------------------------------------------------------
% Generates SECDEDcodewords based on number of data bits needed
%
%% Parse input
%--------------------------------------------------------------------------
if numLetters<numDataBits
    error(['Data bits cannot be more than total bits']);
end
if numDataBits<onBits
    error(['Used bits cannot be more than data bits']);
end        
%% Parse variable input
%--------------------------------------------------------------------------
if nargin > 3
    if (mod(length(varargin), 2) ~= 0 ),
        error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
    end
    parameterCount = length(varargin)/2;
    for parameterIndex = 1:parameterCount,
        parameterName = varargin{parameterIndex*2 - 1};
        parameterValue = varargin{parameterIndex*2};
        switch parameterName
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end
    
%% Main Function
%-------------------------------------------------------------------------

if numLetters ~= 12
    % numDataBits = 11; %Number of data bits.
    nMsg = 2^numDataBits; %Number of messages from the number of data bits.
    uncodedwords = rem(floor([0:2^numDataBits-1]'*pow2(-(numDataBits-1):0)),2); %Generate a matrix of uncoded words where the rows are the uncoded words
    hammingcodewords = encode(uncodedwords,numLetters-1,numDataBits,'hamming/binary'); %Compute Hamming codewords
    SECDEDcodewords = zeros(size(hammingcodewords,1),numLetters);
    for n = 1:length(hammingcodewords)
        SECDEDcodewords(n,:) = [mod(sum(hammingcodewords(n,:), 2), 2) hammingcodewords(n,:)]; %Compute SECDED codewords
    end               
    if ~isempty(onBits)
        usedSECDEDcodewords = SECDEDcodewords(sum(SECDEDcodewords,2)==onBits,:);
    else
        usedSECDEDcodewords= SECDEDcodewords;
        veryShortOrVeryLong = (sum(SECDEDcodewords,2) <=4 | sum(SECDEDcodewords,2) >= numLetters-4);
        usedSECDEDcodewords(veryShortOrVeryLong,:) = [];
    end
%% 

else
% Make 12-Bit SECDED with 4 on bits

    numDataBits = 11;
    numLetters = 16;

    nMsg = 2^numDataBits; %Number of messages from the number of data bits.
    uncodedwords = rem(floor([0:2^numDataBits-1]'*pow2(-(numDataBits-1):0)),2); %Generate a matrix of uncoded words where the rows are the uncoded words
    hammingcodewords = encode(uncodedwords,numLetters-1,numDataBits,'hamming/binary'); %Compute Hamming codewords
    SECDEDcodewords = zeros(size(hammingcodewords,1),numLetters);
    for n = 1:length(hammingcodewords)
        SECDEDcodewords(n,:) = [mod(sum(hammingcodewords(n,:), 2), 2) hammingcodewords(n,:)]; %Compute SECDED codewords
    end 

    twelveBitSECDED = SECDEDcodewords(1:16:end,1:12);
    [nWords nLetters] = size(twelveBitSECDED);
    hD_short = zeros(nWords);
    for j=1:nWords
    for i=1:nWords  % i = 2;
        hD_short(i,j) = sum(xor(twelveBitSECDED(i,:), twelveBitSECDED(j,:)));
    end
    end

    hD_short(hD_short==0) = inf;
    min(hD_short);
    bits12on4 = twelveBitSECDED(sum(twelveBitSECDED,2)==onBits,:);

    [nWords nLetters] = size(bits12on4);
    hD_4onbits = zeros(nWords); 
    for j=1:nWords
    for i=1:nWords  % i = 2;
        hD_4onbits(i,j) = sum(xor(bits12on4(i,:), bits12on4(j,:)));
    end
    end
    % figure(2); clf; imagesc(hD_4onbits); colorbar;
    usedSECDEDcodewords = bits12on4;
end


