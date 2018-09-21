function [sbsInfo, parameters] = LoadSplitByteStreamHeader(filePath)
% ------------------------------------------------------------------------
% sbsInfo = LoadSplitByteStreamHeader(filePath)
% This function loads a split bytestream header file/structure. See
% SaveSplitByteStream();
%
%--------------------------------------------------------------------------
% Necessary Inputs
% filePath/A string to a bytestream file
%--------------------------------------------------------------------------
% Outputs
% sbsInfo/The reconstituted byteStream header. 
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 21, 2017
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2018.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Parse necessary inputs
%--------------------------------------------------------------------------
if any(nargin < 1)
    error('matlabFunctions:invalidArguments', 'A file path is required');
end

%--------------------------------------------------------------------------
% Load split byte stream info
%--------------------------------------------------------------------------
sbsInfo = LoadByteStream(filePath);

%--------------------------------------------------------------------------
% Check validity
%--------------------------------------------------------------------------
if ~isstruct(sbsInfo) || ~isfield(sbsInfo, 'isByteStream') || ~strcmp(sbsInfo.isByteStream, 'Y')
    error('matlabFunctions:invalidFile', 'The specified matb file does not represent a split byte stream');
end
