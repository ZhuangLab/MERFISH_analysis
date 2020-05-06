function outputTargetRegionsTable(filePath, tR, varargin)

% Output to table
% Want table with columns:
% [geneName    isoformID    numTargetRegions]
% For all genes that make it to this point
% 
% Should work on any TargetRegions object
%
% Inputs :   filePath - path to export file
%            tR - TargetRegions object to dump to file
%            varargin - option additional argument for filtMeth as first entry



if nargin >= 4
    minRegs = varargin{2};
else
    minRegs = -1;
end

if nargin >= 3
    filtMeth = varargin{1};
else
    filtMeth = 'unspecified';
    minRegs = -1;
end

exportArray = cell(length(tR), 3);

exportArray(:,1) = {tR.geneName}';
exportArray(:,2) = {tR.id}';
exportArray(:,3) = {tR.numRegions}';

exportArray = sortrows(exportArray, 1);

fID = fopen(filePath, 'w+'); 
fprintf(fID, '# TargetRegions export\n');
fprintf(fID, '# Exported on : %s\n', datestr(datetime));
fprintf(fID, '# Filter method : %s\n', filtMeth);
fprintf(fID, '# Min num regions : %d\n', minRegs);
fprintf(fID, 'geneName\tisoformID\tnumRegions\n');
for k = 1:length(exportArray)
    fprintf(fID, '%s\t%s\t%d\n', exportArray{k,:});
end
fclose(fID);