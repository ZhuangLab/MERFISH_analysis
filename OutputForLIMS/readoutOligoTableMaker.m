% Make some dummy output for a readoutOligo tracking file
% Run convertMERFISHOutputToLIMS.m first

% WorkingStockName | SequenceName | Sequence | Fluor | \lambda | ConcentrateName

outputFile = 'ReadoutOligosTable.csv';

% From nameConversion file:
% WorkingStockName = nameConversion{2}_\lambda.00
% SequenceName = nameConversion{2}
% Sequence = nameConversion{3}
% Fluor = AlexaFluor+nameConversion{5}
% \lambda = nameConversion{5}
% ConcentrateName = nameConversion{2}_\lambda

fID = fopen(outputFile, 'w+');

fprintf(fID, '# ReadoutOligosTable\n');
fprintf(fID, '# Generated %s\n', datetime);
fprintf(fID, 'ReadoutOligoName,SequenceName,Sequence,Fluor,emission_wavelength,ConcentrateName,tube_name,dateCreated\n');

onLine = 2;

for k = 2:(length(nameConversion{1})-1)
    
    fluorsHere = strtrim(strsplit(nameConversion{5}{onLine}, ','));
    
    for f = 1:length(fluorsHere)
    
        fprintf(fID, '%s,%s,%s,%s,%s,%s,%d,%s\n', ...
            strcat(nameConversion{2}{onLine}, '_', fluorsHere{f}, '.00'), ...
            nameConversion{2}{onLine}, ...
            nameConversion{3}{onLine}, ...
            strcat('AlexaFluor', fluorsHere{f}), ...
            strcat(fluorsHere{f}), ...
            strcat(nameConversion{2}{onLine}, '_', fluorsHere{f}),...
            0,...
            (datetime('today', 'format', 'yyyyMMdd') + caldays(randi([-20, 0], 1, 1))));
        
    end
    
    onLine = onLine + 1;
    
end

fclose(fID);
