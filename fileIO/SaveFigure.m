function [savedFilePaths, parameters] = SaveFigure(figHandle, varargin)
% ------------------------------------------------------------------------
% savedFilePaths = SaveFigure(figHandle, varargin)
% This function saves the figure specified by the provided handle to the
% path set by SetFigureSavePath().  
%--------------------------------------------------------------------------
% Necessary Inputs
% figHandle/handle: Handle to a valid figure. 
%
%--------------------------------------------------------------------------
% Outputs
% savePath/string: Path to the saved figure. 
%
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
% 'name'/([]): The name for the figure. The default is to use the figures
%    current name or number.
% 'overwrite'/boolean(false): Determines if the figure will be overwritten. 
%   If the figure already exists, and overwrite is false, a unique index
%   will be appended to the figure, and it will be saved.
% 'formats'/cell array of extensions ({'fig'}): The types of files that
%   will be saved.
% 'appendFormats'/cell array ({}): A list of types to add to the default.
% 'subFolder'/string ({}): A string describing a relative path to a
%    subfolder for saving the specific figure. If empty the figure is saved 
%    in the path set by SetFigureSavePath. 
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 21, 2017
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2018.
%--------------------------------------------------------------------------
% Dependencies:
%   This function uses the export_fig package
%   (http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig)
%   which improves the export quality of certain formats. This package in
%   turn uses ghostscript (http://www.ghostscript.com/) and Xpdf
%   (http://www.foolabs.com/xpdf/).  To turn off this functionality simply
%   pass false with the 'useExportFig' flag.
%   Update: export_fig now comes from https://github.com/ojwoodford/export_fig
%--------------------------------------------------------------------------

% ------------------------------------------------------------------------
% Define default values
% ------------------------------------------------------------------------
global figureSavePath;

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'addGit', 'boolean', false};
defaults(end+1,:) = {'gitState', 'struct', []}; 
defaults(end+1,:) = {'name', 'string', []};
defaults(end+1,:) = {'overwrite', 'boolean', false};
defaults(end+1,:) = {'verbose', 'boolean', true};
defaults(end+1,:) = {'formats', 'cell', {'fig'}};
defaults(end+1,:) = {'appendFormats', 'cell', []};
defaults(end+1,:) = {'useExportFig', 'boolean', true};
defaults(end+1,:) = {'subFolder', 'string', []};
defaults(end+1,:) = {'makeDir', 'boolean', true}; %If a sub folder is specified, is it created if it does not exist
defaults(end+1,:) = {'savePath', 'string', figureSavePath};
defaults(end+1,:) = {'transparent', 'boolean', false};
defaults(end+1,:) = {'closeFig', 'boolean', false};
defaults(end+1,:) = {'saveData', 'boolean', true}; % optionally don't export
% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% ------------------------------------------------------------------------
% Check validity of figure handle
% ------------------------------------------------------------------------
if nargin < 1 || (~(ishghandle(figHandle) && strcmp(get(figHandle, 'type'), 'figure')))
    error('matlabFunctions:invalidHandle', 'Provided figure handle is not valid');
end

% ------------------------------------------------------------------------
% Check for empty savePath --> implies the use of the default
% ------------------------------------------------------------------------
if isempty(parameters.savePath)
    parameters.savePath = figureSavePath;
end

% ------------------------------------------------------------------------
% Check for subFolder and make if necessary
% ------------------------------------------------------------------------
if ~isempty(parameters.subFolder)
    if parameters.subFolder(end) ~= filesep
        parameters.subFolder(end+1) = filesep;
    end
    if ~(exist([parameters.savePath parameters.subFolder]) == 7)
        status = false;
        if parameters.makeDir
            status = mkdir([parameters.savePath parameters.subFolder]);
        end
        if ~status
            error('matlabFunctions:invalidPath', 'Provided subfolder path does not exist or cannot be created');
        end
    end
end

% ------------------------------------------------------------------------
% Check to see if the original directory exists, if makeDir is true
% ------------------------------------------------------------------------
if parameters.makeDir
    if ~(exist(parameters.savePath) == 7)
        mkdir(parameters.savePath);
    end
end

% ------------------------------------------------------------------------
% Compose figure name
% ------------------------------------------------------------------------
if isempty(parameters.name)
    name = get(figHandle, 'Name');
    if isempty(name)
        name = 'unnamed';
    end
else
    name = parameters.name;
end

% ------------------------------------------------------------------------
% Determine types of files to save
% ------------------------------------------------------------------------
if ~isempty(parameters.appendFormats)
    for i=1:length(parameters.appendFormats)
        parameters.formats{end+1} = parameters.appendFormats{i};
    end
    parameters.formats = unique(parameters.formats);
end

% ------------------------------------------------------------------------
% Check if files exist
% ------------------------------------------------------------------------
if ~parameters.overwrite
    done = false;
    count = 0;
    newName = name;
    while ~done
        done = true;
        foundFinalName = true;
        for i=1:length(parameters.formats)
            if exist([parameters.savePath parameters.subFolder newName '.' parameters.formats{i}])
                count = count + 1;
                done = false;
                foundFinalName = false;
                break;
            end
        end
        if ~foundFinalName && count > 0
            newName = [name '(' num2str(count) ')'];
        end
    end
    name = newName;

    if count > 0 && parameters.verbose
        display(['Found existing files. Appending ' num2str(count) ' to all file names.']);
    end
end

% ------------------------------------------------------------------------
% Save files
% ------------------------------------------------------------------------
if parameters.saveData
    for i=1:length(parameters.formats)
        savedFilePath = [parameters.savePath parameters.subFolder name '.' parameters.formats{i}];
        switch parameters.formats{i}
            case {'fig', 'ai'}
                saveas(figHandle, savedFilePath, parameters.formats{i});
            case {'eps', 'png'}
                if parameters.useExportFig
                    try
                        if parameters.transparent
                            export_fig(figHandle, savedFilePath,'-transparent');
                        else
                            export_fig(figHandle, savedFilePath);
                        end
                    catch er
                        warning(er.getReport); 
                        saveas(figHandle, savedFilePath, parameters.formats{i});
                        warning('matlabFunctions:SaveFigure', ['Error using export_fig. Figure saved using matlab defaults. Confirm that you have ghostscript']);
                    end
                else
                    saveas(figHandle, savedFilePath, parameters.formats{i});
                end

            otherwise
                warning('matlabFunctions:SaveFigure', ['Unrecognized/unsupported format ' parameters.formats{i}]);
        end

        if parameters.verbose
            display(['Saved: ' savedFilePath]);
        end
    end
end

if parameters.closeFig
    close(figHandle);
end