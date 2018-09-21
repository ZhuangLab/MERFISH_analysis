function [savePath, parameters] = SetFigureSavePath(newPath, varargin)
% ------------------------------------------------------------------------
% savePath = SetFigureSavePath(newPath, varargin)
% This function sets the global save path for figures
%--------------------------------------------------------------------------
% Necessary Inputs
% newPath/string: Path to a valid directory. If empty this function
%   returns the current save path
%
%--------------------------------------------------------------------------
% Outputs
% savePath/string: Path to the current save path
%
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
% 'makeDir'/boolean(false): Make the directory if it does not exist.
% 'incrementDir'/boolean(false): If a directory exists, create a new one
%   with a postscript "_N" where N is an integer. 
% 'verbose/boolean(true): Display new directory
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 21, 2017
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2018.

% ------------------------------------------------------------------------
% Define default values
% ------------------------------------------------------------------------
global figureSavePath;

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'makeDir', 'boolean', false}; % Path to codebook
defaults(end+1,:) = {'incrementDir', 'boolean', false}; 
defaults(end+1,:) = {'verbose', 'boolean', true};
defaults(end+1,:) = {'incrementDelimiter', 'string', '_'};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabFunctions:invalidArguments', 'A data path is required.');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% ------------------------------------------------------------------------
% Add filesep if needed
% ------------------------------------------------------------------------
if newPath(end) ~= filesep()
    newPath(end+1) = filesep;
end

% ------------------------------------------------------------------------
% Check for existing folder
% ------------------------------------------------------------------------
doesFolderExist = (exist(newPath) == 7);

% ------------------------------------------------------------------------
% Handle making folder
% ------------------------------------------------------------------------
if ~parameters.makeDir % Only assign to existing folders
    if ~doesFolderExist
        warning('matlabFunctions:invalidPath', ['Provided path does not exist']);
        savePath = [];
    else
        figureSavePath = newPath;
        savePath = newPath;
    end
else % Make the folder if it does not exist
    if ~doesFolderExist % If it does not exist
        status = mkdir(newPath);
        if ~status %If folder construction failed
            error('matlabFunctions:invalidPath', ['Unable to create ' newPath]);
        else
            figureSavePath = newPath;
            savePath = figureSavePath;
        end
    else % If it exists
        if ~parameters.incrementDir % And the user does not want to increment
            figureSavePath = newPath;
            savePath = figureSavePath;
        else % Handle the case in which the folder exists, but the user wants a new unique folder
            count = 1;
            incrementedPath = newPath;
            while (exist(incrementedPath) == 7) % Keep looping until a unique folder name is created
                incrementedPath = [newPath(1:(end-1)) parameters.incrementDelimiter num2str(count) filesep];
                count = count + 1;
            end
            status = mkdir(incrementedPath);
            if ~status
                error('matlabFunctions:invalidPath', ['Unable to create ' incrementedPath]);
            else
                figureSavePath = incrementedPath;
                savePath = figureSavePath;
            end
        end
    end
end

% ------------------------------------------------------------------------
% Display new path
% ------------------------------------------------------------------------
if parameters.verbose
    display(['Current Save Figure Path: ' figureSavePath]);
end