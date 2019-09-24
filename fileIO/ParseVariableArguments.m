function parameters = ParseVariableArguments(variableArgs, defaults, mfileName)
% ------------------------------------------------------------------------
% parameters = ParseVariableArguments(variableArgs, defaults)
% This function parses variable arguments based on the default arguments
% and assigns the provided values to the appropriate field in the
% parameters field. 
%--------------------------------------------------------------------------
% Necessary Inputs
% variableArgs/A cell array of pairs of field names, e.g. 'verbose', and
%   values, e.g. true. 
% defaults/A cell array containing information on the default parameters 
% for a function.  Each entry must be of the form:
%   {'parameterName', 'parameterType', defaultValue}
%   The variable arguments are compared against this list. 
%--------------------------------------------------------------------------
% Outputs
% parameters/ A structure containing a field for each parameter name
% initialized to the specified value either in defaults or in the variable 
% arguments.  
%
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% April 10, 2014
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 3
    mFileName = [];
end

% -------------------------------------------------------------------------
% Validate and parse variable arguments
% -------------------------------------------------------------------------
if (mod(length(variableArgs), 2) ~= 0 ),
    error(['Extra parameters must be passed in pairs to ' mfileName '.']);
end
parameterNames = variableArgs(1:2:end);
parameterValues = variableArgs(2:2:end);
numFlags = length(parameterNames);

% -------------------------------------------------------------------------
% Create default parameters structure and return if no variable args were
% provided
% -------------------------------------------------------------------------
defaultParameters = CreateDefaultParameters(defaults);
parameters = defaultParameters;
if isempty(variableArgs)
    return
end

% -------------------------------------------------------------------------
% Overwrite with provided parameters
% -------------------------------------------------------------------------
parametersInd = find(strcmp( parameterNames, 'parameters'));
if ~isempty(parametersInd)
    if length(parametersInd) > 1
        error('matlabSTORM:invalidParameters', 'Only one parameters flag may be used');
    end
    
    parameters = CheckParameter(parameterValues{parametersInd}, ...
        'struct', 'parameters');
    
    parameterNames = parameterNames(...
        [1:(parametersInd-1) (parametersInd+1):numFlags]);
    parameterValues = parameterValues(...
        [1:(parametersInd-1) (parametersInd+1):numFlags]);
end

% ---------------------------------------------------------------------
% Fill in default parameters fields if missing
% ---------------------------------------------------------------------
defaultFields = fieldnames(defaultParameters);
foundFields = fieldnames(parameters);
missingFields = setdiff(defaultFields, foundFields);
for i=1:length(missingFields)
    parameters.(missingFields{i}) = defaultParameters.(missingFields{i});
end
% -------------------------------------------------------------------------
% Check validity of provided values
% -------------------------------------------------------------------------
for i=1:length(parameterNames)
    % ---------------------------------------------------------------------
    % Find corresponding default value
    % ---------------------------------------------------------------------
    ind = find(strcmp(defaults(:,1), parameterNames{i}));
    if isempty(ind)
        error('matlabSTORM:invalidParameters', ... 
            [parameterNames{i} ' is not a valid parameters flag for ' mfileName]);
    end
    % ---------------------------------------------------------------------
    % Check provided value
    % ---------------------------------------------------------------------
    if iscell(defaults{ind, 2})
        parameters.(parameterNames{i}) = ...
            CheckList(parameterValues{i}, defaults{ind, 2}, parameterNames{i});
    else
        parameters.(parameterNames{i}) = ...
            CheckParameter(parameterValues{i}, defaults{ind, 2}, parameterNames{i});
    end
end
