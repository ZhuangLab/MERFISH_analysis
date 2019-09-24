function defaultParameters = CreateDefaultParameters(defaults)
% ------------------------------------------------------------------------
% defaultParameters = CreateDefaultParameters(defaults)
%--------------------------------------------------------------------------
% Necessary Inputs
% defaults/A Nx3 array containing information on the default parameters 
% for a function.  Each entry must be of the form:
%   {'parameterName', 'parameterType', defaultValue}
%   
%   For example, 
%   {'verbose', 'boolean', true}
%   {'listOption', {'option 1', 'option 2'}, 'option 1'}
%   {'arrayValue', 'array', [1 2 3]}
%
%--------------------------------------------------------------------------
% Outputs
% defaultParameters/ A structure containing a field for each parameter name
% initialized to the specified default value.  
%
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% April 10, 2014
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Create default parameters structure
% -------------------------------------------------------------------------
for i=1:size(defaults, 1)
    defaultParameters.(defaults{i,1}) = defaults{i,3};
end

