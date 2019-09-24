function value = CheckParameter(value, type, name)
%--------------------------------------------------------------------------
% CheckParameter(value, type, name)
% This function returns true/false depending on whether the given parameter
% value satisfies the requirements specified by type
%--------------------------------------------------------------------------
% Inputs:
% value: The value of any variable to be checked
%
% type/string or cell array of strings: A flag to specify the check
%   conditions or the set of conditions to test
%   Valid types: 
%           'positive'
%           'nonnegative'
%           'struct'
%           'cell'
%           'string'
%           'boolean'
%           'array' 
%           'filePath'
%           'fileDir'
%           'colormap'
%           'fraction'
%           'handle'
%           'function'
%           'map'
%           'float'
%           'integer'
%           'freeType'
%           'parallel'
% 
% name/string: The name of the parameter to be checked
%--------------------------------------------------------------------------
% Outputs:
% value: The value of the variable to be checked
%
%--------------------------------------------------------------------------
% Jeffrey Moffitt *  &  Alistair Boettiger #
% * jeffmoffitt@gmail.com, # boettiger.alistair@gmail.com
% October 2013 
% Version 1.3
%--------------------------------------------------------------------------
% Creative Commons Liscence
% Attribution-NonCommercial-ShareAlike 3.0 Unported License
% 2013
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Hardcoded Variables
%--------------------------------------------------------------------------
quiet = 0;

%--------------------------------------------------------------------------
% Parse Required Inputs
%--------------------------------------------------------------------------
if nargin < 3
    error('Not enough inputs');
end

if ~ischar(name)
    error('Invalid parameter name');
end

if ~iscell(type)
    temp = type;
    type = cell(1,1);
    type{1} = temp;
    clear temp;
end

%--------------------------------------------------------------------------
% Check Conditions
%--------------------------------------------------------------------------
for i=1:length(type)
    switch type{i}
        case 'positive'
            if value <= 0
                error([name ' is not positive']);
            end
        case 'nonnegative'
            if value < 0
                error([name ' is not nonnegative']);
            end
        case 'struct'
            if ~isstruct(value)
                error([name ' is not a structure']);
            end
        case 'array'
            if length(value)<= 1
                error([name ' is not an array']);
            end
        case 'boolean'
            if ~islogical(value) && value == 1 && value == 0
                error([name ' is not a boolean']);
            end
        case 'string'
            if ~ischar(value)
                error([name ' is not a string']);
            end
        case 'filePath'
            if ~exist(value) == 2 
                error([name ' is not a valid path to a file']);
            end
        case 'fileDir'
            if ~exist(value) == 7
                error([name ' is not a valid path to a directory']);
            end
        case 'cell'
            if ~iscell(value)
                error([name ' is not a cell']);
            end
        case 'colormap'
            if ~(ischar(value) || size(value,2) == 3)
                error([name ' is not a valid colormap']); 
            end
        case 'fraction'
            if ~(value >= 0 && value <= 1)
                error([name ' is not a valid fraction']);
            end
        case 'handle'
            if ~ishandle(value)
               error([name ' is not a handle']);  
            end
        case 'integer'
            if ~( round(value)==value )
               error([name ' is not an integer']);  
            end
        case 'function'
            if ~strcmp(class(value), 'function_handle')
                error([name ' is not a function handle']);
            end
        case 'map'
            if ~strcmp(class(value), 'containers.Map')
                error([name ' is not a containers.Map object']);
            end
        case 'float'
            if ~isfloat(value)
                error([name 'is not a float']); 
            end
        case 'parallel'
            if ~strcmp(class(value), 'parallel.Pool')
                error([name ' is not a parallel.Pool object']);
            end
        case 'freeType'
        otherwise
            error([type{i} ' is not a valid type']);
    end
end
