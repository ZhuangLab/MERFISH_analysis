function value = CheckList(value, list, name)
%--------------------------------------------------------------------------
% CheckList(value, list, name)
% This function returns true/false depending on whether the given parameter
% is a member of list
%--------------------------------------------------------------------------
% Inputs:
% value: The value of any variable to be checked
%
% list/cell array: A list of allowable values
%
% name/string: The name of the parameter to be checked
%--------------------------------------------------------------------------
% Outputs:
% value: The value of the variable to be checked
%
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% October 3, 2012
% jeffmoffitt@gmail.com
%
% Version 1.0
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

if ~iscell(list)
    temp = list;
    list = cell(1,1);
    list{1} = temp;
    clear temp;
end

%--------------------------------------------------------------------------
% Check Conditions
%--------------------------------------------------------------------------
if ~ismember(value, list)
    disp(['valid options for ',name,':'])
    disp(list);
    error([value ' is not a valid option for ' name]);
end
