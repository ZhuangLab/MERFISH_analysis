function [idx, notIdx] = StringFind(cellString,targetString,varargin)
% idx = StringFind(cellString,targetString)
% returns the index in the cell array cellString containing the
% string contained in targetString.  
%
% if both cellString and targetString are cell arrays, StringFind searches
% cellString for each element in targetString and returns the indices in
% cellString where the target element was found.  If these 
% 
% [idS2inS1, idS2notinS1] = StringFind(S1,S2)
% returns idx, the indices in cell string S1 which contain matches (partial
% matches by default as with strfind) to the strings in cell string S2, and
% notIdx, the indices of cell string S2 which are not found in S1.
% 
% 
% To return perfect matches only declare 
%   StringFind(cell,target,'exactly',true)
% 
% Examples
% cellString = {'a','b','c'}
% targetString = {'e','b','d','c'} 
% idx = StringFind(cellString,targetString) 
% ans = [2,3]; 
%
%  Defaults
% 'exactly', 'boolean', false
% 'cellOutput','boolean',false
% 'boolean','boolean',false 
% 
%--------------------------------------------------------------------------
% Alistair Boettiger
% boettiger.alistair@gmail.com
% July 9th, 2014
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'exactly', 'boolean', false};
defaults(end+1,:) = {'cellOutput','boolean',false};
defaults(end+1,:) = {'boolean','boolean',false}; 

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 2
    error('matlabSTORM:invalidArguments', 'a cellString and targetString are required');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
%% Main Function
% -------------------------------------------------------------------------

notIdx = [];

if parameters.exactly % perfect matches
    if ~iscell(targetString)
        if parameters.boolean
            idx = strcmp(cellString,targetString);
        else 
            idx = find(strcmp(cellString,targetString));
        end
    else
        idx = cell(1,length(targetString)); 
        for i=1:length(targetString)
            idx{i} =find(strcmp(cellString,targetString{i}));
        end
        if ~isempty(cellfun(@length,idx))
            if max( cellfun(@length,idx) ) == 1  && ~parameters.cellOutput % convert to cell arrays of scalars to a vector unless requested not to
                notIdx = cellfun(@isempty,idx);
                idx( notIdx ) = []; % remove empties to allow concatinating scalars 
                idx = cell2mat(idx);
            end
        end
    end
    
else % for non perfect matches
    if ~iscell(targetString)
        if parameters.boolean
            idx = ~cellfun(@isempty, strfind(cellString,targetString));
        else
            idx = find(~cellfun(@isempty, strfind(cellString,targetString)));
        end
    else
        idx = cell(1,length(targetString)); 
        for i=1:length(targetString)
            idx{i} = find(~cellfun(@isempty, strfind(cellString,targetString{i})));
        end
        if max( cellfun(@length,idx) ) == 1   && ~parameters.cellOutput  % convert to cell arrays of scalars to a vector  
            notIdx = cellfun(@isempty,idx);
            idx( notIdx ) = []; % remove empties to allow concatinating scalars 
            idx = cell2mat(idx);
        end
    end
end


if parameters.boolean && ~islogical(idx)
   idxLogical = false(length(cellString),1);
   idxLogical(idx) = true;
   idx = idxLogical;
end

