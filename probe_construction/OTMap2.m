classdef OTMap2 < handle
% ------------------------------------------------------------------------
% OTMap = OTMap2(initialData, varargin)
% This class provides an interface to a key/value storage system.
%
% This class stores key value pairs using a containers.Map object. This 
% approach can be faster/more efficient than a simple 2XN array in some
% situations.
%
% See OTMap and OTTable.
%--------------------------------------------------------------------------
% Necessary Inputs
% initialData -- An 2xN array of key (1) and value (2) pairs. Both must be
% doubles. The key values do not need to be unique
%--------------------------------------------------------------------------
% Methods
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
% None
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% April 20, 2015
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Define properties
% -------------------------------------------------------------------------
properties (GetAccess=private)
    data	% Storge of the key/value pairs
end

% -------------------------------------------------------------------------
% Define methods
% -------------------------------------------------------------------------
methods
    % -------------------------------------------------------------------------
    % Define constructor
    % -------------------------------------------------------------------------
    function obj = OTMap2(initialData)
		% obj = OTMap2(initialData);
		% initialData is a 2XN array of key/value pair
	
        % -------------------------------------------------------------------------
        % Check input
        % -------------------------------------------------------------------------
        if nargin < 1 || isempty(initialData)
            obj.data = containers.Map('KeyType','double','ValueType', 'double');
            return
        end
        if ~isa(initialData, 'double') || size(initialData,1) ~= 2
            error('matlabFunctions:invalidArguments', ...
                ['initialData must be a double of size 2xN']);
        end
        
        % -------------------------------------------------------------------------
        % Find unique keys and accumulate values
        % -------------------------------------------------------------------------         
        [uniqueKeys, ~, ic] = unique(initialData(1,:));
        values = accumarray(ic, initialData(2,:), [])';
        
        % -------------------------------------------------------------------------
        % Set data
        % -------------------------------------------------------------------------         
        obj.data = containers.Map(uniqueKeys, values);
    end
    
    % -------------------------------------------------------------------------
    % AddToMap
    % -------------------------------------------------------------------------
    function AddToMap(obj, newData)
		% Add additional key/value pairs to an existing class
		% obj.AddToMap(newData); 

        % -------------------------------------------------------------------------
        % Check data
        % -------------------------------------------------------------------------         
        if ~isa(newData, 'double') || size(newData,1) ~=2
            error('Invalid data format');
        end
        
        % -------------------------------------------------------------------------
        % Find unique entries and sum as needed
        % -------------------------------------------------------------------------
        [keys, ~, ic] = unique(newData(1,:));
        values = accumarray(ic, newData(2,:), [])';
        
        % -------------------------------------------------------------------------
        % Find overlapping values, sum where needed, and reassign
        % -------------------------------------------------------------------------
        validKeys = obj.data.isKey(num2cell(keys));
        if any(validKeys)
            updatedValues = obj.data.values(num2cell(keys(validKeys)));
            updatedValues = [updatedValues{:}] + values(validKeys);
            oldKeys = keys(validKeys);
        else
            oldKeys = [];
            updatedValues = [];
        end
        newKeys = keys(~validKeys);
        newValues = values(~validKeys);
        newMap = containers.Map([oldKeys newKeys], [updatedValues newValues]);
        
        obj.data = [obj.data; newMap];
    end
        
    % -------------------------------------------------------------------------
    % Return key values
    % -------------------------------------------------------------------------
    function values = GetValues(obj, keys)
		% Return values for specified keys
		% values = obj.GetValues(keys)
    
		% -------------------------------------------------------------------------
        % Prepare output
        % -------------------------------------------------------------------------         
        values = zeros(1, length(keys));
        
        % -------------------------------------------------------------------------
        % Find keys
        % -------------------------------------------------------------------------         
        validKeys = obj.data.isKey(num2cell(keys));
        if any(validKeys)
            values(validKeys) = cell2mat(obj.data.values(num2cell(keys(validKeys))));
        end
    end
    
    % -------------------------------------------------------------------------
    % Return Table
    % -------------------------------------------------------------------------
    function data = GetTable(obj)
		% Return values of the internal containers.Map object as a 2XN array
		% data = obj.GetTable();
    
		keys = cell2mat(obj.data.keys());
        values = cell2mat(obj.data.values());
        data = [keys; values;];
    end
    
    % -------------------------------------------------------------------------
    % Return keys
    % -------------------------------------------------------------------------
    function data = keys(obj)
		% Return all keys from map
		% keys = obj.keys();

        data = cell2mat(obj.data.keys());
    end

    % -------------------------------------------------------------------------
    % Return Values
    % -------------------------------------------------------------------------
    function data = values(obj)
		% Return all values from map
		% values = obj.values();

        data = cell2mat(obj.data.values());
    end
    
    % -------------------------------------------------------------------------
    % Return length
    % -------------------------------------------------------------------------
    function numEntries = length(obj)
		% Return the number of key/value pairs
		% numEntries = obj.length();
		% numEntries = length(obj);

		numEntries = obj.data.Count;
    end
    
end
end