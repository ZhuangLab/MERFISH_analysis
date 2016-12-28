classdef OTMap < handle
% ------------------------------------------------------------------------
% OTMap = OTMap(initialData, varargin)
% This class provides an interface to a key/value storage system.
%
% This class stores key value pairs as an 2xN array and performs 
% addition operations via intersection with new data. This approach can 
% be faster/more efficient than that utilized by OTMap2 in some situations.
%
% See OTMap2 and OTTable.
%--------------------------------------------------------------------------
% Necessary Inputs
% initialData -- An 2xN array of key (1) and value (2) pairs. Both must be
% doubles. The key values do not need to be unique.
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
    data	% Storage of the key/value pairs
end

% -------------------------------------------------------------------------
% Define methods
% -------------------------------------------------------------------------
methods
    
    % -------------------------------------------------------------------------
    % Define constructor
    % -------------------------------------------------------------------------
    function obj = OTMap(initialData)
		% obj = OTMap(initialData);
		% initialData is a 2XN array of key/value pairs
		
        % -------------------------------------------------------------------------
        % Check input
        % -------------------------------------------------------------------------
        if nargin < 1
            initialData = zeros(2,0);
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
        obj.data = [uniqueKeys; values];
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
        % Find overlapping values, sum where needed, and reassign
        % -------------------------------------------------------------------------         
        obj.data = cat(2,obj.data, newData);
        [uniqueKeys, ~, ic] = unique(obj.data(1,:));
        values = accumarray(ic, obj.data(2,:), [])';
        obj.data = [uniqueKeys; values;];
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
        [~, ia, ib] = intersect(obj.data(1,:), keys);
        
        % -------------------------------------------------------------------------
        % Return values
        % -------------------------------------------------------------------------         
        values(ib) = obj.data(2,ia);
    end
    
    % -------------------------------------------------------------------------
    % Return Table
    % -------------------------------------------------------------------------
    function data = GetTable(obj)
		% Return the internal data array
		% data = obj.GetTable();
		
        data = obj.data;
    end
    
    % -------------------------------------------------------------------------
    % Return keys
    % -------------------------------------------------------------------------
    function data = keys(obj)
		% Return all keys from map
		% keys = obj.keys();
		
        data = obj.data(1,:);
    end

    % -------------------------------------------------------------------------
    % Return Values
    % -------------------------------------------------------------------------
    function data = values(obj)
		% Return all values from map
		% values = obj.values();
		
        data = obj.data(2,:);
    end
    
    % -------------------------------------------------------------------------
    % Return length
    % -------------------------------------------------------------------------
    function numEntries = length(obj)
		% Return the number of key/value pairs
		% numEntries = obj.length();
		% numEntries = length(obj);
		
        numEntries = size(obj.data,2);
    end
end
end