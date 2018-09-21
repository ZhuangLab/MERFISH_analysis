classdef Vrykolakas < handle
% ------------------------------------------------------------------------
% vryObj = Vrykolakas(varargin)
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 21, 2017
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2018.
%--------------------------------------------------------------------------
% This class is used to keep a ssh shell alive

% -------------------------------------------------------------------------
% Define properties
% -------------------------------------------------------------------------
properties
   verbose = true       % Control the verbosity of the class
end

properties (SetAccess=protected)
    % Iternal workings
    statusTimer = [];           % A timer to automatically check status
    timerPeriod = 10*60;        % The period at which the status should be automatically polled
    numStatusChecks = 0;        % The number of status checks
    aliveString = 'I''m alive!' % String to display upon timer completion
end

% -------------------------------------------------------------------------
% Public Methods
% -------------------------------------------------------------------------
methods
    
    % -------------------------------------------------------------------------
    % Define constructor
    % -------------------------------------------------------------------------
    function obj = Vrykolakas(varargin)
        % This class is an object designed to keep an ssh shell alive
        
        % -------------------------------------------------------------------------
        % Parse variable inputs
        % -------------------------------------------------------------------------
        % Define defaults
        defaults = cell(0,3); 
        
        % Properties for automatic updating
        defaults(end+1,:) = {'timerPeriod', ...
            'nonnegative', 10*60};
        defaults(end+1,:) = {'aliveString', ...
            'string', 'I''m alive!'};
        
        % Create parameters structure
        parameters = ParseVariableArguments(varargin, defaults, mfilename);
        
        % Pass fields to new object
        foundFields = fields(parameters);
        for f=1:length(foundFields)
            obj.(foundFields{f}) = parameters.(foundFields{f});
        end
        
        % -------------------------------------------------------------------------
        % Start the timer
        % -------------------------------------------------------------------------
        obj.StartTimer();
        
    end

    % -------------------------------------------------------------------------
    % Start timer for automatic update of status
    % -------------------------------------------------------------------------
    function StartTimer(obj)
        % Start the automatic update timer
        % obj.StartTimer()
        
        obj.statusTimer = timer;
        obj.statusTimer.TimerFcn = @(~,~)disp(obj.aliveString);
        obj.statusTimer.Period = obj.timerPeriod;
        obj.statusTimer.ExecutionMode = 'fixedRate';
        obj.statusTimer.StartDelay = obj.timerPeriod;
        start(obj.statusTimer);
    end
    
    % -------------------------------------------------------------------------
    % Stop the status timer
    % -------------------------------------------------------------------------
    function StopTimer(obj)
        % Stop the automatic update timer
        % obj.StopTimer()
        if ~isempty(obj.statusTimer)
            
            % Stop timer
            stop(obj.statusTimer);
            
        end

    end

    % -------------------------------------------------------------------------
    % Delete 
    % -------------------------------------------------------------------------
    function delete(obj)
        % Delete the object
        % delete(obj)
        
        % Stop the timer
        obj.StopTimer();
        
    end
    
end
    
end
