
function waitforfreecpu(varargin)
%--------------------------------------------------------------------------
% waitforfreecpu
%   pauses Matlab and waits until CPU useage drops below the default load
%   before proceeding.  This prevents matlab functional calls from eating
%   up 100% of available cpu and causing system lags/freezes
%
% Windows only.  
%-------------------------------------------------------------------------
% Optional Inputs 
% 'MaxLoad' / double / 75
%               -- Total CPU usage must drop below 
% 'RefreshTime' / double /10
%               -- number of seconds computer waits before checking cpu
%               memory again.  Time for check is not zero, so rapid checks
%               will slow down your process.  
%-------------------------------------------------------------------------
% Outputs - None
% 
%------------------------------------------------------------------------- 
% Examples
% waitforfreecpu('MaxLoad',90)
%   waits for cpu load to drop below 90% before returning control 
% waitforfreecpu('MaxLoad',90,'RefreshTime',10,'verbose',false);

%--------------------------------------------------------------------------
% Default Parameter Values
%--------------------------------------------------------------------------
MaxLoad = 75;
RefreshTime = 10; 
verbose = true;

%--------------------------------------------------------------------------
% Parse Variable Input Arguments
%--------------------------------------------------------------------------

if (mod(length(varargin), 2) ~= 0 ),
    error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
end
parameterCount = length(varargin)/2;

for parameterIndex = 1:parameterCount,
    parameterName = varargin{parameterIndex*2 - 1};
    parameterValue = varargin{parameterIndex*2};
    switch parameterName   
        case 'MaxLoad'
            MaxLoad = CheckParameter(parameterValue, 'positive', 'MaxLoad');
        case 'RefreshTime'
            RefreshTime = CheckParameter(parameterValue, 'positive', 'RefreshTime');
        case 'verbose'
            verbose = CheckParameter(parameterValue,'boolean','verbose');
        otherwise
            error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.' '  See help ' mfilename]);
    end
end



%--------------------------------------------------------------------------
%% Main Function
%--------------------------------------------------------------------------
if ispc
    gotpaused = false;
    [~,load] = dos('wmic cpu get loadpercentage');  % the key line. 
    load = str2double(regexp(load,'[0-9]+','match'));
    while load> MaxLoad
        if verbose
            disp('waiting for free cpu...')
        end
        gotpaused = true;
        pause(RefreshTime);
        [~,load] = dos('wmic cpu get loadpercentage');
        load = str2double(regexp(load,'[0-9]+','match'));
    end
    if verbose && gotpaused;
        disp('now running..');
    end
else
    warning('waitforfreecpu only works for windows.')
end