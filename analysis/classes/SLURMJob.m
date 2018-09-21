classdef SLURMJob < handle
% ------------------------------------------------------------------------
% jobObj = SLURMJob(varargin)
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 21, 2017
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2018.
% -------------------------------------------------------------------------
% This class is a wrapper around a job submitted via SLURM to a cluster

% -------------------------------------------------------------------------
% Define properties
% -------------------------------------------------------------------------
properties
   verbose = true       % Control the verbosity of the class
   veryverbose = false  % Display all interactions with slurm
end

properties (SetAccess=protected)
    % Properties that describe the job
    jobID = ''          % The job id assigned at submission
    user = ''           % The username associated with the job
    name = 'slurm'      % The name associated with the job
    
    % Job execution properties
    memoryLimit = 6000          % Specify the amount of memory for the job
    timeLimit = 60              % Time allocated for job in minutes
    partition = 'general'       % The partition to which to submit the job
    ntasks = 1                  % The number of tasks to assign to a job
    cpusPerTask = 1             % The number of cpus per task
    constraint = ''             % Job constraint
    numberResubmit = 0          % The number of times a failed job will be resubmitted before registering as a failure
    exclude = ''                % Nodes to exclude
    
    % Reporting
    errorLog = ''               % Location of the error log file
    outputLog = ''              % Location of the output file
    openMode = 'append'         % Overwrite or append to existing log/output files ('append' or 'truncate')
    mailType = 'NONE'           % When emails will be sent
    mailUser = ''               % The email to which mail will be sent
    
    % SLURM commands
    scriptPath = ''             % The path where where the slurm script will be written
    uniqueID = ''               % A unique ID associated with each job
    scriptText = {};            % The text of the script
    scriptName = '';            % Script name
    
    % External job complete function
    completeFcn = [];           % A handle to an external function that determines if the job is actually complete
    preCheck = false;           % Check the job completion status before submitting?
    
    % Status
    submitted = false;                  % Has the job been submitted?
    numFailures = 0;                    % The number of failures
    failed = false;                     % Has the job failed?
    completed = false;                  % Has the job completed?
    user_canceled = false;              % Has the user canceled the job?
    state = '';                         % The job state
    duration = 0;                       % The duration of job in seconds
    startTime = '';                     % The start time (when the job was submitted first)
    endTime = '';                       % The end time (when the job finally registered as completed)
    history = cell(0,2);                % A history of the job state and time 
    command_line_history = cell(0,4);   % A history of all commands sent to command line and responses
    
    % Iternal workings
    statusTimer = [];           % A timer to automatically check status
    timerPeriod = 0;            % The period at which the status should be automatically polled
    numStatusChecks = 0;        % The number of status checks
    jobTimer                    % The timer used to keep track of job duration
end

% -------------------------------------------------------------------------
% Define events
% -------------------------------------------------------------------------
events
    JobSubmitted                % Signal the job as submitted
    JobComplete                 % Signal the job as completed
    JobFailed                   % Signal the job as failed
end

% -------------------------------------------------------------------------
% Public Methods
% -------------------------------------------------------------------------
methods
    
    % -------------------------------------------------------------------------
    % Define constructor
    % -------------------------------------------------------------------------
    function obj = SLURMJob(scriptText, varargin)
        % This class is a wrapper for a single job submitted to SLURM
        
        % -------------------------------------------------------------------------
        % Parse variable inputs
        % -------------------------------------------------------------------------
        % Define defaults
        defaults = cell(0,3); 
        
        % Properties for job execution
        defaults(end+1,:) = {'name', ...                  
            'string', 'slurm'};
        defaults(end+1,:) = {'memoryLimit', ...
            'positive', 6000}; % In MB
        defaults(end+1,:) = {'timeLimit', ...
            'positive', 60}; % In minutes
        defaults(end+1,:) = {'partition', ...
            'string', 'general'};
        defaults(end+1,:) = {'ntasks', ...
            'positive', 1};
        defaults(end+1,:) = {'cpusPerTask', ...
            'positive', 1};
        defaults(end+1,:) = {'constraint', ...
            'string', ''};
        defaults(end+1,:) = {'exclude', ...
            'string', ''};
        defaults(end+1,:) = {'numberResubmit', ...
            'nonnegative', 0};
        
        % Properties for reporting
        defaults(end+1,:) = {'errorLog', ...
            'string', ''};
        defaults(end+1,:) = {'outputLog', ...
            'string', ''};
        defaults(end+1,:) = {'openMode', ...
            {'append', 'truncate'}, 'append'};
        defaults(end+1,:) = {'mailType', ...
            'string', 'NONE'};
        defaults(end+1,:) = {'mailUser', ...
            'string', ''};
        
        % Script properties
        defaults(end+1,:) = {'scriptPath', ...
            'fileDir', ''};
        defaults(end+1,:) = {'scriptName', ...
            'string', ''};
       
        % Properties for automatic updating
        defaults(end+1,:) = {'timerPeriod', ...
            'nonnegative', 0};
        
        % Properties for debugging/status updates
        defaults(end+1,:) = {'verbose', ...
            'boolean', true};
        defaults(end+1,:) = {'veryverbose', ...
            'boolean', false};
        
        % Properties for confirming completion
        defaults(end+1,:) = {'completeFcn', ...
            'function', @(~)true};
        defaults(end+1,:) = {'preCheck', ...
            'boolean', true};
        
        % Create parameters structure
        parameters = ParseVariableArguments(varargin, defaults, mfilename);
        
        % Pass fields to new object
        foundFields = fields(parameters);
        for f=1:length(foundFields)
            obj.(foundFields{f}) = parameters.(foundFields{f});
        end
        
        % Return empty object          
        if nargin < 1
            return;
        end
        
        % -------------------------------------------------------------------------
        % Check necessary input
        % -------------------------------------------------------------------------
        if isempty(scriptText) || ~iscell(scriptText) 
            error('matlabFunctions:invalidArgument', 'A set of commands for the job must be provided.');
        end
        
        % -------------------------------------------------------------------------
        % Generate a unique ID for the job
        % -------------------------------------------------------------------------
        obj.uniqueID = char(java.util.UUID.randomUUID);
        
        % -------------------------------------------------------------------------
        % Parse script text (in future versions)
        % -------------------------------------------------------------------------
        obj.scriptText = scriptText;        
        
        % -------------------------------------------------------------------------
        % Write the associated slurm script
        % -------------------------------------------------------------------------
        obj.WriteScript();
        
    end
    
    % -------------------------------------------------------------------------
    % Write the slurm script associated with the job
    % -------------------------------------------------------------------------
    function WriteScript(obj, varargin)
        % Write the slurm script for this job
        % obj.WriteScript();
        
        % Define the script filename
        if isempty(obj.scriptName)
            obj.scriptName = [obj.uniqueID '.slurm'];
        end
        
        % Make directory if required
        if ~exist(obj.scriptPath, 'dir')
            mkdir(obj.scriptPath);
        end
        
        % Attempt to open file
        fid = fopen([obj.scriptPath obj.scriptName], 'W');
        
        if fid < 1
            error('matlabFunctions:invalidFile', 'Unable to open script file');
        end
        
        % Write the file header
        fprintf(fid, '%s\n', '#!/bin/bash');
        fprintf(fid, '#SBATCH -n %s\n', num2str(obj.ntasks)); % Number of tasks
        fprintf(fid, '#SBATCH -N %s\n', '1');  % Place all tasks on a single core
        fprintf(fid, '#SBATCH -c %s\n', num2str(obj.cpusPerTask)); % Number of cpus per task
        fprintf(fid, '#SBATCH -p %s\n', obj.partition); % Request a partition
        fprintf(fid, '#SBATCH -t %s\n', num2str(obj.timeLimit)); % Time limit in minutes
        fprintf(fid, '#SBATCH --mem=%s\n', num2str(obj.memoryLimit)); % Memory liimt in MB
        if ~isempty(obj.mailUser)
            fprintf(fid, '#SBATCH --mail-user=%s\n', obj.mailUser);
            if ~isempty(obj.mailType)
                fprintf(fid, '#SBATCH --mail-type=%s\n', obj.mailType);
            end
        end        
        if ~isempty(obj.outputLog)
            fprintf(fid, '#SBATCH -o %s\n', obj.outputLog);
        end
        if ~isempty(obj.errorLog)
            fprintf(fid, '#SBATCH -e %s\n', obj.errorLog);
        end
        if ~isempty(obj.openMode)
            fprintf(fid, '#SBATCH --open-mode=%s\n', obj.openMode);
        end
        if ~isempty(obj.constraint)
            fprintf(fid,'#SBATCH --constraint=%s\n', obj.constraint);
        end
        if ~isempty(obj.exclude)
            fprintf(fid,'#SBATCH --exclude=%s\n', obj.exclude);
        end
        
        % Write the file contents
        fprintf(fid, '\n');
        if ischar(obj.scriptText)
            fprintf(fid, '%s\n', obj.scriptText);
        else
            for e=1:length(obj.scriptText)
                fprintf(fid, '%s\n', obj.scriptText{e});
            end
        end
        
        % Close fid
        fclose(fid);

    end
    
    % -------------------------------------------------------------------------
    % Resubmit a job
    % -------------------------------------------------------------------------
    function Resubmit(obj, varargin)
        % Resubmit a job that has already been submitted
        % success = obj.Resubmit();
        
        % Clear previous submission
        obj.jobID = '';
        
        % Display progress
        if obj.verbose
            disp(['Resubmitting job: ' obj.name]);
        end
        
        % Submit
        obj.Submit();
    end
    
    % -------------------------------------------------------------------------
    % Resubmit a job
    % -------------------------------------------------------------------------
    function Requeue(obj, varargin)
        % Requeue a job
        
        %%% UNDER CONSTRUCTION

    end
    
    
    % -------------------------------------------------------------------------
    % Submit the job to the queue
    % -------------------------------------------------------------------------
    function success = Submit(obj, varargin)
        % Submit the job to the queue
        % success = obj.Submit();
        
        % -------------------------------------------------------------------------
        % Check to see if the job has already been submitted
        % -------------------------------------------------------------------------
        if ~isempty(obj.jobID)
           success = false;
           warning('This job has already been submitted. Use Resubmit to submit again.');
           return;
        end
        
        % -------------------------------------------------------------------------
        % Check job completion status (if requested) prior to submission
        % -------------------------------------------------------------------------
        if obj.preCheck
            if obj.completeFcn()
                
                % Update state
                obj.state = 'COMPLETED';
                
                % Update history
                obj.history(end+1,:) = {'ALREADY COMPLETED', datetime('now')};

                % Mark as submitted
                obj.submitted = true;
                
                % Start job timer
                obj.jobTimer = tic;
                
                % Record start time
                obj.startTime = datetime('now');
                
                % Send submitted signal
                notify(obj, 'JobSubmitted');
                                
                % Display progress
                if obj.verbose
                    disp(['Completed: ' obj.name]);
                end

                % Mark as complete
                obj.completed = true;

                % Stop timer
                if ~isempty(obj.statusTimer)
                    stop(obj.statusTimer);
                end

                % Update duration
                obj.duration = toc(obj.jobTimer);

                % Send notification
                notify(obj, 'JobComplete'); 
                
                % Leave this function
                return;
                
            end
        end
        
        % Display progress
        if obj.verbose
            disp(['Submitting job: ' obj.name]);
        end
        
        % -------------------------------------------------------------------------
        % Compose job command
        % -------------------------------------------------------------------------
        jobCommand = ['sbatch ' obj.scriptPath obj.scriptName];
        
        % -------------------------------------------------------------------------
        % Submit
        % -------------------------------------------------------------------------
        maxAttempts = 3;
        for i=1:maxAttempts
            try
                % Submit command to command line
                [status, cmdout] = system(jobCommand);
                success = ~status;
            catch
                success = false;
                cmdout = '';
                status = 'FAILED SYSTEM CALL';
            end

            % Archive the submission and response
            obj.command_line_history(end+1,:) = {jobCommand, status, cmdout, datetime('now')};

            % Display to user
            if obj.veryverbose
                display(['Submitted: ' jobCommand]);
                display(['Returned: ' cmdout]);
            end
            
            % Leave the attempt loop if the command was submitted properly
            if success
                break;
            end
            
            % Display the problem
            if obj.verbose
                disp(['Encountered a problem on job submission. Attempt ' num2str(i)]);
            end
            
            % Pause if there was a problem
            pause(5);
        end
        
        % -------------------------------------------------------------------------
        % Parse output to get job information
        % -------------------------------------------------------------------------
        if success
            structuredOutput = regexp(cmdout, 'Submitted batch job (?<jobID>[0-9]+)', 'names');
            if ~isempty(structuredOutput) && ~isempty(structuredOutput(1).jobID)
                obj.jobID = str2num(structuredOutput(1).jobID);
            else
                success = false;
            end
        end
        
        % -------------------------------------------------------------------------
        % Update job status
        % -------------------------------------------------------------------------
        if success
            % Mark as submitted
            obj.submitted = true;
            % Start job timer
            obj.jobTimer = tic;
            % Record start time
            obj.startTime = datetime('now');
            % Add to job history
            obj.history(end+1,:) = {'SUBMITTED', datetime('now')};
            % Send submitted signal
            notify(obj, 'JobSubmitted');
        end

        % -------------------------------------------------------------------------
        % Start job status timer
        % -------------------------------------------------------------------------
        if obj.timerPeriod > 0
            obj.StartTimer();
        end
        
    end

    % -------------------------------------------------------------------------
    % Start timer for automatic update of status
    % -------------------------------------------------------------------------
    function StartTimer(obj)
        % Start the automatic update timer
        % obj.StartTimer()
        
        % Handle the case that the status timer already exists
        if ~isempty(obj.statusTimer)
            obj.StopTimer();
        end
        
        % Configure the timer and start it
        obj.statusTimer = timer;
        obj.statusTimer.TimerFcn = @(~,~)obj.UpdateStatus();
        obj.statusTimer.Period = obj.timerPeriod;
        obj.statusTimer.ExecutionMode = 'fixedSpacing';
        obj.statusTimer.StartDelay = obj.timerPeriod;
        obj.statusTimer.Name = ['Timer: ' obj.name];
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
            
            % Update history
            obj.history(end+1,:) = {'TIMER_STOPPED', datetime('now')};
        else
        
            % Update history
            obj.history(end+1,:) = {'TIMER_ALREADY_STOPPED', datetime('now')};
        end

    end
    
    % -------------------------------------------------------------------------
    % Check status
    % -------------------------------------------------------------------------
    function UpdateStatus(obj)
        % Update the status of a submitted job
        % obj.UpdateStatus();
        
        % Compose command
        command = ['sacct -j ' num2str(obj.jobID) ...
            ' --format=JobID,State' ...
            ' -n' ... % No header
            ' -p' ... % Parsable output
            ' --delimiter=,']; % Delimiter
        % Define the attempt loop
        maxAttempt = 3;
        for m=1:maxAttempt
            % Submit  command
            try
                % Submit to command line
                [status, cmdout] = system(command);

                % Archive the submission and response
                obj.command_line_history(end+1,:) = {command, status, cmdout, datetime('now')};

                if obj.veryverbose
                    disp(['Issued: ' command]);
                    disp(['Returned: ' cmdout]);
                end
                success = ~status;
            catch
                success = false;
            end

            % Update number of checks
            obj.numStatusChecks = obj.numStatusChecks + 1;

            % Handle whether or not the queury was successful
            if ~success
                if obj.verbose
                    disp(['Error in status query for job: ' obj.name]);
                end

                % Mark state as corrupted to handle this situation
                state = 'CORRUPTED';
            else
                % Parse the state provided by SLURM
                output = regexp(cmdout, ...
                    '(?<jobID>[0-9]+),(?<state>\w+)', 'names');

                % Handle the case that the returned job id does not match
                if str2num(output.jobID) ~= obj.jobID
                    if obj.verbose
                        disp(['Error in status query for job: ' obj.name]);
                        disp(['Requested status of job id ' num2str(obj.jobID) ' and received a response for job id ' num2str(output.jobID)]);
                    end
                    state = 'CORRUPTED STATUS REQUEST';
                    success = false; % Mark the status request as a failure
                else % Handle the case that everything is correct
                    state = output.state;
                end                        
            end
            
            % Break the status inquiry loop if it was a success
            if success
                break;
            end
            
        end
                
        % Handle the state
        obj.HandleState(state);
        
    end

    % -------------------------------------------------------------------------
    % Cancel job
    % -------------------------------------------------------------------------
    function success = Cancel(obj, varargin)
        % Cancel the job
        % obj.Cancel();
        
        % Check to see if the job was submitted and is running
        if obj.submitted && ~obj.completed
        
            % Compose the command
            command = ['scancel ' num2str(obj.jobID)];

            % Send the command
            [status, cmdout] = system(command);
            
            % Archive the submission and response
            obj.command_line_history(end+1,:) = {command, status, cmdout, datetime('now')};

            % Display to user
            if obj.veryverbose
                disp(['Issued: ' command]);
                disp(['Returned: ' cmdout]);
            end        
            success = ~status;
            
            % Handle successful cancel
            if success
                % Mark the job as canceled by the user
                obj.user_canceled = true;
                % Display status of cancel if appropriate
                if obj.verbose
                    disp(['Canceled ' obj.name ': ' num2str(obj.jobID)]);
                end
            else
                % Display status of cancel if appropriate
                if obj.verbose
                        disp(['Failed to canceled ' obj.name ': ' num2str(obj.jobID)]);
                end
            end
        end
    end
    
    % -------------------------------------------------------------------------
    % Display file
    % -------------------------------------------------------------------------
    function Display(obj, fileToDisplay)
        % Display a text file associated with this class
        % obj.Display('script'); % The slurm script
        % obj.Display('error'); % The error file
        % obj.Display('output'); % The output file
        
        % Structure output
        PageBreak();
        
        % Switch to desired file type
        switch fileToDisplay
            case {'script', 'scriptFile', 's'}
                disp(['Displaying script']);
                filePath = [obj.scriptPath obj.scriptName];
            case {'output', 'outputFile', 'out', 'o'}
                disp(['Display output']);
                filePath = [obj.outputLog];
            case {'error', 'errorFile', 'err', 'e'}
                disp(['Display error']);
                filePath = [obj.errorLog];
            otherwise
                error('matlabFunctions:invalidArguments', 'Unknown file type');
        end
        
        % Display file name
        disp(['File: ' filePath]);
        PageBreak();
        
        % Open file
        fid = fopen(filePath, 'r');
        
        if fid < 0
            error('matlabFunctions:invalidFile', 'Could not open the file');
        end
        
        % Load and print file line by line
        while ~feof(fid)
            line = fgetl(fid);
            if line ~= -1
                disp(line);
            end
        end
        
        % Close file id
        fclose(fid);
    end
    
end

methods (Access=protected)
    % -------------------------------------------------------------------------
    % Handle the state update of the object
    % -------------------------------------------------------------------------
    function HandleState(obj, state)
        % Handle the state of the object
                
        % Update job state
        obj.state = state;
        
        % Update history
        obj.history(end+1,:) = {state, datetime('now')};
        
        % Switch based on state
        switch state
            case 'RUNNING'
                
                % If an output file has been requested check to see if it
                % was created
                if ~isempty(obj.outputLog)
                    if ~exist(obj.outputLog, 'file') % If it does not exist, there is a problem
                        obj.Cancel(); % Cancel the job
                        obj.HandleState('CORRUPTED'); % Mark the job as corrupted
                        return; 
                    end
                end            
            
            case {'FAILED', 'TIMEOUT', 'CORRUPTED'} % Mark as failed
                % Display progress
                if obj.verbose
                    disp(['Failed: ' obj.name]);
                end

                % Update number of failures
                obj.numFailures = obj.numFailures + 1;

                % Stop update timer
                obj.StopTimer();

                % Check for requested resubmit
                if obj.numFailures <= obj.numberResubmit
                    obj.Resubmit();
                else
                    % Mark job as failed
                    obj.failed = true;
                    % Update duration
                    obj.duration = toc(obj.jobTimer);
                    % Mark end time
                    obj.endTime = datetime('now');
                    % Send notification
                    notify(obj, 'JobFailed');
                end 
                
            case 'COMPLETED' % Job completed
                if obj.completeFcn() % Check to see if valid completion
                    % Display progress
                    if obj.verbose
                        disp(['Completed: ' obj.name]);
                    end

                    % Mark as complete
                    obj.completed = true;

                    % Stop timer
                    obj.StopTimer();
                    
                    % Update duration
                    obj.duration = toc(obj.jobTimer);

                    % Send notificaty
                    notify(obj, 'JobComplete'); 
                else % Handle the case that it failed the complete function
                    if obj.verbose
                        disp(['Failed complete function: ' obj.name]);
                        obj.history(end+1, :) = {'FCN_FAILED', datetime('now')};
                    end
                    obj.HandleState('FAILED');
                end

            case 'CANCELED' % Job canceled
                
                % Handle the case that this cancel was requested by the
                % user
                if obj.user_canceled
                    % Display progress
                    if obj.verbose
                        disp(['Canceled: ' obj.name]);
                    end

                    % Stop timer
                    obj.StopTimer();

                    % Update duration
                    obj.duration = toc(obj.jobTimer);

                    % Send notification
                    notify(obj, 'JobFailed'); 
                else % This job was canceled via SLURM
                    obj.HandleState('FAILED'); % Mark this as a failure and handle appropriately
                end
                
            otherwise % Handle all other cases (included corrupted status request)
                % Do nothing
                
                
        end
        
    end 
end
    
end
