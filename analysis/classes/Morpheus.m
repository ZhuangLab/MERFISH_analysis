classdef Morpheus < handle
% ------------------------------------------------------------------------
% obj = Morpheus(varargin)
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% lmoffitt@mcb.harvard.edu
% September 21, 2017
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2018.
% -------------------------------------------------------------------------
% This class handles sending email messages from various SLURM jobs
% To configure matlab to send these emails, one must run the following
% commands at least once (where the appropriate values are provided)
% setpref('Internet', 'E_mail', email_account);
% setpref('Internet', 'SMTP_Server',smpt_server);
% setpref('Internet', 'SMTP_Username',email_username);
% setpref('Internet', 'SMTP_Password',email_password);


% -------------------------------------------------------------------------
% Define properties
% -------------------------------------------------------------------------
properties
   verbose = true       % Control the verbosity of the class
   name = ''            % The name of this morpheus instance
end

properties (SetAccess=protected)
    
    % Listeners 
    successListeners    % Listeners for success of a class
    errorListeners      % Listeners for failure/error of a class    
    
    % Properties associated with email recipients
    recipient = ''      % Email recipient
    
    % Properties associated with errors
    maxNumErrors = Inf  % The maximum number of errors for which messages will be sent
    numErrors = 0       % The number of accumulated error signals
    
    % Sent messages archive
    sentMessages = cell(0,5)   % An archive of the messages sent
end


% -------------------------------------------------------------------------
% Public Methods
% -------------------------------------------------------------------------
methods
    
    % -------------------------------------------------------------------------
    % Define constructor
    % -------------------------------------------------------------------------
    function obj = Morpheus(recipient, varargin)
        % Create a messenger class
        
        % -------------------------------------------------------------------------
        % Parse variable inputs
        % -------------------------------------------------------------------------
        % Define defaults
        defaults = cell(0,3); 
        
        % Properties for class behavior
        defaults(end+1,:) = {'verbose', ...                  
            'boolean', true};
        defaults(end+1,:) = {'name', ...
            'string', ''};
        defaults(end+1,:) = {'maxNumErrors', ...
            'nonnegative', Inf};
        
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
        
        % Define recipient
        obj.recipient = recipient;
        
        % Start email server
        % Gmail server.
        props = java.lang.System.getProperties;
        props.setProperty('mail.smtp.auth', 'true');
        props.setProperty('mail.smtp.starttls.enable', 'true');  
        props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory' );
        props.setProperty('mail.smtp.socketFactory.port', '465');

        
    end
    
    % -------------------------------------------------------------------------
    % Send message
    % -------------------------------------------------------------------------
    function SendMessage(obj, subject, body, sourceObj)
        
        % Handle no source object provided
        if ~exist('sourceObj', 'var')
            sourceObj = [];
        end
        
        % Handle multi-line bodies
        body = cellstr(body);
        
        % Append the morpheus name to the body
        if ~isempty(obj.name)
            newBody = cell(1, length(body)+1);
            newBody(2:end) = body(:);
            newBody{1} = ['Message from morpheus: ' obj.name];
            body = newBody;
        end

        % Format body for sendmail
        body = sprintf('%s\n', body{:});
        
        % Send a message if the recipient has been defined
        if ~isempty(obj.recipient)
            % Send the message
            sendmail(obj.recipient, subject, body);
            
            % Archive the message
            obj.sentMessages(end+1,:) = {obj.recipient, subject, body, datetime('now'), sourceObj};
            
        end
        
    end
    
    % -------------------------------------------------------------------------
    % Add success listener
    % -------------------------------------------------------------------------
    function AddSuccessListener(obj, sourceObj, eventName)
        
        if isempty(obj.successListeners)
            obj.successListeners = addlistener(sourceObj, eventName, @(source,event)obj.HandleSuccessMessageRequest(source, event));
        else
            obj.successListeners(end+1) = addlistener(sourceObj, eventName, @(source,event)obj.HandleSuccessMessageRequest(source, event));
        end
    end

    % -------------------------------------------------------------------------
    % Add error listener
    % -------------------------------------------------------------------------
    function AddErrorListener(obj, sourceObj, eventName)
        
        if isempty(obj.errorListeners)
            obj.errorListeners = addlistener(sourceObj, eventName, @(source,event)obj.HandleErrorMessageRequest(source, event));
        else
            obj.errorListeners(end+1) = addlistener(sourceObj, eventName, @(source,event)obj.HandleErrorMessageRequest(source, event));
        end
    end
    
    % -------------------------------------------------------------------------
    % Handle success message request
    % -------------------------------------------------------------------------
    function HandleSuccessMessageRequest(obj, sourceObj, event)
        % Compose message from source object
        
        [header, body] = sourceObj.Status();
                        
        % Send mail
        obj.SendMessage(header, body, sourceObj);
        
        if obj.verbose
            disp(['Sending message from ' sourceObj.name ' an instance of ' class(sourceObj)]);
        end
    end

    % -------------------------------------------------------------------------
    % Handle error message request
    % -------------------------------------------------------------------------
    function HandleErrorMessageRequest(obj, sourceObj, event)
        
        % Compose message from source object
        [header, body] = sourceObj.Status();
        
        % Check accumulated number of error messages
        obj.numErrors = obj.numErrors + 1;
        
        % Handle case that the number of accumulated errors is less than
        % the maximum
        if obj.numErrors < obj.maxNumErrors
        
            % Send mail
            obj.SendMessage(header, body, sourceObj);

            if obj.verbose
                disp(['Sending message from ' sourceObj.name ' an instance of ' class(sourceObj)]);
            end
        else % Handle the case that it is not
            if obj.numErrors == obj.maxNumErrors % Send a message when the maximum is reached
                header = 'Morpheus Maximum Number of Errors';
                
                % Create custom body
                body = cell(0,1);
                body{1} = ['Morpheus has received the maximum number of error signals: ' num2str(obj.maxNumErrors) '.'];
                body{2} = ['No further error messages will be sent.'];
                body = char(body);
                
                obj.SendMessage(header, body);
            end
        end
    end
    
    

end
    
end
