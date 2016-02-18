function BatchLaunchOligoArray(oligoArrayCommandtemp,geneFasta,varargin)
%--------------------------------------------------------------------------
% Output
%--------------------------------------------------------------------------
% BatchLaunchOligoArray(oligoArrayCommandtemp,geneFasta,varargin)
% Outputs GENENAME_oligos.txt files for each gene in geneFasta into the
% specified save folder.  This text files contains probes that target the
% gene GENENAME and don't cross-hybridize to any other probes in the
% BLASTdatabase. The probes also meet the specified TM, GC, and secondary
% structure constraints (see Optional Input). 
% 
%--------------------------------------------------------------------------
% Required Input
%--------------------------------------------------------------------------
% oligoArrayCommandtemp - a specialized string produced by OligoArrayCmd 
%              see the OligoArrayCmd function
% geneFasta - a fasta file containing all the genes you want to design
%              probes for. 
% 
%--------------------------------------------------------------------------
% Optional Input
%--------------------------------------------------------------------------
% {'savePath', 'string', ''};
% {'headerSeps', 'string',''};
% {'maxFragment', 'positive',1E3};
% {'batchsize', 'positive',4 };
% {'maxTime', 'positive',30 };%  in minutes
% {'maxBlastTime', 'positive', 15};%  in minutes
% {'maxCPU', 'positive', 75};  %   Percent of total CPU at which to halt launching jobs. 
% {'refreshTime', 'positive',20 }; % time to wait inbetween checking for CPU (in seconds)
% {'runExternal', 'boolean',true };
% {'verbose', 'boolean',true };
% {'runSbatch', 'boolean',false }; % For harvard RC server only 
% {'remotePath', 'string','/n/home05/boettiger/'}; % For harvard RC server only
% 
%--------------------------------------------------------------------------
% Alistair Boettiger
% boettiger.alistair@gmail.com
% June 07 2014
%--------------------------------------------------------------------------
% Copyright Presidents and Fellows of Harvard College, 2016.


% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'savePath', 'string', ''};
defaults(end+1,:) = {'headerSeps', 'string',''};
defaults(end+1,:) = {'maxFragment', 'positive',1E3};
defaults(end+1,:) = {'batchsize', 'positive',4 };
defaults(end+1,:) = {'maxTime', 'positive',30 };%  in minutes
defaults(end+1,:) = {'maxBlastTime', 'positive', 15};%  in minutes
defaults(end+1,:) = {'maxCPU', 'positive', 75};  %   Percent of total CPU at which to halt launching jobs. 
defaults(end+1,:) = {'refreshTime', 'positive',20 }; % time to wait inbetween checking for CPU (in seconds)
defaults(end+1,:) = {'runExternal', 'boolean',true };
defaults(end+1,:) = {'verbose', 'boolean',true };
defaults(end+1,:) = {'runSbatch', 'boolean',false }; % For harvard RC server only 
defaults(end+1,:) = {'remotePath', 'string','/n/home05/boettiger/'}; % For harvard RC server only

if nargin < 2
    error('matlabSTORM:invalidArguments', 'requires an oligoArrayCommand string and geneGasta');
end

parameters = ParseVariableArguments(varargin, defaults, mfilename);

%% parsing input a bit more
if isempty(parameters.savePath) 
    s = strfind(oligoArrayCommandtemp,'-o ');
    e = strfind(oligoArrayCommandtemp,'-r ');
    parameters.savePath = extractpath(oligoArrayCommandtemp(s+3:e-3)); 
end

if parameters.runSbatch % if running remotely, can't use windows system batch operations
    parameters.runExternal = false;
end

cmdLogFile = [parameters.remotePath,'log.txt'];
cmdErrorFile = [parameters.remotePath,'err.txt'];
remoteLogFile = [parameters.remotePath,'sbatchLog.txt'];
remoteErrorFile = [parameters.remotePath,'sbatchErr.txt'];



%% Launch jobs

% Go back and loop over 
Gs = length(geneFasta);
prc = cell(Gs,1);
for g=1:Gs
    locusName = geneFasta(g).Header;
    sequence = geneFasta(g).Sequence;
    
    % Trim gene_name and remove illegal characters
    locusName = regexprep(locusName,...
        {'\.','\:','\s',',','\\','\/'},{'P','C','_','',']','['}); 
    if ~isempty(parameters.headerSeps)
        sps = strfind(geneFasta(g).Header,parameters.headerSeps);
        locusName = locusName(sps(1)+length(parameters.headerSeps):sps(2)-1);
    end
        
    % Decompose region into smaller fragments
    seq_length = length(sequence); 
    locusSaveName =  locusName;
    gene_fragment_seqs = cell(100,1);
    gene_fragment_names = cell(100,1); 
   subfragment_start = 1;
   subfragment_end = min(seq_length,parameters.maxFragment);
   f = 0; % fragment counter
   while  subfragment_start < seq_length 
        f = f+1;
        gene_fragment_seqs{f} = ...
        sequence(subfragment_start:subfragment_end);
        subfragment_start = subfragment_start + parameters.maxFragment;
        subfragment_end = min(subfragment_end + parameters.maxFragment,seq_length);
        gene_fragment_names{f} = [locusName,'_pt',num2str(f)];
   end

    hasdata = logical(1-cellfun(@isempty,gene_fragment_seqs));
    gene_fragment_names = gene_fragment_names(hasdata);
    gene_fragment_seqs = gene_fragment_seqs(hasdata);
    Gene_fasta.Sequence = gene_fragment_seqs;
    Gene_fasta.Header = gene_fragment_names;

    fileout = [parameters.savePath,locusSaveName,'.fasta'];
    WriteFasta(fileout,Gene_fasta,[],'Append',false,'Warnings',false); 
    

    disp(locusSaveName)
    % Replace 'fastainput' with the fasta name for this gene.
    % Then send OligoArrayCommand to Queue.  
    oligoArrayCommand = regexprep(oligoArrayCommandtemp,'genename',locusSaveName);

    if parameters.verbose
 display('-----------------------------------------------------------------');
 display(['Running file ',num2str(g),' of ',num2str(Gs),':'])
 display(['     ' oligoArrayCommand]);
 display('-----------------------------------------------------------------');
    end

    if ~isempty(parameters.maxCPU) && ~parameters.runSbatch
        waitforfreecpu('MaxLoad',parameters.maxCPU,...
                     'RefreshTime',parameters.refreshTime,...
                     'verbose',parameters.verbose);
    end
    
    if parameters.runExternal
       prc{g} = SystemRun([oligoArrayCommand,' &'],'Hidden',true);
       Nrunning = inf;
       t=0;
       while Nrunning >= parameters.batchsize 
           prcLaunched = prc(~cellfun(@isempty, prc));
           hasexited = cellfun(@(x) x.HasExited,prcLaunched);
           Nrunning = length(prcLaunched)-sum(hasexited);

             % Measure run times for all processes
            % Kill processes that have been running too long. 
            startTimes = cell2mat(cellfun(@(p) double(p.StartTime.Day*60*24 +...
                p.StartTime.Hour*60 + p.StartTime.Minute),...
                prcLaunched,'ErrorHandler',@cellfunempty,'UniformOutput',false));
            startTimes(hasexited) = inf; 
            currTime = double(System.DateTime.Now.Day*60*24 + ...
                System.DateTime.Now.Hour*60 + System.DateTime.Now.Minute);
            ranTooLong = find(currTime - startTimes > parameters.maxTime);
            for k = 1:length(ranTooLong)
                pID = prcLaunched{ranTooLong(k)}.Id;
                stopProcess = ['wmic process ',num2str(pID),' delete '];
                disp(['Analysis exceeded ',num2str(parameters.maxTime), 'min. ' 'Aborting process ',pID]); 
                dos(stopProcess);
            end 

             pause(.05); 
            t = t + 1;
           if t/.05 >= parameters.maxBlastTime
               ProcessTimeout('blastall.exe','maxTime',parameters.maxBlastTime*60);
           end
       end 
     
    elseif parameters.runSbatch
        batFile = [parameters.remotePath,'Run',locusSaveName,'.bat'];
        fid = fopen(batFile,'w+');       
        fprintf(fid,'#!/bin/bash','');          fprintf(fid,'%s\n',''); %#ok<CTPCT>
        fprintf(fid,['#SBATCH -n ','1'],'');    fprintf(fid,'%s\n','');
        fprintf(fid,['#SBATCH -t ',num2str(2*parameters.maxTime)],''); fprintf(fid,'%s\n','');
        fprintf(fid,['#SBATCH -p ','serial_requeue'],''); fprintf(fid,'%s\n','');
        fprintf(fid,['#SBATCH --mem','=3000'],''); fprintf(fid,'%s\n','');
        fprintf(fid,['#SBATCH -o ',remoteLogFile],'');         fprintf(fid,'%s\n','');
        fprintf(fid,['#SBATCH -e ',remoteErrorFile],'');       fprintf(fid,'%s\n','');
        fprintf(fid,[oligoArrayCommand,' > ',cmdLogFile,' 2> ',cmdErrorFile]);
        fclose(fid);
        system(['sbatch ',batFile]); 
        
        MaxRemoteJobs(parameters.batchsize,'currTask',g,'totTasks',Gs,'verbose',parameters.verbose); 
        CleanupDaemon(pwd); 
        
    else
        system(oligoArrayCommand); 
    end
    
end

if parameters.runExternal
    % Wait till final processes finish
    Nrunning = length(prcLaunched)-sum(hasexited);
    while Nrunning > 0
       prcLaunched = prc(~cellfun(@isempty, prc));
       hasexited = cellfun(@(x) x.HasExited,prcLaunched);
       Nrunning = length(prcLaunched)-sum(hasexited);
       pause(5);
       ProcessTimeout('blastall.exe','maxTime',parameters.maxBlastTime*60);
    end
end