%% Startup Script
% -------------------------------------------------------------------------
% This script initializes the matlab workspace and defines useful paths and
% global variables for STORM analysis and MERFISH analysis
% If you already have a startup script add the following code to this
% script. To be functional, there are specific paths that must be set based
% on the local machine.  
% -------------------------------------------------------------------------
%% Clear Existing Workspace
close all;
clear all;
clc;
display('------------------------------------------------------------------');
warning off all
restoredefaultpath; % Clear previous paths
warning on all

%% Define global Variables and Paths
global scratchPath;       % place for matlab-storm to save temporary files
global pythonPath;        % path to Python 2.7
global matlabStormPath;   % path to matlab-storm
global stormAnalysisPath; % path to storm-analysis

homePath          = '/n/home06/lsepulvedaduran';
basePath          = [homePath '/Software/new_merfish_pipeline/'];
scratchPath       = [homePath '/scratch/']; 
pythonPath        = [homePath '/.conda/envs/merfish_analysis/lib/python2.7/']; 
matlabStormPath   = [basePath 'matlab-storm/'];  
stormAnalysisPath = [basePath 'storm-analysis/'];  


%% Call the matlab-storm startup script
display('Adding matlab-storm');
addpath([matlabStormPath,'Startup/']);
matlabstorm_startup;    % Configure matlab-storm

%% Enable execution of daoSTORMExe
display('Adding path to daoSTORM');
daoSTORMexe = [homePath  '/.conda/envs/merfish_analysis/bin/python ' stormAnalysisPath '3d_daostorm/mufit_analysis.py'];
display(['    daoSTORMexe = ' daoSTORMexe]);
	
%% Add merfish_analysis: downloaded from https://github.com/ZhuangLab/MERFISH_analysis
display('------------------------------------------------------------------');
display('Adding MERFISH_analysis');
MERFISHAnalysisPath = [basePath 'MERFISH_analysis/']; % Path the folder in which you installed this software package
addpath(genpath(MERFISHAnalysisPath), '-begin');
display(['    ' MERFISHAnalysisPath]);
display(['      And all enclosed paths']);

%% Add export_fig
display('------------------------------------------------------------------');
display('Adding export_fig');
addpath([ basePath 'export_fig/']); % Path to export_fig
