%--------------------------------------------------------------------------
% MERFISH_analysis installation
%--------------------------------------------------------------------------
% Leonardo Sepulveda, 03/01/2019
%
% Below are the general steps used to install the pipeline in the Harvard 
% RC  cluster. This installation requieres to have the following packages 
% already installed:
% - git
% - anaconda3 (for python 2.7) 
% - matlab
% - fftw
% - scons
%
% In the installation below, these packages are called with the "module 
% load" command. 
% 
% The runMERFISH.m script is the main script used to launch the pipeline. 
% It contains the paths to several files required for the analysis, and
%  also calls the MERFISHscheduler.m script, that is the main script to 
% run the decoding. The MERFISH_analysis/example_scripts/decoding folder 
% also contains 3 files required to seamlessly run the pipeline: 
% - L26E1_codebook.csv: file containing the mapping between gene names 
%   and their barcodes for the L26 library designed to show the 
%   performance of MERFISH in U-2 OS cells (Chen et al., Science (2015), 
%   Moffit et al.,PNAS(2016a); Moffit et al.,PNAS(2016b)).  
% - FPKMDataPublished.matb: File containing the expression level the 
%   genes targeted in the L26 library, measured by RNAseq (average of 3 
%   replicates, Moffit et al., PNAS (2016a))
% - data_organization.csv: file containing information on the 
%   organization of the data files, for example, the names of the files, 
%   the frames in a multi-Z image associated to each bit for each round 
%   of hybridization, the frames containing the images for segmentation 
%   and image registration, etc.  
%--------------------------------------------------------------------------
INSTALLATION INSTRUCTIONS:
%--------------------------------------------------------------------------

% Create folders
mkdir new_merfish_pipeline
cd new_merfish_pipeline

% download repositories from github 
git clone https://github.com/ZhuangLab/matlab-storm.git
git clone https://github.com/ZhuangLab/MERFISH_analysis.git
git clone https://github.com/altmany/export_fig.git
git clone https://github.com/ZhuangLab/storm-analysis.git

% return storm-analysis to an earlier state
cd storm-analysis
git reset --hard 73866042aad3e3ceb1f5f4a766f0ca9ef54eb890

% generate python environment
module load Anaconda3/5.0.1-fasrc01
conda create -n merfish_analysis python=2.7
source activate merfish_analysis

% Installing dependencies
pip install matplotlib
conda install tifffile -c conda-forge
pip install h5py
pip install pillow
pip install scipy
conda update --all

% install storm analysis
module load scons
module load fftw
cd storm-analysis
chmod +x compile_all_linux.sh
./compile_all_linux.sh

% test storm-analysis
cd tests
python run_tests.py

% Add a new alias to the .bashrc file to easy loading of merfish environment
# load merfish related modules, test implementation
alias loadMF='module load centos6/0.0.1-fasrc01; module load Anaconda3/5.0.1-fasrc01; source activate merfish_analysis; module load fftw; module load matlab/R2017a-fasrc02; export MATLAB_USE_USEWORK=1'

% create folder for MATLAB startup, 
cd ~
mkdir Documents
cd Documents
mkdir MATLAB % folder for MATLAB startup
cp startup.m ~/Documents/MATLAB/

% create folder for scratch files
cd ~
mkdir scratch % scratch folder for matlab-storm

% Move to a cluster node and run the merfish pipeline
srun -p zhuang --mem 8000 -t 0-48:00 --pty /bin/bash
loadMF
matlab
>> runMERFISH
