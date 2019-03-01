MERFISH_analysis installation
Leonardo Sepulveda, 01/18/19

Below is the steps I used to install the pipeline in the Harvard cluster. It assumes that you have already installed git, anaconda3, matlab, fftw and scons in your cluster (They are called using the module load command). Besides that, I noticed that some of the files in the package need small modifications, so I copied them to the folder modifiedFunctions. You will need to replace the version that comes from git for the ones provided. 

I provide a master script called runMERFISH that initializes variables and calls the main MERFISH script (MERFISHScheduler). I also provide a dataOrganization.csv file for a multi-z experiment. 

I tested the pipeline and I get the same results I got with the old pipeline. 

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
