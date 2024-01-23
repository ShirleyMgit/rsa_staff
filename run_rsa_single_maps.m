function run_rsa_single_maps(nsb,dirT,nFiles,nameSearchL,dirL)
% experiment details: numPile= 10; nMap = 4; 
nSess = 4;

runN  = 'single_map_rsa';
% path definitions:
rsa_tool_path = 'rsatoolbox';

% enter relevant path:
root = '/home/fs0/abaram/scratch/shirley/alon';
addpath(genpath(fullfile(root,'code')));
addpath(genpath('/home/fs0/abaram/scratch/MATLAB/spm12'))

data_path = '/data/holly-host/smark/fmri_sub_preproc_dir/';
cleaned_data_path = '/data/holly-host/smark/fmri_sub_preproc_dir/fsl2spmFix_BasisSetExp/sub-';
addpath(genpath(fullfile(data_path,rsa_tool_path)))

addpath(spm_path)
addpath(data_path)

% directory for temporary files:
Datafile = [data_path,'tempFiles/'];

subjects = cell(1,28);
% load SPM.mat:
load(fullfile('path to spm', 'SPM.mat'));

areaM = 'wholeBrainMcorrect';
nameDirO = [data_path,dirT,'/searchLightAllEvTim_replicate/sub',num2str(vnpar(nsb)),'/',areaM,'/'];
mkdir(nameDirO);
% the names of the images created by the searchlight:
outFiles = cell(1,nFiles);
for f=1:nFiles
    outFiles{f} = [nameDirO,nameSearchL,num2str(f),'.nii'];
end


% loading images:
inFiles = [];
for s = 1:nSess
    allfiles = cellstr(spm_select('FPList', [cleaned_data_path,num2str(vnpar(nsb)),'/sess0',num2str(s),'/'], '^r_cleaned_smoothed_.*.nii$'));%not warp!
    inFiles = [inFiles; allfiles];
end
% loading searchlight definitions:
nameDir = [data_path,'NativeSearchLightDefinitions/',areaM,'/',dirL,'/s',num2str(vnpar(nsb))];
load([nameDir,'/',nameSearchL,'.mat']);
runSearchlight_ShirleyRunNtempDdifName(L,inFiles,outFiles,@TimNormLoo_RSAtoolboxMatrixEoutNanInfAllEv,nsb,runN,Datafile,'optionalParams',{SPM});






