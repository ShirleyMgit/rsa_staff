function runRSAinBlob(sb,dirT,clusterThreshold)

% participants number:
vnpar = [51,50,49,48,46,45,44,43,42,40,38:-1:34,32:-1:28,26:-1:24,22:-1:18];

rsa_tool_path = 'rsatoolbox';
based = '/home/smark';%/fMRI_ana';
spm_path = '/data/smark/spm';
data_path = '/data/smark/fmri_sub_preproc_dir/';
cleaned_data_path = '/data/smark/fmri_sub_preproc_dir/fsl2spmFix_BasisSetExp/';
addpath(spm_path)
addpath(genpath(fullfile(data_path,rsa_tool_path)))

blobDataDir = fullfile(cleaned_data_path,'blobsData');

maskName = ['proHexsameMap_noncol_maskT',num2str(clusterThreshold),'overlapB'];

newDir = fullfile(blobDataDir,[maskName,'_maskBrain']);
subDir = fullfile(newDir,['sub',num2str(vnpar(sb))]);

load(fullfile(subDir,'blobD.mat'));

% load SPM.mat:
load([data_path,'sub',num2str(vnpar(sb)),'/stats/',dirT,'/SPM.mat']);

H = TimNormLoo_RSAtoolboxMatrixEoutNanInfExpVarTryNoiseNormP(DataM,SPM);

save(fullfile(subDir,'TimProj.mat'),'H');