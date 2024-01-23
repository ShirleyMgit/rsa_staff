function runPileSimilarityCorBetweenMaps(nsb,areaM)
% The function takes pile similarity within map searchlight result and
% correlate it between maps
runN  = 5;

% specific calculation definition:
nFiles = 16;
nameOF = 'corPileSimBetMaps_L200_Melements';
funPointer = @pileSimilarityCorBetMaps;

% path definitions:
dirT = 'allCatchIn1PileRegPaddVspNnative040519';

% path definitions:
based = '/home/smark/fMRI_ana';
spm_path = '/data/smark/spm';
data_path = '/data/smark/fmri_sub_preproc_dir/';
firstLsearchL_path = [data_path,'PileNregNative300419All/searchLightAll/'];
addpath(spm_path)
addpath(data_path)
rsa_tool_path = 'rsatoolbox';
addpath(genpath(fullfile(data_path,rsa_tool_path)))

% name of searchlight maps results: starts from 5
SL_file_name = 'SimDistPileSpN_rz';

% directory for temporary files:
Datafile = [data_path,'tempFiles/'];
% participants number:
vnpar = [43,24,22,21,20,25,26,28:30,18,31,32,34:38,40,42];

nameDirO = [data_path,dirT,'All/searchLightAll/seconeLeverlSL/sub',num2str(vnpar(nsb)),'/',areaM,'/'];
mkdir(nameDirO);

% the names of the images created by the searchlight:
outFiles = cell(1,nFiles);
for n=1:nFiles
    outFiles{n} = [nameDirO,nameOF,num2str(n),'.nii'];
end

% loading first level searchlight images:

inFiles = [];
for n = 5:8
    file1 = cellstr(spm_select('FPList', [firstLsearchL_path,'sub',num2str(vnpar(nsb)),'/wholeBrain/',SL_file_name,num2str(n),'.nii']));
    inFiles = [inFiles;file1];
end


%loading searchlight definitions:
nameDir = [data_path,'NativeSearchLightDefinitions/',areaM,'/s',num2str(vnpar(nsb))];
load([nameDir,'L200.mat']);
runSearchlight_ShirleyRunNtempD(L,inFiles,outFiles,funPointer,nsb,runN,Datafile,'optionalParams',{});


