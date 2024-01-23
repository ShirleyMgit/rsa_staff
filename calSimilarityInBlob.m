function calSimilarityInBlob(dirT,maskName,blobName)

spm_path = '/data/smark/spm';
data_path = '/data/smark/fmri_sub_preproc_dir/';
cleaned_data_path = '/data/smark/fmri_sub_preproc_dir/fsl2spmFix_BasisSetExp';

addpath(spm_path);

% data folder:
if isempty(blobName)
    blobDataDir = fullfile(cleaned_data_path,'blobsData');
else
    blobDataDir = fullfile(cleaned_data_path,'blobsData',blobName);
end
% maskName = 'Ttmap_proBothSameMap_noncolHippo_T3';
% dirT = 'PileNregNativecleaned_MotionCSFonly_newMask';

% participants number:
vnpar = [51,50,49,48,46,45,44,43,42,40,38:-1:34,32:-1:28,26:-1:24,22:-1:18];
%subjRSA={'Ma_S51_s','Yue_S50_s','AnV_S49_s','UBh_S48_s','DeC_S46_s','WuT_S45_s','ZiL_S44_s','JoS_S43_s','HnJ_S42_s','AnL_S40_s','RaL_S38_s','YiL_S37_s','NeV_S36_s','zaK_S35_s','BeT_S34_s','JaC_S32_s','ANo_S31_s','JYe_S30_s','SVa_S29_s','Lu_S28_s','SWe_S26_s','NSi_S25_s','EBe_S24_s','zMa_S22_s','ZJL_S21_s','RSe_S20_s','SFe_S19_s','CHu_S18_s'};%,'CHu_S18_s'};%NSi_S25_s'

sb = 10;
subDir = fullfile(blobDataDir,maskName,['sub',num2str(vnpar(sb))]);

load(fullfile(subDir,'blobD.mat'));
load(fullfile(subDir,'vlenS.mat'));
% load SPM.mat:
load(fullfile(data_path,['sub',num2str(vnpar(sb))],'stats',dirT,'SPM.mat'));
% whitening
[betaAll00,resMS,Sw_hat,beta_hat,shrinkage,trRR]=noiseNormalizeBeta(DataM,SPM);

numPile = 10;
zeroBetas =  5 +2 +7;

v = ones(1,numPile*4);
bn0 =  [zeros(1,5) v zeros(1,numPile) zeros(1,zeroBetas)];%[zeros(1,5) v zeros(1,5)];%zeros in the end - ,the catch, the start and the movement regressors + ventricles regressors
bn = [bn0 bn0 bn0 bn0];
bnb = [bn 0 0 0 0];
nwb = length(bn);
beta_num = sum(bn);
b_idx0 = 1:nwb;
b_idx = b_idx0(bn==1);
betaAll1 = betaAll00(b_idx,:);

DistMs = creatDistStNew();

np = 10;
lvvoxels = length(betaAll1(1,:));
betaAll0 = reshape(betaAll1,[40,4,lvvoxels]);
betaAllLooAv = zeros(40,4,lvvoxels);

nr = 1:4;
%% LOO averages:
for r = 1:4
    betaAllLooAv(:,r,:) = mean(betaAll0(:,nr~=r,:),2);
end

ch1all = zeros(10,10);
ch2all = zeros(10,10);
cc1all = zeros(10,10);
cc2all = zeros(10,10);

for r = 1:4
    %averages:
    %hex 1
    conSHex1av = squeeze(betaAllLooAv(1:10,r,:));
    conSHex1av = conSHex1av - repmat(mean(conSHex1av),np,1);%substract the mean over condition
    % hex 2
    conSHex2av = squeeze(betaAllLooAv(11:20,r,:));
    conSHex2av = conSHex2av - repmat(mean(conSHex2av),np,1);
    % cluster 1:
    conSC1av = squeeze(betaAllLooAv(21:30,r,:));
    conSC1av = conSC1av - repmat(mean(conSC1av),np,1);
    % cluster 2:
    conSC2av = squeeze(betaAllLooAv(31:40,r,:));
    conSC2av = conSC2av - repmat(mean(conSC2av),np,1);
    
    % run patterns:
    %hex 1
    conSHex1r = betaAll1((r-1)*40+1:(r-1)*40+10,:);
    conSHex1r =  conSHex1r - repmat(mean(conSHex1r),np,1);
    %hex 2
    conSHex2r = betaAll1((r-1)*40+11:(r-1)*40+20,:);
    conSHex2r = conSHex2r - repmat(mean(conSHex2r),np,1);
    %cluster 1:
    conSC1r = betaAll1((r-1)*40+21:(r-1)*40+30,:);
    conSC1r = conSC1r - repmat(mean(conSC1r),np,1);
    %cluster 2:
    conSC2r = betaAll1((r-1)*40+31:(r-1)*40+40,:);
    conSC2r = conSC2r - repmat(mean(conSC2r),np,1);
    
    % correlation matrices:
    % hex 1:
    Ch1  = conSHex1av*conSHex1r';
    n1 = sqrt(sum(conSHex1av.*conSHex1av,2));
    n2 = sqrt(sum(conSHex1r.*conSHex1r,2));
    normM1 = n1*n2';%sqrt(diag(Ch1 )*diag(Ch1)');% not sure about this normaization..
    Ch1 = atanh(Ch1./normM1);
    ch1all = ch1all + Ch1;
    
    % hex 2:
    
    Ch2  = conSHex2av*conSHex2r';
    n1 = sqrt(sum(conSHex2av.*conSHex2av,2));
    n2 = sqrt(sum(conSHex2r.*conSHex2r,2));
    normM2 = n1*n2';%sqrt(diag(Ch2 )*diag(Ch2)');
    Ch2 = atanh(Ch2./normM2);
    ch2all = ch2all + Ch2;
    
    % cluster 1:
    CC1  = conSC1av*conSC1r';
    n1 = sqrt(sum(conSC1av.*conSC1av,2));
    n2 = sqrt(sum(conSC1r.*conSC1r,2));
    normM3 = n1*n2';%sqrt(diag(CC1)*diag(CC1)');
    CC1 = atanh(CC1./normM3);
    cc1all = cc1all + CC1;
    
     % cluster 2:
    CC2  = conSC2av*conSC2r';
    n1 = sqrt(sum(conSC2av.*conSC2av,2));
    n2 = sqrt(sum(conSC2r.*conSC2r,2));
    normM4 = n1*n2';%sqrt(diag(CC2)*diag(CC2)');
    CC2 = atanh(CC2./normM4);
    cc2all = cc2all + CC2;
end

ch1all = tanh(ch1all / 4);
ch2all = tanh(ch2all / 4);
cc1all = tanh(cc1all / 4);
cc2all = tanh(cc2all / 4);

[r1,p1] = corr(ch1all(:),DistMs.Hex(:),'type','Spearman');
[r2,p2] = corr(ch2all(:),DistMs.Hex(:),'type','Spearman');
[cr1,cp1] = corr(cc1all(:),DistMs.ClustBig(:),'type','Spearman');
[cr2,cp2] = corr(cc2all(:),DistMs.ClustSmall(:),'type','Spearman');
rVec = [r1,r2,cr1,cr2];
pVec = [p1,p2,cp1,cp2]; 
SimResult = struct('rVec',rVec,'pVec',pVec,'lvoxels',lvvoxels,'ch1all',ch1all,'ch2all',ch2all,'cc1all',cc1all,'cc2all',cc2all);
save(fullfile(blobDataDir,maskName,'SimResult.mat'),SimResult);