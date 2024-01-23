function makeBlobRSA(sb)

based = '/home/smark';%/fMRI_ana';
spm_path = '/data/smark/spm';
data_path = '/data/smark/fmri_sub_preproc_dir/';
cleaned_data_path = '/data/smark/fmri_sub_preproc_dir/fsl2spmFix_BasisSetExp/sub-';
addpath(spm_path)

clusterThreshold = 2.0;

dirT = 'PileNregNativecleaned_MotionCSFonly_newMaskAllvarcleanedF';
searchL = 'L100_FSLcleaned280920b';
dirR = 's6L100_PileNregNativecleaned_MotionCSFonly_newMaskvarcleanedFnSub28b';

statFolder = [data_path,dirT,'\searchLightResults\wholeBrainMcorrect\',searchL,'\',dirR,'\'];
maskName = 'ECblob_proHexsameMap_noncol.nii';

% number of voxels in the mask:
con_tmp = spm_vol(maskName); % SPM function to get image info
[ROI_dat1,XYZ1] = spm_read_vols(con_tmp,0);
vmask = ROI_dat1(:);
nmaskV = sum(vmask>0);

% save experiment length:
vlSes = zeros(4,1);
inFiles = [];
for s = 1:nSess
    allfiles = cellstr(spm_select('FPList', [cleaned_data_path,num2str(vnpar(sb)),'/sess0',num2str(s),'/'], '^r_cleaned_smoothed_.*.nii$'));%not warp!
    lenS = length(allfiles);
    vlSes(s) = lenS;
    %inFiles = [inFiles; allfiles];
end

save([cleaned_data_path,num2str(vnpar(sb)),'vlenS.mat'],'vlSes');
DataM = zeros(sum(vlSes),nmaskV);
% load data and mask it:
for s = 1:nSess
    allfiles = cellstr(spm_select('FPList', [cleaned_data_path,num2str(vnpar(sb)),'/sess0',num2str(s),'/'], '^r_cleaned_smoothed_.*.nii$'));%not warp!
    lenS = length(allfiles);
    for f = 1:lenS
        fmaskN = [allfiles{f},'_mainB.nii'];
        maskImMine(maskName, allfiles{f},fmaskN ,0.5);
        con_tmp = spm_vol(fmaskN); % SPM function to get image info
        [ROI_dat1,XYZ1] = spm_read_vols(con_tmp,0);
        D1 = ROI_dat1(:);
        Dblob = D1(D1>0);
        if s>1
            DataM (vlSes(s-1)+f,:) = Dblob;
        else
            DataM (f,:) = Dblob;
        end
    end
end



