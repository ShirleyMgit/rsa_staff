function groupLevel_rsa(root,subjects,slName)
% get a tstat for group level using a t-test (for a quick parametric test -
% stats for the paper will be performed non-parametrically with PALM)
rsa_toolbox_path = fullfile('from_git', 'util','rsaToolbox');
addpath(rsa_toolbox_path);

nFiles = 16;
smoothKernel = 6;
matFname = [slName '_rsa???'];%% I DONT KNOW HOW YOU HAVE CALLED THE FILES

niiDims = [91,109,91]; % dims of standard brain
projMat = zeros(niiDims(1),niiDims(2),niiDims(3),nFiles,length(subjects)); 
for iSub=1:length(subjects)
    sub = subjects{iSub};
    subspaceGener_dir = fullfile(root,'subspaceGener',sub);
    projMatFiles = cell(nFiles,1);
    for f=1:nFiles
        projMatFiles{f} = [subspaceGener_dir,'/smth',num2str(smoothKernel) 'mni_' matFname num2str(f),'.nii'];
        V           = spm_vol(projMatFiles{f}); % heatder info      
        [projMat_currElement,~] = spm_read_vols(V,0);    
        projMat(:,:,:,f,iSub) = projMat_currElement;
    end
end

%% diffrences between z-score distances:
% Hex:
Hex1_1to2 = squeeze(projMat(:,:,:,1,:) - projMat(:,:,:,2,:));
Hex2_1to2 = squeeze(projMat(:,:,:,4,:) - projMat(:,:,:,5,:));
Hex1_1to3 = squeeze(projMat(:,:,:,1,:) - projMat(:,:,:,3,:));
Hex2_1to3 = squeeze(projMat(:,:,:,4,:) - projMat(:,:,:,6,:));
all_hex_1to2 = Hex1_1to2 + Hex2_1to2;
all_hex_1to3 = Hex1_1to3 + Hex2_1to3;

% Cluster:
Cl1_1to2 = squeeze(projMat(:,:,:,7,:) - projMat(:,:,:,8,:));
Cl2_1to2 = squeeze(projMat(:,:,:,10,:) - projMat(:,:,:,11,:));
Cl1_1to3 = squeeze(projMat(:,:,:,7,:) - projMat(:,:,:,9,:));
Cl2_1to3 = squeeze(projMat(:,:,:,10,:) - projMat(:,:,:,12,:));
all_Cl_1to2 = Cl1_1to2 + Cl2_1to2;
all_Cl_1to3 = Cl1_1to3 + Cl2_1to3;

%% kandel tau (b) - but it seems the same thing as caluclated by the toolbox:
Hex1_kandel = squeeze(projMat(:,:,:,13,:));
Hex2_kandel  = squeeze(projMat(:,:,:,14,:));
Cl1_kandel  = squeeze(projMat(:,:,:,15,:));
Cl2_kandel  = squeeze(projMat(:,:,:,16,:));
Hex1_2_kandel  = 0.5*(Hex1 + Hex2);
Cl1_2_kandel  = 0.5*(Cl1 + Cl2);
all_kandel = 0.5*(Hex1_2_kandel + Cl1_2_kandel);

%% stats
nVox = prod(niiDims);
Hex1_2_kandel_linVox = reshape(Hex1_2_kandel,[nVox,length(subjects)]);
Hex1_kandel_linVox = reshape(Hex1_kandel,[nVox,length(subjects)]);
Hex2_kandel_linVox = reshape(Hex2_kandel,[nVox,length(subjects)]);
All_kandel_linVox = reshape(all_kandel,[nVox,length(subjects)]);
all_hex_1to2_contrast_linVox = reshape(all_hex_1to2,[nVox,length(subjects)]);
all_hex_1to3_contrast_linVox = reshape(all_hex_1to3,[nVox,length(subjects)]);
tStats_linVox = zeros(nVox,1);
rank_p_stats = zeros(nVox,1);
rank_p_stats_All = zeros(nVox,1);
rank_p_statsHex1 = zeros(nVox,1);
rank_p_statsHex2 = zeros(nVox,1);
tStats_all_hex_1to2_linVox = zeros(nVox,1);
tStats_all_hex_1to3_linVox = zeros(nVox,1);
for iVox = 1:nVox
    [~,~,~,stats] = ttest(Hex1_2_kandel_linVox(iVox,:),0,'tail','right');
    tStats_linVox(iVox) = stats.tstat;
    rank_p_stats(iVox) = rsa.util.signrank_onesided(Hex1_2_kandel_linVox(iVox,:));
    [~,~,~,stats_hex1to2] = ttest(all_hex_1to2_contrast_linVox(iVox,:),0,'tail','right');
    tStats_all_hex_1to2_linVox(iVox) = stats_hex1to2.tstat;
    [~,~,~,stats_hex1to3] = ttest(all_hex_1to3_contrast_linVox(iVox,:),0,'tail','right');
    tStats_all_hex_1to3_linVox(iVox) = stats_hex1to3.tstat;
    rank_p_statsHex1(iVox) = rsa.util.signrank_onesided(Hex1_kandel_linVox(iVox,:));
    rank_p_statsHex2(iVox) = rsa.util.signrank_onesided(Hex2_kandel_linVox(iVox,:));
    rank_p_stats_All(iVox)  = rsa.util.signrank_onesided(All_kandel_linVox(iVox,:));
end
tStats = reshape(tStats_linVox,niiDims);
p_rank_Stats = reshape(rank_p_stats,niiDims);
p_rank_StatsHex1 = reshape(rank_p_statsHex1,niiDims);
p_rank_StatsHex2 = reshape(rank_p_statsHex2,niiDims);
p_rank_Stats_All = reshape(rank_p_stats_All,niiDims);
tStats_all_hex_1to2_Vox = reshape(tStats_all_hex_1to2_linVox,niiDims);
tStats_all_hex_1to3_Vox = reshape(tStats_all_hex_1to3_linVox,niiDims);

save4Dnii(tStats,V,fullfile(root,'subspaceGener','groupStats','tStat_hex1and2_kandel.nii'))
save4Dnii(p_rank_Stats,V,fullfile(root,'subspaceGener','groupStats','hex1and2_p_rank_Stats.nii'))
save4Dnii(p_rank_StatsHex1,V,fullfile(root,'subspaceGener','groupStats','hex1_p_rank_Stats.nii'))
save4Dnii(p_rank_StatsHex2,V,fullfile(root,'subspaceGener','groupStats','hex2_p_rank_Stats.nii'))
save4Dnii(p_rank_Stats_All,V,fullfile(root,'subspaceGener','groupStats','all_p_rank_Stats.nii'))
save4Dnii(tStats_all_hex_1to2_Vox,V,fullfile(root,'subspaceGener','groupStats','tStat_all_hex_1to2.nii'))
save4Dnii(tStats_all_hex_1to3_Vox,V,fullfile(root,'subspaceGener','groupStats','tStat_all_hex_1to3.nii'))

Cl1_2_kandel_contrast_linVox = reshape(Cl1_2_kandel,[nVox,length(subjects)]);
all_Cl_1to2_contrast_linVox = reshape(all_Cl_1to2,[nVox,length(subjects)]);
all_Cl_1to3_contrast_linVox = reshape(all_Cl_1to3,[nVox,length(subjects)]);
Cl1_kandel_linVox = reshape(Cl1_kandel,[nVox,length(subjects)]);
Cl2_kandel_linVox = reshape(Cl2_kandel,[nVox,length(subjects)]);
rank_p_statsCl1 = zeros(nVox,1);
rank_p_statsCl2 = zeros(nVox,1);
tStats_linVox = zeros(nVox,1);
rank_p_stats = zeros(nVox,1);
tStats_all_Cl_1to2_linVox = zeros(nVox,1);
tStats_all_Cl_1to3_linVox = zeros(nVox,1);
for iVox = 1:nVox

    [~,~,~,stats] = ttest(Cl1_2_kandel(iVox,:),0,'tail','right');
    tStats_linVox(iVox) = stats.tstat;
    rank_p_stats(iVox) = rsa.util.signrank_onesided(Cl1_2_kandel_contrast_linVox);

    [~,~,~,stats_Cl1to2] = ttest(all_Cl_1to2_contrast_linVox(iVox,:),0,'tail','right');
    tStats_all_Cl_1to2_linVox(iVox) = stats_Cl1to2.tstat;
    [~,~,~,stats_Cl1to3] = ttest(all_Cl_1to3_contrast_linVox(iVox,:),0,'tail','right');
    tStats_all_hex_1to3_linVox(iVox) = stats_Cl1to3.tstat;

    rank_p_statsCl1(iVox) = rsa.util.signrank_onesided(Cl1_kandel_linVox(iVox,:));
    rank_p_statsCl2(iVox) = rsa.util.signrank_onesided(Cl2_kandel_linVox(iVox,:));

end
tStats = reshape(tStats_linVox,niiDims);
p_rank_Stats = reshape(rank_p_stats,niiDims);
p_rank_StatsCl1 = reshape(rank_p_statsCl1,niiDims);
p_rank_StatsCl2 = reshape(rank_p_statsCl2,niiDims);
tStats_all_Cl_1to2_Vox = reshape(tStats_all_Cl_1to2_linVox,niiDims);
tStats_all_Cl_1to3_Vox = reshape(tStats_all_Cl_1to3_linVox,niiDims);

save4Dnii(tStats,V,fullfile(root,'subspaceGener','groupStats','tStat_Cl1and2_kandel.nii'))
save4Dnii(p_rank_Stats,V,fullfile(root,'subspaceGener','groupStats','Cl1and2_p_rank_Stats.nii'))
save4Dnii(tStats_all_Cl_1to2_Vox,V,fullfile(root,'subspaceGener','groupStats','tStat_all_Cl_1to2.nii'))
save4Dnii(tStats_all_Cl_1to3_Vox,V,fullfile(root,'subspaceGener','groupStats','tStat_all_Cl_1to3.nii'))
save4Dnii(p_rank_StatsCl1,V,fullfile(root,'subspaceGener','groupStats','Cl1_p_rank_Stats.nii'))
save4Dnii(p_rank_StatsCl2,V,fullfile(root,'subspaceGener','groupStats','Cl2_p_rank_Stats.nii'))
