function H = pileSimilarityConSpN_RSAtoolbox(RowData,SPM)
% RowData: raw data from searchlight

% pre - whitening and calculate beta whights:
% default: normalize using all runs (overall):
% if wanted to be changed should add in the options: 'normmode':   
% 'overall': Does the multivariate noise normalisation overall (default)
% 'runwise': Does the multivariate noise normalisation by run
[betaAll00,resMS,Sw_hat,beta_hat,shrinkage,trRR]=noiseNormalizeBeta(RowData,SPM);%,varargin);

% % structure of distance matrix:
DistMs = creatDistSpNode();

% which betas: (for later analysis)
numPile= 10;
nMap = 4;
v = ones(1,numPile*nMap);
zeroBetas = 10 + 5 +2 +7;
bn0 =  [zeros(1,5) v zeros(1,zeroBetas)];%zeros in the end - betas for 5 maps,the catch, the start and the movement regressors + ventricles regressors
bn = [bn0 bn0 bn0 bn0];
nwb = length(bn);
b_idx0 = 1:nwb;
b_idx = b_idx0(bn==1);

betaAll0 = betaAll00(b_idx,:);

lvvoxels = length(betaAll0(1,:));
betaAll1 = reshape(betaAll0,[40,4,lvvoxels]);
betaAllLooAv = zeros(40,4,lvvoxels);

nr = 1:4;
%% LOO averages:
for r = 1:4
    betaAllLooAv(:,r,:) = mean(betaAll1(:,nr~=r,:),2);
end

c1L = zeros(4,4);
c2L = zeros(4,4);
c3L = zeros(4,4);
c4L = zeros(4,4);

nPile = 10;
for r = 1:4
    
    %averages:
    %hex 1
    conSHex1av = squeeze(betaAllLooAv(1:10,r,:));
    conSHex1av = conSHex1av - repmat(mean(conSHex1av),nPile,1);
    conSHex1av = conSHex1av - repmat(mean(conSHex1av,2),1,lvvoxels);
    
    % hex 2
    conSHex2av = squeeze(betaAllLooAv(11:20,r,:));
    conSHex2av = conSHex2av - repmat(mean(conSHex2av),nPile,1);
    conSHex2av = conSHex2av - repmat(mean(conSHex2av,2),1,lvvoxels);
    
    % cluster 1:
    conSC1av = squeeze(betaAllLooAv(21:30,r,:));
    conSC1av = conSC1av - repmat(mean(conSC1av),nPile,1);
    conSC1av = conSC1av - repmat(mean(conSC1av,2),1,lvvoxels); 
    
    % cluster 2:
    conSC2av = squeeze(betaAllLooAv(31:40,r,:));
    conSC2av = conSC2av - repmat(mean(conSC2av),nPile,1);  
    conSC2av = conSC2av - repmat(mean(conSC2av,2),1,lvvoxels); 
    
    % run patterns:
    %hex 1
    conSHex1r = betaAll0((r-1)*40+1:(r-1)*40+10,:);
    conSHex1r =  conSHex1r - repmat(mean(conSHex1r),nPile,1);
    conSHex1r =  conSHex1r - repmat(mean(conSHex1r,2),1,lvvoxels);
    
    %hex 2
    conSHex2r = betaAll0((r-1)*40+11:(r-1)*40+20,:);
    conSHex2r = conSHex2r - repmat(mean(conSHex2r),nPile,1);
    conSHex2r =  conSHex2r - repmat(mean(conSHex2r,2),1,lvvoxels);
    
    %cluster 1:
    conSC1r = betaAll0((r-1)*40+21:(r-1)*40+30,:);
    conSC1r = conSC1r - repmat(mean(conSC1r),nPile,1);
    conSC1r = conSC1r - repmat(mean(conSC1r,2),1,lvvoxels);
    
    %cluster 2:
    conSC2r = betaAll0((r-1)*40+31:(r-1)*40+40,:);
    conSC2r = conSC2r - repmat(mean(conSC2r),nPile,1);
    conSC2r = conSC2r - repmat(mean(conSC2r,2),1,lvvoxels);
    
    %correlations:
    % hex 1:
    
    CC1  = conSHex1av*conSHex1r';
    normM1 = sqrt(diag(CC1 )*diag(CC1)');
    CC1 = atanh(CC1./normM1);
    
    c1L(1,r) = mean(CC1(DistMs.HexBig==1));
    c1L(2,r) = mean(CC1(DistMs.HexBig==2));
    c1L(3,r) = mean(CC1(DistMs.HexBig==3));
    c1L(4,r) = mean(CC1(DistMs.HexBig==4));
    
    
   % To hex 2:
  
    CC2  = conSHex2av*conSHex2r';
    normM2 = sqrt(diag(CC2 )*diag(CC2)');
    CC2 = atanh(CC2./normM2);
    
    c2L(1,r) = mean(CC2(DistMs.HexSmal==1));
    c2L(2,r) = mean(CC2(DistMs.HexSmal==2));
    c2L(1,r) = mean(CC2(DistMs.HexSmal==3));
    c2L(1,r) = mean(CC2(DistMs.HexSmal==4));
    
    
    % cluster 1:
   
    CC3  = conSC1av*conSC1r';
    normM3 = sqrt(diag(CC3)*diag(CC3)');
    CC3 = atanh(CC3./normM3);
    
    c3L(1,r) = mean(CC3(DistMs.Clust==1));
    c3L(2,r) = mean(CC3(DistMs.Clust==2));
    c3L(3,r) = mean(CC3(DistMs.Clust==3));
    c3L(4,r) = mean(CC3(DistMs.Clust==4));
    
    % cluster 2
    CC4  = conSC2av*conSC2r';
    normM4 = sqrt(diag(CC4)*diag(CC4)');
    CC4 = atanh(CC4./normM4);
    
    c4L(1,r) = mean(CC4(DistMs.Clust==1));
    c4L(2,r) = mean(CC4(DistMs.Clust==2));
    c4L(3,r) = mean(CC4(DistMs.Clust==3));
    c4L(4,r) = mean(CC4(DistMs.Clust==4));   
end

ms1 = mean(c1L,2);
ms2 = mean(c2L,2);
ms3 = mean(c3L,2);
ms4 = mean(c4L,2);


H = [ms1;ms2;ms3;ms4];

