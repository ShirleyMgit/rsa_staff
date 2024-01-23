function H = pileSimilarityCorBetMaps(pileSimData)% Bet in the name is an error
% Here the input to the function is center of searchligh for pile
% similarity within each maps, the output is correlection between the
% similarity of piles between maps.

% RowData: raw data from searchlight, here it is searchlight pile
% similarity results from each map. so I should have 4 numbers for each
% voxel (This is an averaged quentety thereofore no runs).
disp(size(pileSimData));

corMat = pileSimData*pileSimData';
%corMat = pileSimData'*pileSimData;
n1 = sqrt(diag(corMat));
normM =n1*n1';
corMatallR = corMat./normM;

disp(size(corMatallR ));

% All entries of the matrix (4x4):
H = corMatallR(:);