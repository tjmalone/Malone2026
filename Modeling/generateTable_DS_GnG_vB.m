%% postProcess_DS
% Process data using a final model determined elsewhere. This
% version processes both training and test data, glm, and LOOCV for common,
% grid, and non-grid cells. The relative training and test fits are plotted
% for comparison.
%

%% Load inputs and define final model

clear; close all; clc

% set folder
p1 = 'D:\AD_Project\Behavior';
cd(p1)

% define model parameters days
useDays = 2:11;
idxY = 1;
nTrials = 19;

% define cell types
cellTypes = {'common','grid','nongrid'};
nCTypes = length(cellTypes);

% define genotype groups
load('groupIDs.mat','groups','sexes')
genoGroups = {{'0','1'},'0','1',{'0','1'}};
useGroup = 1;

% define table use columns
% interLapCorr, spatialSelectivity, speedScoreAbs, dfofSig, 1, 2, 3, 4
% decodeIORaw, dfofIORaw, interMouseCorr: 5, 6, 7
% Fixed categorical (sex, geno): 8, 9
% Response: 10
% Random categorical (mouse, day): 11, 12

% define model variables (with decode)
colsFxd = [1:6 8:9];
colY = 10;
colsFull = [colsFxd, colY];


%% Determine formula and NaN rows

for cc = 1:nCTypes
    % load data
    load(['data/model_vB/corrDataDS_' cellTypes{cc} '_vB.mat'],'corrData','corrLabels')

    % run model for all subsamples
    nRnd = size(corrData,4);

    % generate subsample table
    tbl = genModelTable(corrData(:,:,:,1),corrLabels,idxY,nTrials,groups,sexes,useDays);

    % define full table and remove nan rows
    useGroupIdx = ismember(tbl.geno,genoGroups{useGroup});
    tblFull = tbl(useGroupIdx,colsFull);
    keepIdx = all(~ismissing(tblFull),2);

    if cc==1
        % initialize NaN rows
        globalIdxKeep = keepIdx;
    else
        % update NaN rows
        globalIdxKeep = globalIdxKeep & keepIdx;
    end
end


%% Run model for all cell types

allTable = struct();

for cc = 1:nCTypes
    % load data
    load(['data/model_vB/corrDataDS_' cellTypes{cc} '_vB.mat'],'corrData','corrLabels')

    % run model for all subsamples
    nRnd = size(corrData,4);
    allTable.(cellTypes{cc}) = cell(nRnd,1);
    for rr = 1:nRnd
        % generate subsample table
        allTable.(cellTypes{cc}){rr} =...
            genModelTable(corrData(:,:,:,rr),corrLabels,idxY,nTrials,groups,sexes,useDays);
    end
end

save('data/model_vB/tableFull_DS_GnG_vB.mat','allTable','globalIdxKeep','cellTypes')

