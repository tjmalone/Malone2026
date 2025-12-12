%% cmpBAModelAnalysis
% Perform analysis on the behavior/activity modeling

clear; close all; clc

% set folder
p1 = 'D:\AD_Project\Behavior';
cd(p1)

% load genotypes and sexes
load('groupIDs.mat')

% define use days
useDays = 2:11;

% define genotype groups
useGroups = {{'0','1'},'0','1'};
useGroupsName = {'All','WT','PS19'};

cellTypes = {'ste','pyr'};

for cc = 1:length(cellTypes)
    % load data
    load(['data/model_vB/corrData_' cellTypes{cc} '_PS_vB.mat'],'corrData','corrLabels')

    % define group variables
    mouseIdent = repmat((1:size(corrData,1))',1,size(corrData,2));
    mouseSex = zeros(size(mouseIdent));
    mouseSex(sexes{3},:) = 1;
    mouseGeno = zeros(size(mouseIdent));
    mouseGeno(groups{2},:) = 1;
    mouseDay = repmat(1:size(corrData,2),size(corrData,1),1);

    % set regression parameters
    paramsFxd = 4:8;
    paramsY = 1;
    nFxd = length(paramsFxd);
    nY = length(paramsY);

    % define factor names
    varNamesFxd = corrLabels(paramsFxd);
    varNamesY = {'y'};

    % load current data
    dataFxd = reshape(corrData(:,useDays,paramsFxd),[],nFxd);
    dataY = reshape(corrData(:,useDays,paramsY),[],nY);

    % normalize data and y
    dataFxdNorm = normalize(dataFxd);
    nTrials = 19;
    dataYCount = dataY/100*nTrials;

    % define categorical variables
    dataM = categorical(reshape(mouseIdent(:,useDays),[],nY));
    dataSex = categorical(reshape(mouseSex(:,useDays),[],nY));
    dataGeno = categorical(reshape(mouseGeno(:,useDays),[],nY));
    dataDay = categorical(reshape(mouseDay(:,useDays),[],nY));

    % define categorical labels
    catLabsFxd = {'fSex','fGeno'};
    catLabsRnd = {'mouse','day'};

    % generate data table
    tblA = array2table(dataFxdNorm,'VariableNames',varNamesFxd);
    tblB = table(dataSex,dataGeno,dataYCount,dataM,dataDay,...
        'VariableNames',[catLabsFxd,varNamesY(:),catLabsRnd]);
    tbl = [tblA tblB];
    varNamesAll = tbl.Properties.VariableNames;
    save(['data/model_vB/table_' cellTypes{cc} '_PS_vB.mat'],'tbl','nTrials')
end
