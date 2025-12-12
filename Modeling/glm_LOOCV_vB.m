%% cmpBAModelAnalysis
% Perform analysis on the behavior/activity modeling

clear; close all; clc

% set folder
p1 = 'D:\AD_Project\Behavior';
cd(p1)

% load genotypes and sexes
load('groupIDs.mat')

% load data
load('data/model_vB/corrData_vB.mat','corrData','corrLabels')

% define use days
useDays = 2:11;

% define genotype groups
useGroups = {{'0','1'},'0','1'};
useGroupsName = {'All','WT','PS19'};

% define group variables
mouseIdent = repmat((1:size(corrData,1))',1,size(corrData,2));
mouseSex = zeros(size(mouseIdent));
mouseSex(sexes{3},:) = 1;
mouseGeno = zeros(size(mouseIdent));
mouseGeno(groups{2},:) = 1;
mouseDay = repmat(1:size(corrData,2),size(corrData,1),1);

% set regression parameters
paramsFxd = 3:25;
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
dataMDay = dataM.*dataDay;

% define categorical labels
catLabsFxd = {'fSex','fGeno'};
catLabsRnd = {'mouse','sex','geno','day','mouseDay'};

% generate data table
tblA = array2table(dataFxdNorm,'VariableNames',varNamesFxd);
tblB = table(dataSex,dataGeno,dataYCount,dataM,dataSex,dataGeno,dataDay,dataMDay,...
    'VariableNames',[catLabsFxd,varNamesY(:),catLabsRnd]);
tbl = [tblA tblB];
varNamesAll = tbl.Properties.VariableNames;
save('data/model_vB/table_allRows_vB.mat','tbl','nTrials')

% define table use columns
% meanVel: 1
% interDayCorr, interMouseCorr, interLap Corr, corrNear, corrFar: 2, 3, 4, 5, 6
% spatialSelectivity, fldNum, fldWidth, fldCoverage, spatialInfo: 7, 8, 9, 10, 11
% speedScoreRaw, speedScoreRaw, selfRunCorr (NaN), dfofSig: 12, 13, 14, 15 
% decoding (NaN), decodeIORaw (NaN), decodeIONorm (NaN): 16, 17, 18
% ampDiffPerNL, ampDiffIORaw, ampDiffIONorm: 19, 20, 21
% dfofIORaw, dfofIOAbs: 22, 23
% Fixed categorical (sex, geno): 24, 25
% Response: 26
% Random categorical (mouse, day, mouse/day): 27, 30, 31

% define stable columns
defCols = 24:25;
defInt = 0;

% useCols = [1 2:4 7 10 12 13 15 17 22 defCols 26];

% remove field coverage
useCols = [1 2:4 7 12 13 15 17 22 defCols 26];

defColsUse = defCols(ismember(defCols,useCols));

if length(defCols)~=2 || defInt~=1
    useInt = 0;

    % define lower bound
    if ~isempty(defColsUse)
        formLower = [varNamesAll{useCols(end)} ' ~ ' strjoin(varNamesAll(defColsUse),' + ')];
    else
        formLower = 'y ~ 1';
    end

    % define upper bound
    formUpper = [varNamesAll{useCols(end)} ' ~ ' strjoin(varNamesAll(useCols(1:end-1)),' + ')];
else
    useInt = 1;

    % define lower bound
    formLower = [varNamesAll{useCols(end)} ' ~ ' strjoin(varNamesAll(defColsUse),'*')];

    % define upper bound
    formUpper = [varNamesAll{useCols(end)} ' ~ ' strjoin(varNamesAll(useCols(1:end-3)),' + ')...
        ' + ' strjoin(varNamesAll(defColsUse),'*')];
end

nUseCols = length(useCols);
nItr = size(corrData,1);
varsFinal = zeros(nItr,nUseCols-1);
mdlWeights = zeros(nItr,nUseCols-1);
mdlR2 = zeros(nItr,1);
useGroup = 1;

% define parameter columns
paramCols = useCols;
paramCols(useCols==26) = [];


%% Run analysis

for ii = 1:nItr
    disp(ii)

    % load table
    curTbl = tbl;
    curM = dataM;

    % keep only mice in current group
    idxRemoveG = ~ismember(curTbl.geno,useGroups{useGroup});
    curTbl(idxRemoveG,:) = [];
    curM(idxRemoveG,:) = [];
    preN = length(curM);

    % remove "leave-one-out" mouse
    idxRemoveLOO = ismember(curM,num2str(ii));
    if sum(idxRemoveLOO)==0
        varsFinal(ii,:) = NaN;
        continue
    end
    curTbl(idxRemoveLOO,:) = [];
    curM(idxRemoveLOO,:) = [];

    % set use columns
    curTbl = curTbl(:,useCols);
    nFxd = find(strcmp(curTbl.Properties.VariableNames,'y'))-1;

    % % run full model selection without forced confounds
    % mdlFinal = stepwiseglm(curTbl,formUpper,'Upper',formUpper,...
    %     'Distribution','binomial','BinomialSize',nTrials,'Verbose',1);

    % run full model selection with forced confounds
    mdlFinal = stepwiseglm(curTbl,formUpper,'Upper',formUpper,...
        'Lower',formLower,'NSteps',nUseCols-2,...
        'Distribution','binomial','BinomialSize',nTrials,'Verbose',1);

    % % build full model selection without forced confounds
    % mdlFinal = stepwiseglm(curTbl,'constant','Upper',formUpper,...
    %     'Distribution','binomial','BinomialSize',nTrials,'Verbose',1);

    % % build full model selection with forced confounds
    % mdlFinal = stepwiseglm(curTbl,formLower,'Upper',formUpper,...
    %     'Lower',formLower,...
    %     'Distribution','binomial','BinomialSize',nTrials,'Verbose',1);

    % identify final variables
    idxGlob = find(mdlFinal.Formula.InModel);
    if useInt==1 && sum(mdlFinal.Formula.Terms(end,:))==2
        idxGlob(end+1) = nUseCols;
    end

    % record final model variables
    varsFinal(ii,idxGlob) = 1;

    % record final model weights
    mdlWeights(ii,idxGlob) = mdlFinal.Coefficients{2:end,"Estimate"};

    % record final model R2
    mdlR2(ii) = mdlFinal.Rsquared.Adjusted;
end

varsMean = mean(varsFinal,1,'omitnan')*100;
weightsMean = mean(mdlWeights,1,'omitnan');
R2Mean = mean(mdlR2,'omitnan');

fprintf('Final Factor Prevalence:\n');
for ii = 1:nUseCols-1
    fprintf('%d. Factor %s - %.1f%%\n',ii,varNamesAll{paramCols(ii)},varsMean(ii));
end
