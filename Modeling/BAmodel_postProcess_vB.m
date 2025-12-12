%% BAmodel_postProcess
% Process data using a final model determined elsewhere. This
% generic version processes both training and test data, using lme of
% glm, and using LOOCV or LTOCV. BAmodel_plotGen should be run after to
% plot results.
% 

%% Load inputs and define final model

clear; close all; clc

% set folder
p1 = '/MATLAB Drive/FY2025/Behavior';
cd(p1)

% load groups
load('groupIDs.mat','groups')

% load and normalize data
load('/MATLAB Drive/FY2025/Behavior/data/model_vB/table_allRows_vB.mat','tbl','nTrials')
dataM = double(tbl.mouse);

% define model type
modelType = 'glm'; % lme or glm
if strcmp(modelType,'lme')
    modelParam = 'ML';
elseif strcmp(modelType,'glm')
    modelParam = nTrials;
end

% define genotype groups
useGroup = 1;
groupSaveName = {'All-LOO','WT-LOO','PS19-LOO','All-LTO'};
genoGroups = {{'0','1'},'0','1',{'0','1'}};
useGroupIdx = ismember(tbl.geno,genoGroups{useGroup});

% define leave out type
if useGroup==4
    leaveType = 'LTO';
else
    leaveType = 'LOO';
end

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

% define final SG variables
colsFinal = [4 7 13 15 17 22 24 25]; % SG
colsFull = [1:4 7 12 13 15 17 22 24:25 26]; % SG

svSupp = '_SG';

% process column indices
nFinal = length(colsFinal);
nFull = length(colsFull);
idxFinal = find(ismember(colsFull,colsFinal));
idxOther = setdiff(1:nFull-1,idxFinal);

% define full table and remove nan rows
tblFull = tbl(useGroupIdx,colsFull);
keepIdx = all(~ismissing(tblFull),2);
tblFull = tblFull(keepIdx,:);

% generate NaN info table
tblNaN = tbl(useGroupIdx,[24:25 27 30]);
tblNaN.keepIdx = keepIdx;

% get mouse indices for use group
useTblM = dataM(useGroupIdx);
useTblM = useTblM(keepIdx);

% define iteration indices indices
if strcmp(leaveType,'LOO')
    idxMouse = unique(useTblM);
    nItr = length(idxMouse);
elseif strcmp(leaveType,'LTO')
    idxWT = repmat(groups{1}',length(groups{2}),1);
    idxAD = repelem(groups{2}',length(groups{1}),1);
    idxMouse = [idxWT idxAD];
    nItr = size(idxMouse,1);
end

% % initialize parallel pool if none exists
% if isempty(gcp('nocreate'))
%     parpool;
% end


%% Calculate sub-model statistics

% define data output
idxUse = cell(1,nItr);
mdl = cell(1,nItr);
R2 = zeros(nItr,1);
coeff = zeros(nItr,nFinal+1);
yTrain = cell(nItr,1);
yPredTrain = cell(nItr,1);
yPredRemTrain = cell(nItr,1);
yPredExpTrain = cell(nItr,1);
yTest = cell(nItr,1);
yPredTest = cell(nItr,1);
yPredRemTest = cell(nItr,1);
yPredExpTest = cell(nItr,1);
domR2Train = cell(nItr,1);
yPredTier = cell(nItr,1);

% define loop progress
nPer = 10;
progressMarks = round(linspace(1/nPer,1,nPer)*nItr);

for ii = 1:nItr
    if ismember(ii,progressMarks)
        fprintf('Progress: %3.0f%% (%d/%d)\n', ii/nItr*100, ii, nItr);
    end

    % define remove indices
    mRemove = idxMouse(ii,:);
    idxRemove = ismember(useTblM,mRemove);
    
    % define train table, y, and indices
    tblTrain = tblFull(~idxRemove,:);
    tblTest = tblFull(idxRemove,:);

    % run sub/expansion models
    mdlCur = stepwiseModel_post(tblTrain,tblTest,idxFinal,idxOther,nFull,modelType,nTrials,0);

    % perform dominance analysis for training data
    [domDataTrain,yPredTierCur] = dominanceAnalysis(tblTrain,mdlCur.formula,modelParam,modelType,tblTest);

    % store data
    idxUse{ii} = ~idxRemove;
    mdl{ii} = mdlCur.mdl;
    R2(ii) = mdlCur.mdl.Rsquared.Ordinary;
    coeff(ii,:) = mdlCur.mdl.Coefficients.Estimate;
    yTrain{ii} = mdlCur.yTrain;
    yPredTrain{ii} = mdlCur.yPredTrain;
    yPredRemTrain{ii} = mdlCur.yPredRemTrain;
    yPredExpTrain{ii} = mdlCur.yPredExpTrain;
    yTest{ii} = mdlCur.yTest;
    yPredTest{ii} = mdlCur.yPredTest;
    yPredRemTest{ii} = mdlCur.yPredRemTest;
    yPredExpTest{ii} = mdlCur.yPredExpTest;
    domR2Train{ii} = domDataTrain;
    yPredTier{ii} = yPredTierCur;
end

% store structure data
LOOData = struct();
LOOData.idxMouse = idxMouse;
LOOData.tblM = useTblM;
LOOData.modelType = modelType;
LOOData.leaveType = leaveType;
LOOData.idxUse = idxUse;
LOOData.mdl = mdl;
LOOData.R2 = R2;
LOOData.coeff = coeff;
LOOData.yTrain = yTrain;
LOOData.yPredTrain = yPredTrain;
LOOData.yPredRemTrain = yPredRemTrain;
LOOData.yPredExpTrain = yPredExpTrain;
LOOData.yTest = yTest;
LOOData.yPredTest = yPredTest;
LOOData.yPredRemTest = yPredRemTest;
LOOData.yPredExpTest = yPredExpTest;
LOOData.domR2Train = domR2Train;
LOOData.yPredTier = yPredTier;

fprintf('\nMean R2 = %.2f\n',mean(LOOData.R2(:,1)));


%% Process sub/expansion data

dataProcessed = struct();
dataTypes = {'Train','Test'};

for dt = 1:length(dataTypes)
    % get current results
    curDType = dataTypes{dt};
    curYData = LOOData.(['y' curDType]);
    curYPred = LOOData.(['yPred' curDType]);
    curYPredRem = LOOData.(['yPredRem' curDType]);
    curYPredExp = LOOData.(['yPredExp' curDType]);

    % calculate R2 values by iteration
    R2Pred = cellfun(@(x,y) findR2(x,y),curYData,curYPred,'UniformOutput',false);
    R2Pred = cat(1,R2Pred{:});
    dataCur.byMouse.R2Pred = R2Pred;
    R2Rem = cellfun(@(x,y) findR2(x,y),curYData,curYPredRem,'UniformOutput',false);
    R2Rem = cat(1,R2Rem{:});
    dataCur.byMouse.R2Rem = R2Rem;
    R2Exp = cellfun(@(x,y) findR2(x,y),curYData,curYPredExp,'UniformOutput',false);
    R2Exp = cat(1,R2Exp{:});
    dataCur.byMouse.R2Exp = R2Exp;

    % calculate r values by iteration
    rPred = cellfun(@(x,y) corr(x,y,'rows','complete'),curYData,curYPred,'UniformOutput',false);
    rPred = cat(1,rPred{:});
    dataCur.byMouse.rPred = rPred;
    rRem = cellfun(@(x,y) corr(x,y,'rows','complete'),curYData,curYPredRem,'UniformOutput',false);
    rRem = cat(1,rRem{:});
    dataCur.byMouse.rRem = rRem;
    rExp = cellfun(@(x,y) corr(x,y,'rows','complete'),curYData,curYPredExp,'UniformOutput',false);
    rExp = cat(1,rExp{:});
    dataCur.byMouse.rExp = rExp;

    % calculate p values
    [~,dataCur.R2.pRem] = ttest(repmat(R2Pred,1,size(R2Rem,2)),R2Rem,'Tail','right');
    [~,dataCur.R2.pExp] = ttest(R2Exp,repmat(R2Pred,1,size(R2Exp,2)),'Tail','right');
    [~,dataCur.r.pRem] = ttest(repmat(rPred,1,size(rRem,2)),rRem,'Tail','right');
    [~,dataCur.r.pExp] = ttest(rExp,repmat(rPred,1,size(rExp,2)),'Tail','right');

    % concatenate predictions
    dataCur.yData = cat(1,curYData{:});
    dataCur.yPred = cat(1,curYPred{:});
    dataCur.yPredRem = cat(1,curYPredRem{:});
    dataCur.yPredExp = cat(1,curYPredExp{:});

    % calculate R2 values
    dataCur.R2.Pred = findR2(dataCur.yData,dataCur.yPred);
    dataCur.R2.Rem = findR2(dataCur.yData,dataCur.yPredRem);
    dataCur.R2.Exp = findR2(dataCur.yData,dataCur.yPredExp);

    % calculate r value
    dataCur.r.Pred = corr(dataCur.yData,dataCur.yPred,'rows','complete');
    dataCur.r.Rem = corr(dataCur.yData,dataCur.yPredRem,'rows','complete');
    dataCur.r.Exp = corr(dataCur.yData,dataCur.yPredExp,'rows','complete');

    % store data
    dataProcessed.(lower(curDType)) = dataCur;
end


%% Calculate training dominance

dTypes = {'domInd','domPart','domInter','domTot'};
domTrain = struct();

for dt = 1:length(dTypes)
    % get current dom field name
    curDType = dTypes{dt};

    % concatnate and take mean of dominance
    domCat = cellfun(@(x) x.(curDType),LOOData.domR2Train,'UniformOutput',false);
    domMean = mean(cat(1,domCat{:}),1,'omitnan');

    % store dominance mean
    domTrain.R2.(curDType) = domMean;
end

dataProcessed.domTrain = domTrain;


%% Calculate test dominance

% concatenate dominance tiers
yPredTierAll = LOOData.yPredTier;
nJ = numel(yPredTierAll{1});

yPredTierCat = cell(size(yPredTierAll{1}));

for ii = 1:nItr
    for jj = 1:nJ
        nK = size(yPredTierAll{ii}{jj},1);

        if ii==1
            yPredTierCat{jj} = cell(nK,1);
        end

        for kk = 1:nK
            szCur = size(yPredTierAll{ii}{jj}{kk,1},1);
            yPredTierCat{jj}{kk}(end+1:end+szCur,:) = cat(2,yPredTierAll{ii}{jj}{kk,:});
        end
    end
end

% calculate dominance R2 and r values
domR2 = cell(size(yPredTierCat));
domr = cell(size(yPredTierCat));
yTest = dataProcessed.test.yData;

for jj = 1:nJ
    nK = size(yPredTierCat{jj},1);

    for kk = 1:nK
        if mod(jj,size(yPredTierCat,1))~=1
            domR2{jj}(kk,1) = -diff(findR2(yTest,yPredTierCat{jj}{kk}));
            domr{jj}(kk,1) = -diff(corr(yTest,yPredTierCat{jj}{kk},'rows','complete'));
        else
            domR2{jj}(kk,1) = findR2(yTest,yPredTierCat{jj}{kk}(:,1));
            domr{jj}(kk,1) = corr(yTest,yPredTierCat{jj}{kk}(:,1),'rows','complete');
        end
    end
end

domTest = struct();

% calculate dominance tier means
domTest.R2.domTier = cellfun(@mean,domR2);
domTest.r.domTier = cellfun(@mean,domr);

% store individual dominance
domTest.R2.domInd = domTest.R2.domTier(1,:);
domTest.r.domInd = domTest.r.domTier(1,:);

% store partial dominance
domTest.R2.domPart = mean(domTest.R2.domTier(2:end-1,:),1,'omitnan');
domTest.r.domPart = mean(domTest.r.domTier(2:end-1,:),1,'omitnan');

% store interactional dominance
domTest.R2.domInter = domTest.R2.domTier(end,:);
domTest.r.domInter = domTest.r.domTier(end,:);

% store total dominance
domTest.R2.domTot = mean(domTest.R2.domTier,1);
domTest.r.domTot = mean(domTest.r.domTier,1);

dataProcessed.domTest = domTest;


%% Save processed data

dataProcessed.coeff = mean(LOOData.coeff,1);
dataProcessed.labsFinal = string(tblFull.Properties.VariableNames(idxFinal));
dataProcessed.labsOther = string(tblFull.Properties.VariableNames(idxOther));
save(['/MATLAB Drive/FY2025/Behavior/data/model_vB/LOOCVplot_' modelType '_' groupSaveName{useGroup} svSupp '_vB.mat'],...
    'LOOData','dataProcessed','nTrials','tblNaN')

