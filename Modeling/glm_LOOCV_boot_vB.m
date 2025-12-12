%% cmpBAModelAnalysis
% Perform analysis on the behavior/activity modeling

clear; close all; clc

% set folder
p1 = '/MATLAB Drive/FY2025/Behavior';
cd(p1)

% load groups
load('groupIDs.mat','groups')

% load and normalize data
load('/MATLAB Drive/FY2025/Behavior/data/model_vB/table_allRows_vB.mat','tbl','nTrials')
dataM = double(tbl.mouse);
varNamesAll = tbl.Properties.VariableNames;

% define genotype groups
useGroup = 1;
groupSaveName = {'All-LOO','WT-LOO','PS19-LOO','All-LTO'};
genoGroups = {{'0','1'},'0','1',{'0','1'}};
useGroupsName = {'All','WT','PS19'};

useGroupIdx = ismember(tbl.geno,genoGroups{useGroup});

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
useCols = [1:4 7 12 13 15 17 22 24:25 26]; % SG

% define lower bound
defCols = [24 25];
defCols = defCols(ismember(defCols,useCols));
if ~isempty(defCols)
    formLower = [varNamesAll{useCols(end)} ' ~ ' strjoin(varNamesAll(defCols),' + ')];
else
    formLower = 'y ~ 1';
end

nUseCols = length(useCols);
useTbl = tbl(useGroupIdx,useCols);
nFxd = find(strcmp(useTbl.Properties.VariableNames,'y'))-1;

% remove nans and normalize
keepIdx = all(~ismissing(useTbl),2);
useTbl = useTbl(keepIdx,:);
useY = useTbl.y;

% get mouse indices for use group
useTblM = dataM(useGroupIdx);
useTblM = useTblM(keepIdx);
useMVals = unique(useTblM);
nMVals = length(useMVals);
nRows = length(useTblM);

% define bootstrap parameters
rng(42)
nRnd = 1000;
rndIdx = zeros(nRows,nRnd);
for rr = 1:nRnd
    rndIdx(:,rr) = randperm(nRows);
end

newPlot = 0;

% define data output
if newPlot==1 || ~isfile('data/bootDataTemp.mat')
    YTestFinal = cell(nMVals,nRnd+1);
    YTestPredFinal = cell(nMVals,nRnd+1);
    R2Final = zeros(nMVals,nRnd+1);
    rFinal = zeros(nMVals,nRnd+1);
    nVFinal = zeros(nMVals,nRnd+1);
    complete = false(nMVals,nRnd+1);
else
    load('data/bootDataTemp.mat')
end

% initialize parallel pool if none exists
% if isempty(gcp('nocreate'))
%     parpool;
% end


%% Run analysis

tic
for ii = 1:nMVals
    disp(ii)

    % define remove indices
    mRemove = useMVals(ii);
    idxRemove = ismember(useTblM,mRemove);
    
    % define train table, y, and indices
    trainTbl = useTbl(~idxRemove,:);
    trainY = useTbl.y(~idxRemove);
    trainRnd = [rndIdx(~idxRemove,:) nan(sum(~idxRemove),1)];

    % define test table, y, and indices
    testTbl = useTbl(idxRemove,:);
    testY = useTbl.y(idxRemove);
    testRnd = [rndIdx(idxRemove,:) nan(sum(idxRemove),1)];

    for rr = 1:nRnd+1
        % skip completed comparisons
        if complete(ii,rr)
            continue
        end

        % shuffle behavior
        if rr<=nRnd
            trainTbly = useY(trainRnd(:,rr));
            testTbly = useY(testRnd(:,rr));
        else
            trainTbly = trainY;
            testTbly = testY;
        end

        % create full train and test table
        trainTblFull = trainTbl;
        trainTblFull.y = trainTbly;
        testTblFull = testTbl;
        testTblFull.y = testTbly;

        % run full model selection
        mdlFinal = stepwiseglm(trainTblFull,'linear','Upper','linear',...
            'Lower',formLower,'NSteps',nUseCols-2,...
            'Distribution','binomial','BinomialSize',nTrials,'Verbose',0);

        % identify final variables
        idxGlob = find(mdlFinal.Formula.InModel);

        % test final model on test data
        yTestPred = predict(mdlFinal,testTblFull)*nTrials;
        curR2 = findR2(testTbly,yTestPred);
        curr = corr(testTbly,yTestPred);

        % record test R2
        YTestFinal{ii,rr} = testTbly;
        YTestPredFinal{ii,rr} = yTestPred;
        R2Final(ii,rr) = curR2;
        rFinal(ii,rr) = curr;
        nVFinal(ii,rr) = length(idxGlob);

        % save temporary file
        complete(ii,rr) = true;
        if mod(rr,100)==0
            save('data/bootDataTemp.mat','YTestFinal','YTestPredFinal','R2Final','rFinal','nVFinal','complete')
        end
    end
end
toc

% process concatenated R2
R2Cat = zeros(nRnd+1,1);
for rr = 1:nRnd+1
    YTestCat = cat(1,YTestFinal{:,rr});
    YTestPredCat = cat(1,YTestPredFinal{:,rr});
    R2Cat(rr) = findR2(YTestCat,YTestPredCat);
end
R2CatShuff = R2Cat(1:nRnd);
R2CatData = R2Cat(end);

R2CatLess = sum(R2CatShuff>R2CatData);
R2CatEqual = sum(R2CatShuff==R2CatData);
R2CatTile = (R2CatLess+R2CatEqual/2+1)/(nRnd+1);

% process individual R2
R2Mean = mean(R2Final,1,'omitnan');
R2Shuff = R2Mean(1:nRnd);
R2Data = R2Mean(end);

R2Less = sum(R2Shuff>R2Data);
R2Equal = sum(R2Shuff==R2Data);
R2Tile = (R2Less+R2Equal/2+1)/(nRnd+1);

% process individual r
rMean = mean(rFinal,1,'omitnan');
rShuff = rMean(1:nRnd);
rData = rMean(end);

rLess = sum(rShuff>rData);
rEqual = sum(rShuff==rData);
rTile = (rLess+rEqual/2+1)/(nRnd+1);


%% Save results

% initialize structure
bootResults = struct();
bootResults.tbl = useTbl;
bootResults.vars = useTbl.Properties.VariableNames;
bootResults.group = useGroupsName{useGroup};

% store all shuffle/data results
bootResults.nVars = nVFinal;
bootResults.R2.all = R2Final;
bootResults.r.all = rFinal;

% store mean results
bootResults.R2Cat.mean = R2Cat;
bootResults.R2.mean = R2Mean;
bootResults.r.mean = rMean;

% store p-value
bootResults.R2Cat.p = R2CatTile;
bootResults.R2.p = R2Tile;
bootResults.r.p = rTile;

% save results
save('data/model_vB/LOOCV-boot_vB.mat','bootResults')


%% Plot p-value figure

clear; close all; clc

cd('D:\AD_Project\Behavior')
load('data/model_vB/LOOCV-boot_vB.mat','bootResults')

pVar = 'R2Cat';
dataShuff = bootResults.(pVar).mean(1:end-1);
dataShuff(dataShuff==-Inf) = NaN;
dataReal = bootResults.(pVar).mean(end);
pUse = bootResults.(pVar).p;

figure; hold on
if strcmp(pVar,'r')
    ksdensity(dataShuff,-1:0.01:1);
else
    ksdensity(dataShuff);
end
xline(dataReal)
title(['Boostrap histogram: ' pVar ', p=' num2str(pUse,1)])
ylim([0 10])
% xlim([-0.4 0.4])

savefig('data/model_vB/LOOCV-boot_newDecodeIO.fig')

