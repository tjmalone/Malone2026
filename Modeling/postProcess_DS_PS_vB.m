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
idxY = 3;
nTrials = 19;

% define cell types
cellTypes = {'allMorph','ste','pyr'};
nCTypes = length(cellTypes);

% define genotype groups
load('groupIDs.mat','groups','sexes')
genoGroups = {{'0','1'},'0','1',{'0','1'}};
useGroup = 1;

% define cell limit
cellLimit = 5;

% define table use columns
% interLapCorr, spatialSelectivity, speedScoreAbs, dfofSig, 1, 2, 3, 4
% dfofIORaw, interMouseCorr: 5, 6
% Fixed categorical (sex, geno): 7, 8
% Response: 9
% Random categorical (mouse, day): 10, 11

% define model variables (with decode)
colsFxd = [1:5 7:8];
colY = 9;
colsFull = [colsFxd, colY];


%% Determine formula and NaN rows

for cc = 1:nCTypes
    % load data
    load(['data/model_vB/corrDataDS_' cellTypes{cc} '_vB.mat'],'corrData','corrLabels','mouseSubSamples')

    % run model for all subsamples
    nRnd = size(corrData,4);

    % generate subsample table
    tbl = genModelTable(corrData(:,:,:,1),corrLabels,idxY,nTrials,groups,sexes,useDays);

    % define full table and remove nan rows
    useGroupIdx = ismember(tbl.geno,genoGroups{useGroup});
    tblFull = tbl(useGroupIdx,colsFull);
    keepIdx = all(~ismissing(tblFull),2);

    % remove mice with insufficient cell numbers
    rmvMice = categorical(find(mouseSubSamples(:,1)<cellLimit));
    keepIdx = keepIdx & ~ismember(tbl.mouse,rmvMice);

    if cc==1
        % initialize NaN rows
        globalIdxKeep = keepIdx;

        % define model formula
        labsFull = string(tbl.Properties.VariableNames);

        formFinal = labsFull(colY) + ' ~ ' + strjoin(labsFull(colsFxd),' + ');
    else
        % update NaN rows
        globalIdxKeep = globalIdxKeep & keepIdx;
    end
end


%% Run model for all cell types

% initialize output struct
modelData = struct();
lossNames = {'r','R2'};
typeNames = {'Train','Test'};
for ii = lossNames
    for jj = typeNames
        modelData.(ii{:}).(jj{:}) = cell(1,nCTypes);
    end
end

% initialize parallel pool if none exists
if isempty(gcp('nocreate'))
    parpool;
end

for cc = 1:nCTypes
    % load data
    load(['data/model_vB/corrDataDS_' cellTypes{cc} '_vB.mat'],'corrData','corrLabels')

    % run model for all subsamples
    nRnd = size(corrData,4);
    for ii = lossNames
        for jj = typeNames
            modelData.(ii{:}).(jj{:}){cc} = zeros(nRnd,1);
        end
    end

    for rr = 1:nRnd
        % print progress
        progress = ['Progress: ' num2str(cc) '_' num2str(rr)];
        fprintf(repmat('\b', 1,20));  % backspace to clear message
        fprintf('%-20s',progress);  % pad to fixed width (e.g., 20 chars)

        % generate subsample table
        tbl = genModelTable(corrData(:,:,:,rr),corrLabels,idxY,nTrials,groups,sexes,useDays);

        % define full table and remove nan rows
        useGroupIdx = ismember(tbl.geno,genoGroups{useGroup});
        tblFull = tbl(useGroupIdx,colsFull);
        tblFull = tblFull(globalIdxKeep,:);

        % get mouse indices for use group
        dataM = double(tbl.mouse);
        useTblM = dataM(useGroupIdx);
        useTblM = useTblM(globalIdxKeep);
        tblUnq = 1:length(useTblM);

        % define iteration indices indices
        idxMouse = unique(useTblM);
        nItr = length(idxMouse);


        %% Calculate sub-model statistics

        % define data output
        yTrain = cell(nItr,1);
        yPredTrain = cell(nItr,1);
        yTest = cell(nItr,1);
        yPredTest = cell(nItr,1);
        unqTrain = cell(nItr,1);
        unqTest = cell(nItr,1);

        parfor ii = 1:nItr
            % define remove indices
            mRemove = idxMouse(ii,:);
            idxRemove = ismember(useTblM,mRemove);

            % define train table, y, and indices
            tblTrain = tblFull(~idxRemove,:);
            tblTest = tblFull(idxRemove,:);

            % fit model
            mdl = fitglm(tblTrain,formFinal,'Distribution','binomial','BinomialSize',nTrials);

            % find base model correlation
            yTrain{ii} = tblTrain.y;
            yTest{ii} = tblTest.y;
            yPredTrain{ii} = predict(mdl,tblTrain)*nTrials;
            yPredTest{ii} = predict(mdl,tblTest)*nTrials;
            unqTrain{ii} = tblUnq(~idxRemove)';
            unqTest{ii} = tblUnq(idxRemove)';
        end

        curData.yTrain = yTrain;
        curData.yTest = yTest;
        curData.yPredTrain = yPredTrain;
        curData.yPredTest = yPredTest;
        curData.unqTrain = unqTrain;
        curData.unqTest = unqTest;


        %% Calculate concatenated model fit

        % loop through train/test
        for dt = 1:length(typeNames)
            % get current results
            curDType = typeNames{dt};
            curYData = cat(1,curData.(['y' curDType]){:});
            curYPred = cat(1,curData.(['yPred' curDType]){:});
            curUnq = cat(1,curData.(['unq' curDType]){:});

            % average matched predictions
            uUnq = unique(curUnq);
            nUnq = length(uUnq);
            meanYData = zeros(nUnq,1);
            meanYPred = zeros(nUnq,1);
            for qq = 1:nUnq
                meanYData(qq) = mean(curYData(curUnq==uUnq(qq)),'omitnan');
                meanYPred(qq) = mean(curYPred(curUnq==uUnq(qq)),'omitnan');
            end

            % calculate R2 values
            R2 = findR2(meanYData,meanYPred);

            % calculate r value
            r = corr(meanYData,meanYPred,'rows','complete');

            % store prediction values
            modelData.R2.(curDType){cc}(rr) = R2;
            modelData.r.(curDType){cc}(rr) = r;
        end
    end
end
fprintf('\n')


%% Save data

save(['data/model_vB/modelData_DS_PS_cLim-' num2str(cellLimit) '_vB.mat'],'modelData')

