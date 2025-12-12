%% BAmodel_postProcess
% Process data using a final model determined elsewhere. This
% generic version processes both training and test data, using lme of
% glm, and using LOOCV or LTOCV. BAmodel_plotGen should be run after to
% plot results.
%

%% Load inputs and define final model

clear; close all; clc

% set folder
p1 = 'D:\AD_Project\Behavior';
cd(p1)

% define model type and genotype groups
modelType = 'glm';
genoGroups = {'0','1'};
svSupp = '_noSG';
nTrials = 19;

cellTypes = {'common','grid','nongrid','ste','pyr'};
nTypes = length(cellTypes);

% initialize data structure
LOOData = struct();

% set random seed
rng(42)

% initialize parallel pool if none exists
if isempty(gcp('nocreate'))
    parpool;
end


%% Perform modeling

for cc = 1:nTypes
    if cc==1
        % load data
        load('data/model_vB/table_allRows_vB.mat','tbl')

        % define final variables
        colsFinal = [4 7 13 15 17 22];
        colsFull = [1:4 7 12 13 15 17 22 26];

    elseif ismember(cc,[2 3])
        % load data
        load(['data/model_vB/table_' cellTypes{cc} '_GnG_vB.mat'],'tbl','nTrials')

        % define final variables
        colsFinal = 1:6;
        colsFull = [1:6 10];

    elseif ismember(cc,[4 5])
        % load data
        load(['data/model_vB/table_' cellTypes{cc} '_PS_vB.mat'],'tbl','nTrials')

        % define final variables
        colsFinal = 1:5;
        colsFull = [1:5 8];

    else
        error('Invalid cell type')
    end

    % process column indices
    nFull = length(colsFull);
    idxFinal = find(ismember(colsFull,colsFinal));
    idxOther = setdiff(1:nFull-1,idxFinal);

    % define full table and remove nan rows
    tblFull = tbl(:,colsFull);
    keepIdx = all(~ismissing(tbl(:,colsFull)),2);
    tblFull = tblFull(keepIdx,:);
    trueY = tblFull.y;
    nRows = size(tblFull,1);

    % get mouse indices for use group
    dataM = double(tbl.mouse);
    useTblM = dataM(:);
    useTblM = useTblM(keepIdx);

    % define iteration indices indices
    idxMouse = unique(useTblM);
    nItr = length(idxMouse);


    %% Calculate sub-model statistics

    % define shuffle parameters
    nShuff = 200;
    shuffIdx = zeros(nRows,nShuff);
    for shuff = 1:nShuff
        shuffIdx(:,shuff) = randperm(nRows);
    end

    % define data output
    idxUse = cell(nItr,nShuff+1);
    mdl = cell(nItr,nShuff+1);
    R2 = zeros(nItr,nShuff+1);
    yTrain = cell(nItr,nShuff+1);
    yPredTrain = cell(nItr,nShuff+1);
    yTest = cell(nItr,1);
    yPredTest = cell(nItr,nShuff+1);

    % define loop progress
    nPer = 10;
    progressMarks = round(linspace(1/nPer,1,nPer)*nItr);

    for ii = 1:nItr
        % define remove indices
        mRemove = idxMouse(ii,:);
        idxRemove = ismember(useTblM,mRemove);

        % define train table, y, and indices
        trainTbl = tblFull(~idxRemove,:);
        trainY = tblFull.y(~idxRemove);
        trainRnd = [shuffIdx(~idxRemove,:) nan(sum(~idxRemove),1)];

        % define test table, y, and indices
        testTbl = tblFull(idxRemove,:);
        testY = tblFull.y(idxRemove);
        testRnd = [shuffIdx(idxRemove,:) nan(sum(idxRemove),1)];

        parfor shuff = 1:nShuff+1
            % shuffle behavior
            if shuff<=nShuff
                trainTbly = trueY(trainRnd(:,shuff));
                testTbly = trueY(testRnd(:,shuff));
            else
                trainTbly = trainY;
                testTbly = testY;
            end

            % create full train and test table
            trainTblFull = trainTbl;
            trainTblFull.y = trainTbly;
            testTblFull = testTbl;
            testTblFull.y = testTbly;

            % run sub/expansion models
            mdlCur = stepwiseModel_post(trainTblFull,testTblFull,...
                idxFinal,idxOther,nFull,modelType,nTrials,0,1);

            % store data
            idxUse{ii,shuff} = ~idxRemove;
            mdl{ii,shuff} = mdlCur.mdl;
            R2(ii,shuff) = mdlCur.mdl.Rsquared.Ordinary;
            yTrain{ii,shuff} = mdlCur.yTrain;
            yPredTrain{ii,shuff} = mdlCur.yPredTrain;
            yTest{ii,shuff} = mdlCur.yTest;
            yPredTest{ii,shuff} = mdlCur.yPredTest;
        end

        % track progress
        if ismember(ii,progressMarks)
            fprintf('Progress: %3.0f%% (%d/%d)\n', ii/nItr*100, ii, nItr);
        end
    end

    % store structure data
    LOOData.(cellTypes{cc}).idxMouse = idxMouse;
    LOOData.(cellTypes{cc}).tblM = useTblM;
    LOOData.(cellTypes{cc}).idxUse = idxUse;
    LOOData.(cellTypes{cc}).yTest = yTest;
    LOOData.(cellTypes{cc}).yPredTest = yPredTest;
end

% save processed data
save('data/model_vB/LOOCVpredictGeno_noDS_vB.mat','LOOData')

