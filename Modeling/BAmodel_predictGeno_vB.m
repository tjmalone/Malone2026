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

typeSets = {'All Common','Common vs. Grid vs. Non-Grid','Common vs. Stellate vs. Pyramidal'};
typeSetCode = {'common','GnG','PS'};
nTypeSets = length(typeSets);

% initialize data structure
LOOData = struct();

% set random seed
rng(42)

% initialize parallel pool if none exists
if isempty(gcp('nocreate'))
    parpool;
end


%% Perform modeling

for tSet = 1:nTypeSets
    if tSet==1
        % load data
        load('data/model_vB/table_allRows_vB.mat','tbl')
        allTable = struct();
        allTable.common = {tbl};

        % define final variables
        colsFinal = [4 7 13 15 17 22];
        colsFull = [1:4 7 12 13 15 17 22 26];

        % find missing data
        keepIdx = all(~ismissing(tbl(:,colsFull)),2);

        % define cell type
        cellTypes = {'common'};
    elseif tSet==2
        % load data
        load('data/model_vB/tableFull_DS_GnG_vB.mat','allTable','globalIdxKeep','cellTypes')

        % define final variables
        colsFinal = 1:6;
        colsFull = [1:6 10];

        % find missing data
        keepIdx = globalIdxKeep;
    else
        % load data
        load('data/model_vB/tableFull_DS_PS_vB.mat','allTable','globalIdxKeep','cellTypes')

        % define final variables
        colsFinal = 1:5;
        colsFull = [1:5 9];

        % find missing data
        keepIdx = globalIdxKeep;
    end


    %% Loop through cell types

    nCellTypes = length(cellTypes);
    for cc = 1:nCellTypes

        % initialize data structure
        nSamples = min(size(allTable.(cellTypes{cc}),1),100);
        LOOData.(typeSetCode{tSet}).(cellTypes{cc}).idxMouse = cell(nSamples,1);
        LOOData.(typeSetCode{tSet}).(cellTypes{cc}).tblM = cell(nSamples,1);
        LOOData.(typeSetCode{tSet}).(cellTypes{cc}).idxUse = cell(nSamples,1);
        LOOData.(typeSetCode{tSet}).(cellTypes{cc}).yTest = cell(nSamples,1);
        LOOData.(typeSetCode{tSet}).(cellTypes{cc}).yPredTest = cell(nSamples,1);

        % loop through sub samples
        for sample = 1:nSamples

            tblCur = allTable.(cellTypes{cc}){sample};

            % process column indices
            nFull = length(colsFull);
            idxFinal = find(ismember(colsFull,colsFinal));
            idxOther = setdiff(1:nFull-1,idxFinal);

            % define full table and remove nan rows
            tblFull = tblCur(:,colsFull);
            tblFull = tblFull(keepIdx,:);
            trueY = tblFull.y;
            nRows = size(tblFull,1);

            % get mouse indices for use group
            dataM = double(tblCur.mouse);
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
            LOOData.(typeSetCode{tSet}).(cellTypes{cc}).idxMouse{sample} = idxMouse;
            LOOData.(typeSetCode{tSet}).(cellTypes{cc}).tblM{sample} = useTblM;
            LOOData.(typeSetCode{tSet}).(cellTypes{cc}).idxUse{sample} = idxUse;
            LOOData.(typeSetCode{tSet}).(cellTypes{cc}).yTest{sample} = yTest;
            LOOData.(typeSetCode{tSet}).(cellTypes{cc}).yPredTest{sample} = yPredTest;
        end
    end
end

% save processed data
save(['data/model_vB/LOOCVplotMin_' svSupp '_vB.mat'],'LOOData')

