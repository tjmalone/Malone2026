%% collateActivityData_new
% Collates activity data for correlation with behavior

clear; close all; clc

baseFolder = '/MATLAB Drive/FY2025/imagingData';
cd(baseFolder)

% load data
load('foldersLearning.mat')

% load distribution data
load('data/distActivity.mat','distActivity')

% load cell information
load('data/cellSelect.mat','cellSelect')
cellTypes = {'common','grid','nongrid'};

% determine downsample numbers
nSamples = 200;

% grid cell indices
cellsG = cellSelect.learn.allSex.grid.allMorph;
cellsG = combCells(cellsG.WT,cellsG.AD);
NCellsG = cellfun(@length,cellsG);

% non-grid cell indices
cellsNG = cellSelect.learn.allSex.nongrid.allMorph;
cellsNG = combCells(cellsNG.WT,cellsNG.AD);
NCellsNG = cellfun(@length,cellsNG);

% common cell indixes
cellsC = cellSelect.learn.allSex.common.allMorph;
cellsC = combCells(cellsC.WT,cellsC.AD);

% calculate subsample cell number
nCellsMG = NCellsG(1:2:end,:)+NCellsG(2:2:end,:);
nCellsMNG = NCellsNG(1:2:end,:)+NCellsNG(2:2:end,:);

NSubsampleM = min(nCellsMG,nCellsMNG);
nMice = size(NSubsampleM,1);
nFOV = nMice*2;
nDays = size(NSubsampleM,2);

% process subsample number by mouse
mouseSubSamples =  NSubsampleM;

% select subsamples
rng(42)
cellsSmp = cell(1,3);
cellsAll = {cellsC,cellsG,cellsNG};

for ct = 1:3
    cellsSmp{ct} = cell(nFOV,nDays,nSamples);
    cellsCur = cellsAll{ct};

    for dd = 1:nDays
        for mm = 1:nMice
            % skip empty mice
            if NSubsampleM(mm,dd)==0; continue; end

            % generate combined indices
            idx = cell(1,2);
            idxMap = [];
            idxMark = [];
            for ff = 1:2
                curF = mm*2-2+ff;
                idx{ff} = cellsCur{curF,dd};
                idxMap = [idxMap; (1:length(idx{ff}))'];
                idxMark = [idxMark; ff*ones(length(idx{ff}),1)];
            end

            combN = length(idxMap);

            % subsample cells
            for rr = 1:nSamples
                if NSubsampleM(mm,dd)==1
                    if combN==1
                        % use single cell if only one cell exists
                        outIdx = 1;
                    else
                        % sample a single cell if N is 1
                        outIdx = randsample(1:combN,1);
                    end
                else
                    % sample N-1 cells if N is greater than 1
                    outIdx = randsample(1:combN,NSubsampleM(mm,dd)-1);
                end

                % determine FOV indices
                for ff = 1:2
                    curF = mm*2-2+ff;

                    % translate mouse indices to FOV indices
                    curOutIdx = sort(idxMap(outIdx(idxMark(outIdx)==ff)));

                    % store FOV indices
                    cellsSmp{ct}{curF,dd,rr} = idx{ff}(curOutIdx);
                end
            end
        end
    end
end

samplesAll = struct();
samplesAll.C.cellsSmp = cellsSmp{1};
samplesAll.C.nSamples = nSamples;
samplesAll.G.cellsSmp = cellsSmp{2};
samplesAll.G.nSamples = nSamples;
samplesAll.NG.cellsSmp = cellsSmp{3};
samplesAll.NG.nSamples = nSamples;

sampleTypes = fieldnames(samplesAll);
nSTypes = length(sampleTypes);


%% Calculate dfof in-cue/out-cue ratio

% initialize output structure
for ss = 1:nSTypes
    samplesAll.(sampleTypes{ss}).dfofIO = cell(nFOV,nDays,samplesAll.(sampleTypes{ss}).nSamples);
end

% loop through fov and day
for ff = 1:nFOV
    % print progress
    progress = ['Progress: ' num2str(ff)];
    fprintf(repmat('\b', 1,20));  % backspace to clear message
    fprintf('%-20s',progress);  % pad to fixed width (e.g., 20 chars)

    for dd = 1:nDays
        % correct for "true" day
        trueInv = find(trueDays(ff,:)==dd);
        if isempty(trueInv)
            continue
        end

        %% Load data and process template

        % change to current directory
        curFolderName = foldersLearning{ff}{trueInv};
        newFolderName = replace(curFolderName,'D:\AD_Project\imagingData\data','/MATLAB Drive/FY2025/SubData');
        newFolderName = replace(newFolderName,'\','/');
        cd(newFolderName)

        % load binned activity and cue template
        d = dir('dfofaveragesmooth_sig*');
        load(d(1).name)
        load('cueAnalysis_sig/tempRL.mat','tempRL')

        for ss = 1:nSTypes
            nRCur = samplesAll.(sampleTypes{ss}).nSamples;

            for rr = 1:nRCur
                % identify current cells
                curCells = samplesAll.(sampleTypes{ss}).cellsSmp{ff,dd,rr};

                % calculate ridge background ratio
                if isempty(curCells)
                    dfofIOCur = [];
                else
                    % get current activity data
                    curDfof = mean(dfofaveragesmooth_sig(:,curCells),2,'omitnan');

                    % perform quantification
                    dataQuant = distrQuant(curDfof',tempRL,10:20,2,[],1:2);

                    I = dataQuant.CueIn;
                    O = dataQuant.CueOut;
                    dfofIOCur = I/O;
                    if dfofIOCur==inf
                        dfofIOCur = NaN;
                    end
                end

                % store correlation values
                samplesAll.(sampleTypes{ss}).dfofIO{ff,dd,rr} = dfofIOCur;
            end
        end
    end
end

cd(baseFolder)


%% Calculate inter-mouse correlation

% initialize output structure
for ss = 1:nSTypes
    samplesAll.(sampleTypes{ss}).interMouseCorr = zeros(nFOV,nDays,samplesAll.(sampleTypes{ss}).nSamples);
end

% define averaging parameters
fovPerMouse = 2;

for ss = 1:nSTypes
    nRCur = samplesAll.(sampleTypes{ss}).nSamples;

    for rr = 1:nRCur
        % print progress
        progress = ['Progress: ' num2str(ss) '-' num2str(rr)];
        fprintf(repmat('\b', 1,20));  % backspace to clear message
        fprintf('%-20s',progress);  % pad to fixed width (e.g., 20 chars)

        % identify current data and cells
        distCur = distActivity.learn.allRuns.dfof.all;
        curCells = samplesAll.(sampleTypes{ss}).cellsSmp(:,:,rr);

        % initialize quantification struct
        dataInterMouseMean = zeros(nFOV/2,nDays);

        % correct empty cells
        maxIdx = cellfun(@max,curCells,'UniformOutput',false);
        emptyCellSel = cellfun(@isempty,maxIdx);
        emptyData = cellfun(@isempty,distCur);
        fixIdx = find(~emptyCellSel & emptyData);
        refIdx = find(~emptyData,1);
        szXCur = size(distCur{refIdx},2);
        for ii = 1:length(fixIdx)
            distCur{fixIdx(ii)} = nan(maxIdx{fixIdx(ii)},szXCur);
        end

        % define use data
        useData = cellfun(@(x,y) x(y,:),distCur,curCells,'UniformOutput',false);

        % fill in empty cells
        emptyAll = cellfun(@isempty,useData);
        for ff = 1:nFOV
            % get all empty and first full index
            curEmpty = find(emptyAll(ff,:));
            curFull = find(~emptyAll(ff,:),1);

            % skip full or empty rows
            if isempty(curEmpty) || isempty(curFull)
                continue
            end

            % fill empty cells
            for ii = curEmpty
                useData{ff,ii} = nan(size(useData{ff,curFull}));
            end
        end

        % fill in remaining empty cells
        curEmpty = find(cellfun(@(x) size(x,2),useData)==0);
        for ii = curEmpty'
            useData{ii} = zeros(0,szXCur);
        end

        % take mouse mean
        dataMouse = cell(nFOV/fovPerMouse,nDays);
        for ff = 1:fovPerMouse:nFOV
            curFOV = ff:ff+fovPerMouse-1;

            for dd = 1:nDays
                dataMouse{(ff+fovPerMouse-1)/fovPerMouse,dd}...
                    = cat(1,useData{curFOV,dd});
            end
        end

        % take mouse mean
        dataMean = cellfun(@(x) mean(x,1,'omitnan'),dataMouse,'UniformOutput',false);

        % concatenate across mice
        dataCat = zeros(nFOV/fovPerMouse,nDays,120);
        for ff = 1:nFOV/fovPerMouse
            for dd = 1:nDays
                dataCat(ff,dd,:) = dataMean{ff,dd};
            end
        end

        % loop across days
        for dd = 1:nDays
            % get current data
            curDataCat = squeeze(mean(dataCat(:,dd,:),2,'omitnan'));

            % calculate correlations
            corrAll = corr(curDataCat',curDataCat');
            corrAll(eye(size(corrAll,1))==1) = NaN;

            % store data by mouse
            for mm = 1:nFOV/2
                dataInterMouseMean(mm,dd) = mean(corrAll(mm,:),2,'omitnan');
            end
        end

        % re-expand correlation to FOV
        interMouse_byFOV = repelem(dataInterMouseMean,2,1);

        % store correlation values
        samplesAll.(sampleTypes{ss}).interMouseCorr(:,:,rr) = interMouse_byFOV;

    end

    % convert to cell
    samplesAll.(sampleTypes{ss}).interMouseCorr...
        = num2cell(samplesAll.(sampleTypes{ss}).interMouseCorr);
end


%% Load additional activity variables

% define by cell activity structure
cellData = struct();

% load inter-lap data
load('data/intraDayData.mat')
cellData.interLapCorr = cellfun(@(x) x',dataActIntra.full.RBR_sigOthers,'UniformOutput',false);

% load by cell properties
load('data/dataByCell.mat','dataAll')
cellData.spatialSelectivity = dataAll.spatialSelectivity;
cellData.speedScoreAbs = cellfun(@abs,dataAll.speedScore,'UniformOutput',false);

% load dfof
load('data/dfofData.mat','dataDfof')
cellData.dfofSig = cellfun(@(x) x',dataDfof,'UniformOutput',false);

% load decoding properties
load('data/modelGrid_nGrid/decodeInfoGlobal_newDecodeIO.mat','decodingActivity')

% define no process parameters
noSubsample = {'decodeIORaw'};
preSubsample = {'dfofIO','interMouseCorr'};

% load genotypes {[WT],[AD]}
load('groupIDs.mat')
nGroups = length(groups);

activityData = struct();


%% Run collation for all cell types

for ss = 1:nSTypes
    % define current activity
    cellData.decodeIORaw = decodingActivity.decodeIORawGlob.(cellTypes{ss});
    nRCur = samplesAll.(sampleTypes{ss}).nSamples;

    cellData.dfofIO = samplesAll.(sampleTypes{ss}).dfofIO;
    cellData.interMouseCorr = samplesAll.(sampleTypes{ss}).interMouseCorr;

    % define cell data fields
    cellDataFields = fieldnames(cellData);
    nCellFields = length(cellDataFields);

    activityData.(cellTypes{ss}) = struct();


    %% Combine all data

    % cycle through FOV
    for ii = 1:nFOV
        % print progress
        progress = ['Progress: ' num2str(ss) '_' num2str(ii)];
        fprintf(repmat('\b', 1,20));  % backspace to clear message
        fprintf('%-20s',progress);  % pad to fixed width (e.g., 20 chars)

        % get genotype
        for jj = 1:nGroups
            if ismember(ii,groups{jj})
                activityData.(cellTypes{ss})(ii).genotype = groupIDs{jj};
            end
        end

        % get mouse
        curDir = split(foldersLearning{ii}{1},'\');
        curMouse = curDir{5};
        activityData.(cellTypes{ss})(ii).mouse = curMouse;
        activityData.(cellTypes{ss})(ii).day = cell(1,nDays);
        for ff = 1:nCellFields
            activityData.(cellTypes{ss})(ii).(cellDataFields{ff}) = cell(nRCur,nDays);
        end

        % cycle through days
        for jj = 1:nDays
            jTrue = trueDays(ii,jj);
            if jTrue>nDays
                continue
            end

            % get day
            curDir = split(foldersLearning{ii}{jj},'\');
            curDay = curDir{6};
            activityData.(cellTypes{ss})(ii).day{jTrue} = curDay;

            % loop through subsamples
            for rr = 1:nRCur
                % identify current cells
                curCells = samplesAll.(sampleTypes{ss}).cellsSmp{ii,jTrue,rr};

                % store by cell activity
                for ff = 1:nCellFields
                    % calculate current data
                    if ismember(cellDataFields{ff},noSubsample)
                        curCellData = cellData.(cellDataFields{ff}){ii,jTrue};
                    elseif ismember(cellDataFields{ff},preSubsample)
                        curCellData = cellData.(cellDataFields{ff}){ii,jTrue,rr};
                    else
                        curCellData = cellData.(cellDataFields{ff}){ii,jTrue}(curCells);
                    end

                    % store current activity data
                    activityData.(cellTypes{ss})(ii).(cellDataFields{ff}){rr,jTrue} = curCellData;
                end
            end
        end
    end
end


%% Save data

save('/MATLAB Drive/FY2025/Behavior/data/model_vB/activityData_DS_vB.mat','activityData')

