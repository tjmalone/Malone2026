%% collateActivityData_new
% Collates activity data for correlation with behavior

clear; close all; clc

baseFolder = 'D:\AD_Project/imagingData';
cd(baseFolder)

% load data
load('foldersLearning.mat')

% load cell information
load('data/cellSelect.mat','cellSelect')

useDays = 1:11;
nDays = length(useDays);
nFOV = length(foldersLearning);
cellTypes = {'ste','pyr'};

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

% load cue anchoring info (dfof IO)
load('data/corrInfoGlobal.mat','dataActivity')

% define no process parameters
noProcess = {'dfofIORaw'};

% load genotypes {[WT],[AD]}
load('groupIDs.mat')
nGroups = length(groups);

activityData = struct();


%% Run collation for all cell types

for cc = 1:length(cellTypes)
    % define current cells
    cellSelectUse = cellSelect.learn.allSex.common.(cellTypes{cc});
    cellSelectCat = combCells(cellSelectUse.WT,cellSelectUse.AD);

    % define current activity
    dfofIORaw = dataActivity.dfofIO.common.allSex.(cellTypes{cc});
    cellData.dfofIORaw = combCells(dfofIORaw.WT,dfofIORaw.AD);

    % define cell data fields
    cellDataFields = fieldnames(cellData);
    nCellFields = length(cellDataFields);

    activityData.(cellTypes{cc}) = struct();


    %% Combine all data

    % cycle through FOV
    for ii = 1:nFOV

        % get genotype
        for jj = 1:nGroups
            if ismember(ii,groups{jj})
                activityData.(cellTypes{cc})(ii).genotype = groupIDs{jj};
            end
        end

        % get mouse
        curDir = split(foldersLearning{ii}{1},'\');
        curMouse = curDir{5};
        activityData.(cellTypes{cc})(ii).mouse = curMouse;
        activityData.(cellTypes{cc})(ii).day = cell(1,nDays);
        for ff = 1:nCellFields
            activityData.(cellTypes{cc})(ii).(cellDataFields{ff}) = cell(1,nDays);
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
            activityData.(cellTypes{cc})(ii).day{jTrue} = curDay;

            % store by cell activity
            for ff = 1:nCellFields
                if ~ismember(cellDataFields{ff},noProcess)
                    curCellData = cellData.(cellDataFields{ff}){ii,jTrue}(cellSelectCat{ii,jTrue});
                    activityData.(cellTypes{cc})(ii).(cellDataFields{ff}){jTrue} = curCellData;
                else
                    activityData.(cellTypes{cc})(ii).(cellDataFields{ff}){jTrue} =...
                        cellData.(cellDataFields{ff}){ii,jTrue};
                end
            end

        end

        cd(baseFolder)
    end
end


%% Save data

save('D:\AD_Project/Behavior/data/model_vB/activityData_PS_vB.mat','activityData')

