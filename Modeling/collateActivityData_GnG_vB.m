%% collateActivityData_new
% Collates activity data for correlation with behavior

clear; close all; clc

baseFolder = '/MATLAB Drive/FY2025/imagingData';
cd(baseFolder)

% load data
load('foldersLearning.mat')

% load cell information
load('data/cellSelect.mat','cellSelect')

useDays = 1:11;
nDays = length(useDays);
nFOV = length(foldersLearning);
cellTypes = {'grid','nongrid'};

% define by cell activity structure
cellData = struct();

% load inter-mouse consistency
load('data/interMouseData_grid.mat','dataInterMouseMean')
interMouseCorr.grid = num2cell(repelem(dataInterMouseMean,2,1));
load('data/interMouseData_nongrid.mat','dataInterMouseMean')
interMouseCorr.nongrid = num2cell(repelem(dataInterMouseMean,2,1));

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

% load cue anchoring info (dfof IO)
load('data/corrInfoGlobal.mat','dataActivity')

% define no process parameters
noProcess = {'interMouseCorr','decodeIORaw','dfofIORaw'};

% load genotypes {[WT],[AD]}
load('groupIDs.mat')
nGroups = length(groups);

activityData = struct();


%% Run collation for all cell types

for cc = 1:length(cellTypes)
    % define current cells
    cellSelectUse = cellSelect.learn.allSex.(cellTypes{cc}).allMorph;
    cellSelectCat = combCells(cellSelectUse.WT,cellSelectUse.AD);

    % define current activity
    cellData.decodeIORaw = decodingActivity.decodeIORawGlob.(cellTypes{cc});
    dfofIORaw = dataActivity.dfofIO.(cellTypes{cc}).allSex.allMorph;
    cellData.dfofIORaw = combCells(dfofIORaw.WT,dfofIORaw.AD);
    cellData.interMouseCorr = interMouseCorr.(cellTypes{cc});

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

save('/MATLAB Drive/FY2025/Behavior/data/model_vB/activityData_GnG_vB.mat','activityData')

