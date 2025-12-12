%% collateActivityData_new
% Collates activity data for correlation with behavior

clear; close all; clc

baseFolder = '/MATLAB Drive/FY2025/imagingData';
cd(baseFolder)

% load data
load('foldersLearning.mat')

useDays = 1:11;
nDays = length(useDays);
nFOV = length(foldersLearning);

% define by cell activity structure
cellData = struct();

% load inter-day consistency data
load('data/interDaySubset.mat','dataOut')
cellData.interDayCorr = cell(nFOV,nDays);
for ff = 1:nFOV
    for dd = 1:nDays
        if ismember(dd,3:nDays-1)
            % cellData.interDayCorr{ff,dd} = mean([dataOut{ff,dd},dataOut{ff,dd-1}],2);
            cellData.interDayCorr{ff,dd} = mean([dataOut{ff,dd},dataOut{ff,dd-1}],2,'omitnan');
        elseif dd==2
            cellData.interDayCorr{ff,dd} = dataOut{ff,dd};
        elseif dd==nDays
            cellData.interDayCorr{ff,dd} = dataOut{ff,nDays-1};
        else
            cellData.interDayCorr{ff,dd} = nan(size(dataOut{ff,1}));
        end
    end
end

% load inter-mouse consistency
load('data/interMouseData.mat','dataInterMouseMean')
cellData.interMouseCorr = num2cell(repelem(dataInterMouseMean,2,1));

% load inter-lap data
load('data/intraDayData.mat')
cellData.interLapCorr = cellfun(@(x) x',dataActIntra.full.RBR_sigOthers,'UniformOutput',false);

% load by cell properties
load('data/dataByCell.mat','dataAll')
cellData.corrNear = dataAll.corrNear;
cellData.corrFar = dataAll.corrFar;
cellData.spatialSelectivity = dataAll.spatialSelectivity;
cellData.FldNum = dataAll.FWidthNum;
cellData.FldWidth = dataAll.FWidthMean;
cellData.FldCoverage = dataAll.FWidthSum;
cellData.spatialInfo = dataAll.spatialInfo;
cellData.speedScoreRaw = dataAll.speedScore;
cellData.speedScoreAbs = cellfun(@abs,dataAll.speedScore,'UniformOutput',false);

% process self correlation data
cellData.selfRunCorr = cellfun(@(x,y) x-y,dataAll.selfRun_S,dataAll.selfRun_F,'UniformOutput',false);

% load dfof
load('data/dfofData.mat','dataDfof')
cellData.dfofSig = cellfun(@(x) x',dataDfof,'UniformOutput',false);

% load decoding accuracy, Updated
useItr = 21;
load('data/decoding/decodingDataMaster_all.mat','decodingDataItr')
cellData.decoding = num2cell(decodingDataItr(useItr).meanPerCorrect);

% load decoding properties
load('data/decodeInfoGlobal_newDecodeIO.mat','decodingActivity')
cellData.decodeIORaw = decodingActivity.decodeIORawGlob.allType;
cellData.decodeIONorm = decodingActivity.decodeIONormGlob.allType;

% load amplitude difference
load('data/cueIdentitySubset.mat','dataOut')
cellData.ampDiffPerNL = dataOut.cuePer;
cellData.ampDiffIORaw = dataOut.IOraw;
cellData.ampDiffIONorm = dataOut.IOnorm;

% load cue anchoring info (ridge-background)
load('data/corrInfoGlobal.mat','dataActivity')
dfofIO = dataActivity.dfofIO.common.allSex.allMorph;

% define cell data fields
cellDataFields = fieldnames(cellData);
nCellFields = length(cellDataFields);

% define no process parameters
noProcess = {'interDayCorr','interMouseCorr','decoding','decodeIORaw','decodeIONorm',...
    'ampDiffPerNL','ampDiffIORaw','ampDiffIONorm'};

% load genotypes {[WT],[AD]}
load('groupIDs.mat')
nGroups = length(groups);

% combine disjoint inputs
dfofIOFull = cell(1,nGroups);
for gg = 1:nGroups
    curData = dfofIO.(groupIDs{gg});
    curEmpty = cellfun(@isempty,curData);
    curData(curEmpty) = {NaN};
    dfofIOFull{gg} = cell2mat(curData);
end
dfofIOCat = max(dfofIOFull{1},dfofIOFull{2});
dfofIORaw = num2cell(dfofIOCat);
dfofIOAbs = num2cell(abs(dfofIOCat-1));


%% Combine all data

activityData = struct();

% cycle through FOV
for ii = 1:nFOV

    % get genotype
    for jj = 1:nGroups
        if ismember(ii,groups{jj})
            activityData(ii).genotype = groupIDs{jj};
        end
    end

    % get mouse
    curDir = split(foldersLearning{ii}{1},'\');
    curMouse = curDir{5};
    activityData(ii).mouse = curMouse;
    activityData(ii).day = cell(1,nDays);
    for ff = 1:nCellFields
        activityData(ii).(cellDataFields{ff}) = cell(1,nDays);
    end

    % cycle through days
    for jj = 1:nDays
        if trueDays(ii,jj)>nDays
            continue
        end

        % get day
        curDir = split(foldersLearning{ii}{jj},'\');
        curDay = curDir{6};
        activityData(ii).day{trueDays(ii,jj)} = curDay;

        % store by cell activity
        for ff = 1:nCellFields
            if ~ismember(cellDataFields{ff},noProcess)
                curCellData = cellData.(cellDataFields{ff}){ii,trueDays(ii,jj)}(alignsLearning{ii}(:,jj));
                activityData(ii).(cellDataFields{ff}){trueDays(ii,jj)} = curCellData;
            else
                activityData(ii).(cellDataFields{ff}){trueDays(ii,jj)} = cellData.(cellDataFields{ff}){ii,trueDays(ii,jj)};
            end
        end

    end

    % store by FOV activity
    activityData(ii).dfofIORaw = dfofIORaw(ii,:);
    activityData(ii).dfofIOAbs = dfofIOAbs(ii,:);
   
    cd(baseFolder)
end


%% Save data

save('/MATLAB Drive/FY2025/Behavior/data/model_vB/activityDataVB.mat','activityData')
