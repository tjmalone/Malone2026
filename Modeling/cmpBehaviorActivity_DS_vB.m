%% Compare behavior and activity

clear; close all; clc

cd('/MATLAB Drive/FY2025/Behavior')


%% Load data

% load behavior
load('data/allMiceDataB5.mat')

% initialize behavior variables
nBehavior = length(B);
behMouse = cell(nBehavior,1);
behDay = cell(nBehavior,1);
perSuccess = cell(nBehavior,1);
dprimeGlobal = cell(nBehavior,1);
meanVel = cell(nBehavior,1);

% extract behavior
for ii = 1:nBehavior
    behMouse{ii} = B(ii).mouse{1};
    behDay{ii} = B(ii).day;

    perSuccess{ii} = B(ii).perSuccess;
    dprimeGlobal{ii} = B(ii).dprimeGlobal;
    meanVel{ii} = B(ii).avgVelMovingRaw;
end

% load activity
load('data/model_vB/activityData_DS_vB.mat','activityData')
cellTypes = fieldnames(activityData);


%% Loop through cell types

for cc = 1:length(cellTypes)
    %% Load current activity data

    % define current activity
    activityDataUse = activityData.(cellTypes{cc});

    % initialize activity variables
    nActivity = length(activityDataUse);
    actMouse = cell(nActivity,1);
    actDay = cell(nActivity,1);

    % initialize variable size actiivty data
    actFields = fieldnames(activityDataUse);
    actFields = actFields(4:end);
    nActFields = length(actFields);
    nRnd = size(activityDataUse(1).(actFields{1}),1);

    % extract mouse info
    for ii = 1:nActivity
        actMouse{ii} = erase(activityDataUse(ii).mouse,'ID_');
        actDay{ii} = activityDataUse(ii).day;
    end


    %% Parse data by day/mouse/morphology

    % initialize activity data
    dataAct = cell(0,4);

    % cycle through activity FOV
    for ii = 1:nActivity
        % extract mouse information
        curGeno = activityDataUse(ii).genotype;
        curMouse = actMouse{ii};

        % loop through activity days
        for jj = 1:length(actDay{ii})
            % extract day and activity
            cDay = actDay{ii}{jj};

            cAct = cell(nRnd,nActFields);
            % loop through rnd subsamples
            for rr = 1:nRnd
                for af = 1:nActFields
                    cAct{rr,af} = activityDataUse(ii).(actFields{af}){rr,jj};
                end
            end

            % generate data matrix
            dataAct(end+1,1:4) = {curGeno,curMouse,cDay,cAct};
        end
    end

    % initialize behavior data
    dataBeh = cell(0,3);

    % cycle through activity FOV
    for ii = 1:nBehavior
        % extract mouse information
        curMouse = behMouse{ii};

        % cycle through activity days
        for jj = 1:length(behDay{ii})
            % extract day and activity
            cDay = behDay{ii}{jj};
            cPerSuccess = perSuccess{ii}(jj);

            % generate data matrix
            dataBeh(end+1,1:3) = {curMouse,cDay,cPerSuccess};
        end
    end


    %% Combine data by day

    nBeh = size(dataBeh,1);
    % Removed detail titles due to large number of variables
    corrLabels = ['perSuccess'; actFields];

    corrData = nan(nBehavior,length(actDay{1}),1+nActFields,nRnd);

    for ii = 1:nBeh
        % get behavior mouse/day
        curMouse = dataBeh{ii,1};
        mouseIdx = find(strcmp(curMouse,behMouse),1);
        curDay = dataBeh{ii,2};

        % match to activity mouse/day
        mMouse = strcmp(dataAct(:,2),curMouse);
        mDay = strcmp(dataAct(:,3),curDay);
        dayIdx = find(mDay(mMouse==1));

        % get matches
        actMatch = find(mMouse & mDay);

        % skip non matched behavior days
        if isempty(actMatch); continue; end

        if length(actMatch)==4
            if all(isnan(corrData(mouseIdx,dayIdx(1),:)))
                dayIdx = dayIdx(1);
                actMatch = actMatch([1 3]);
            else
                dayIdx = dayIdx(2);
                actMatch = actMatch([2 4]);
            end
        else
            dayIdx = dayIdx(1);
        end

        % calculate final activity correlation data
        mAct = cat(3,dataAct{actMatch,4});
        mMeanAct = zeros(nActFields,nRnd);
        for rr = 1:nRnd
            for af = 1:nActFields
                mActCur = cat(1,mAct{rr,af,:});
                mMeanAct(af,rr) = mean(mActCur,'omitnan');
            end
        end

        % calculate final behavior/mouse correlation data
        mPerSuccess = repmat(dataBeh{ii,3},1,nRnd);
        mGenotype = dataAct{actMatch(1),1};

        % save data
        if ~all(isnan(corrData(mouseIdx,dayIdx,:)))
            error('Overlap')
        end
        corrData(mouseIdx,dayIdx,:,:) = [mPerSuccess;mMeanAct];
    end

    save(['data/model_vB/corrDataDS_' cellTypes{cc} '_vB.mat'],'corrData','corrLabels')
end

