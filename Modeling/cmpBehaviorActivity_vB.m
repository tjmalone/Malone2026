%% Compare behavior and activity

clear; close all; clc

cd('D:\AD_Project\Behavior')


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
load('data/model_vB/activityDataVB.mat','activityData')

% initialize activity variables
nActivity = length(activityData);
actMouse = cell(nActivity,1);
actDay = cell(nActivity,1);

% initialize variable size actiivty data
actFields = fieldnames(activityData);
actFields = actFields(4:end);
nActFields = length(actFields);

% extract mouse info
for ii = 1:nActivity
    actMouse{ii} = erase(activityData(ii).mouse,'ID_');
    actDay{ii} = activityData(ii).day;
end


%% Parse data by day/mouse/morphology

% initialize activity data
dataAct = cell(0,4);

% cycle through activity FOV
for ii = 1:nActivity
    % extract mouse information
    curGeno = activityData(ii).genotype;
    curMouse = actMouse{ii};

    % cycle through activity days
    for jj = 1:length(actDay{ii})
        % extract day and activity
        cDay = actDay{ii}{jj};

        cAct = cell(1,nActFields);
        for af = 1:nActFields
            cAct{af} = activityData(ii).(actFields{af}){jj};
        end

        % generate data matrix
        dataAct(end+1,1:4) = {curGeno,curMouse,cDay,cAct};
    end
end

% initialize behavior data
dataBeh = cell(0,5);

% cycle through activity FOV
for ii = 1:nBehavior
    % extract mouse information
    curMouse = behMouse{ii};

    % cycle through activity days
    for jj = 1:length(behDay{ii})
        % extract day and activity
        cDay = behDay{ii}{jj};
        cPerSuccess = perSuccess{ii}(jj);
        cDPrimeGlobal = dprimeGlobal{ii}(jj);
        cMeanVel = meanVel{ii}(jj);

        % generate data matrix
        dataBeh(end+1,1:5) = {curMouse,cDay,cPerSuccess,cDPrimeGlobal,cMeanVel};
    end
end


%% Combine data by day

nBeh = size(dataBeh,1);
% Removed detail titles due to large number of variables
% corrLabels = {'Percent Success','Global d''','Inter-lap Consistency',...
%     'Spatial Selectivity','# of Fields','Success/Fail Self-Correlation',...
%     'Ridge-Background Ratio (Raw)','Ridge-Background Ratio (Abs. anchor)',...
%     'Success/Fail Self-Correlation/Ridge-Background Ratio',...
%     'Spatial Selectivity/Ridge-Background Ratio','Full Combination'};
corrLabels = ['perSuccess'; 'dPrimeGlobal'; 'meanVel'; actFields];

corrDataBase = nan(nBehavior,length(actDay{1}),3+nActFields);

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
        if all(isnan(corrDataBase(mouseIdx,dayIdx(1),:)))
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
    mAct = cat(1,dataAct{actMatch,4});
    mMeanAct = zeros(1,nActFields);
    for af = 1:nActFields
        if iscell(mAct)
            mActCur = cat(1,mAct{:,af});
        else
            mActCur = mAct;
        end

        mMeanAct(af) = mean(mActCur,'omitnan');
    end

    % calculate final behavior/mouse correlation data
    mPerSuccess = dataBeh{ii,3};
    mDprime = dataBeh{ii,4};
    mMeanVel = dataBeh{ii,5};
    mGenotype = dataAct{actMatch(1),1};

    % save data
    if ~all(isnan(corrDataBase(mouseIdx,dayIdx,:)))
        error('Overlap')
    end
    corrDataBase(mouseIdx,dayIdx,:) = [mPerSuccess,mDprime,mMeanVel,mMeanAct];
end

% normalize behavior variables
corrNorm = zeros(size(corrDataBase));
for ii = 1:size(corrDataBase,3)
    % curArray = corrDataBase(:,:,ii);
    % minVal = min(curArray,[],'all');
    % maxVal = max(curArray,[],'all');
    % corrNorm(:,:,ii) = (curArray-minVal)/(maxVal-minVal);

    curArray = reshape(corrDataBase(:,:,ii),[],1);
    curArrayNorm = normalize(curArray);
    corrNorm(:,:,ii) = reshape(curArrayNorm,size(corrDataBase(:,:,ii)));
end

% combine variables
combVars = {};
corrComb = zeros([size(corrDataBase,[1 2]) length(combVars)]);
labComb = cell(length(combVars),1);

for ii = 1:length(combVars)
    curComb = corrNorm(:,:,abs(combVars{ii}));

    curCombSigned = zeros(size(curComb));
    for jj = 1:length(combVars{ii})
        if combVars{ii}(jj)>0
            curCombSigned(:,:,jj) = curComb(:,:,jj);
        else
            % curCombSigned(:,:,jj) = 1-curComb(:,:,jj);
            curCombSigned(:,:,jj) = -curComb(:,:,jj);
        end
    end

    % add combined activity
    corrComb(:,:,ii) = mean(curCombSigned,3);

    % update labels
    labComb{ii} = strjoin(actFields(abs(combVars{ii})-2),'-');

end

corrData = cat(3,corrDataBase,corrComb);
corrLabels = [corrLabels;labComb];

save('data/model_vB/corrData_vB.mat','corrData','corrLabels')


%% Initialize figure parameters

clear; close all; clc

p1 = 'D:\AD_Project\Behavior';
cd(p1)

% load data
load('data/model_vB/corrData_vB.mat')

nFOV = size(corrData,2);
nDays = size(corrData,2);

% load genotypes {[WT],[AD]}
load('groupIDs.mat')
nGroups = length(groups);
nSexes = length(sexes);
sexIDs = {'allSex','female','male'};

% define legend
legs = [groupIDs 'All'];

% define plot colors
colors = {[0 0 1],[1 0 0]};
colors6 = {[0 0 1],[1 0 0];[102 105 255]/255,[255 102 102]/255;[0 0 175]/255,[175 0 0]/255};

% set save folder name
svFile = [p1 '/Figures/Figures' char(datetime('today','format','yyyyMMdd')) '/activityCorr_vB'];
if ~isfile(svFile)
    mkdir(svFile);
    mkdir([svFile '/data']);
end


%% Plot correlation pairs

close all

% set correlation parameters
params = [];
useB = [1];
useA = [2:25];  % all behavior and combined variables
nB = length(useB);
nA = length(useA);
params(1,:) = repmat(useB,1,nA);
params(2,:) = repelem(useA,1,nB);

% close all figures
clAll = 1;
% curmax = 0;
yMax = 4.5;

% whether to average across days
avgON = [1,2];
avgLabels = {' (by session)',' (by mouse)'};
avgSaves = {'session','mouse'};

% day selections
useDays = {2:11};
% useDays = {2:11,2:11};

corrType = 'Pearson';
% corrType = 'Kendall';
% corrType = 'Spearman';

% initialize statistics
outStats = {};
testNameCorr = 'two-tailed Pearson linear correlation';

% loop through correlation pairs
for ff = 1:size(params,2)
    % loop through average types
    for gg = 1:length(useDays)
        %% Generate single figure for all subtypes

        % initialize figure
        figure
        xLab = corrLabels{params(1,ff)};
        yLab = corrLabels{params(2,ff)};
        t = tiledlayout(1,nSexes);
        sgtitle(['Correlation: ' xLab ' vs. ' yLab avgLabels{avgON(gg)}])

        % initialize save data
        dataSave = cell(1,nSexes);

        for ss = 1:nSexes
            % define current data
            data = corrData(:,useDays{gg},params(:,ff));
            useDaysCur = 1:length(useDays{gg});

            % normalize data
            nData = data(:,:,2);
            nData = normalize(nData(:));
            data(:,:,2) = reshape(nData(:),size(data(:,:,2)));

            % curmax = max(curmax,max(abs(data(:,:,2)),[],'all'));

            % define final use data
            if avgON(gg)==1
                useData = data;
            else
                useData = mean(data(:,useDaysCur,:),2,'omitnan');
            end

            % separate data by genotype
            dataFinal = cell(1,nGroups);
            for ii = 1:nGroups
                curGroup = intersect(groups{ii},sexes{ss});
                df = useData(curGroup,:,:);
                curData = squeeze(reshape(df,[],1,2));
                curData(any(isnan(curData),2),:) = [];
                dataFinal{ii} = curData;
            end

            % calculate linear correlation
            dataComb = cat(1,dataFinal{:});
            [rCorr,pCorr] = corr(dataComb(:,1),dataComb(:,2),'type',corrType);

            % plot scatter
            for ii = 1:nGroups
                xCur = dataFinal{ii}(:,1);
                yCur = dataFinal{ii}(:,2);

                if ss~=1
                    nexttile(ss); hold on
                    scatter(xCur,yCur,25,colors6{ss,ii},'LineWidth',0.75,'MarkerEdgeAlpha',0.5);

                    nexttile(1); hold on
                    scatter(xCur,yCur,25,colors6{ss,ii},'filled','LineWidth',0.75,'MarkerEdgeAlpha',1);
                else
                    nexttile(1); hold on
                    scatter(xCur,yCur,1,'MarkerEdgeAlpha',0);
                end
            end

            nexttile(ss); hold on
            if ss==1
                xLims = xlim;
                yLims = ylim;
                yLims = [-yMax yMax];
                ylim(yLims)
            else
                xlim(xLims)
                ylim(yLims)
            end

            % initialize labels
            labelLeg = cell(1,3);
            h = zeros(1,nGroups);

            % plot trend lines
            for ii = 1:nGroups+1
                if ii<=nGroups
                    xCur = dataFinal{ii}(:,1);
                    yCur = dataFinal{ii}(:,2);
                    curColor = colors6{ss,ii};
                else
                    xCur = dataComb(:,1);
                    yCur = dataComb(:,2);
                    curColor = 'k';
                end

                % calculate trendline
                coefficients = polyfit(xCur,yCur, 1);
                xFit = [min(xlim) max(xlim)];
                yFit = polyval(coefficients, xFit);

                % plot trendline
                h(ii) = plot(xFit,yFit,'-','Color',curColor,'LineWidth',0.5);

                % calculate correlation
                [rCur,pCur] = corr(xCur,yCur,'type',corrType);
                labelLeg{ii} = [legs{ii} ': r=' num2str(rCur,2) ', p=' num2str(pCur,2)];

                % calculate full statistics
                corrCat = [xLab ' vs. ' yLab ': ' sexIDs{ss} '-' legs{ii}];
                outStats(end+1,:) = [corrCat corrEffectSize(xCur,yCur,testNameCorr,'mouse sessions')];
            end

            legend(h,labelLeg,'Location','best','AutoUpdate','off')
            xlabel(xLab)
            ylabel(yLab)
            set(gca,'FontSize',12)
            title([sexIDs{ss} avgLabels{avgON(gg)}...
                ': r=' num2str(rCorr,2) ', p=' num2str(pCorr,2)])
            % axis('square')

            % store data for saving
            dataSave{ss} = dataFinal;
        end

        % save figure
        savefig([svFile '/Corr_' corrLabels{params(1,ff)} '-'...
            corrLabels{params(2,ff)} '_by' avgSaves{avgON(gg)} '.fig'])

        % save data
        save([svFile '/data/Corr_' corrLabels{params(1,ff)} '-'...
            corrLabels{params(2,ff)} '_by' avgSaves{avgON(gg)} '.mat'],'corrLabels','dataSave')

        % close figure
        if clAll
            close
        end
    end
end

