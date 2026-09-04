%% Set plotting parameters
% To skip recalculation, start running from here. All previously calculated
% data is loaded in allMiceData*.mat. This file name must match the current
% save output above.
%

% reclear all variables
clear; close all; clc

p1 = 'D:\AD_Project\Revision\Behavior\Results';
cd(p1)

% load processed behavioral data
load('allMiceDataB5.mat','B','rewCut')

% load genotypes {[WT],[AD]}
load('D:\AD_Project\Behavior\groupIDs.mat')

% define plot colors
colors = {[0 0 1],[1 0 0]};
colors4 = {[102 105 255]/255,[255 102 102]/255;[0 0 175]/255,[175 0 0]/255};


%% Plot success rate across windows sizes

% define fields, titles, and labels
fields = 'perSuccess';
ttls = 'Percent successful runs';
ylabs = 'Success Rate (%)';

% define sexes to plot
useSexes = 1:3;
nSexes = length(useSexes);

% define useCuts
useCuts = 1:5;
nCuts = length(useCuts);

% define sessions
Xname = ['FE' string(1:10)];
useSess = 1:11;
nSess = length(useSess);

% define statistics groups
dayCats = {1,2:11};

% plot lines for individual mice
indON = 0;

% adjust for copying
plotNice = 0;

% extract all data
curData = zeros(length(B),nSess,length(rewCut));
for ii = 1:length(B)
    curD = B(ii).(fields)([B(ii).type{1}; B(ii).type{2}]);
    curData(ii,:,:) = cat(2,curD{:})';
end

% define figure
figure; tiledlayout(nSexes,nCuts)
sgtitle(ttls)

% initialize final statistics
outStats = {};

% loop through sexes
for ss = 1:nSexes
    % current sexes
    curSexes = sexes{useSexes(ss)};
    curSexID = sexIDs{useSexes(ss)};

    for cc = 1:nCuts
        % current cut
        curCutIdx = useCuts(cc);
        curCutVal = rewCut(curCutIdx);

        % initialize data matrices
        data = cell(1,2);

        % sort data by group
        for g = 1:2
            data{g} = curData(intersect(groups{g},curSexes),useSess,curCutIdx);
        end

        % calculate p values
        pMC = nan(nSess,1);
        pAnova = nan(nSess,2);

        for ii = 1:length(dayCats)
            % find intersect of plotting x values with day category
            curCat = dayCats{ii};
            lenCur = length(curCat);
            if lenCur>1
                % calculate anova p values with multiple comparisons
                [pACur,pMCCur] = anovaRM2W_full_BH...
                    (data{1}(:,curCat),data{2}(:,curCat),1);
                pAnova(curCat,:) = pACur([1 3])'.*ones(lenCur,2);
            else
                % calculat single day p values
                [~,pMCCur] = ttest2(data{1}(:,curCat),data{2}(:,curCat));
            end

            pMC(curCat) = pMCCur;
        end

        % initialize figure
        nexttile(cc+(ss-1)*nCuts); hold on

        % plot line graph
        plotErrorSig(useSess,data{1},data{2},groupIDs,pAnova,pMC,colors,indON)

        % set figure labels
        legend('off')
        xlabel('Session')
        ylabel(ylabs)
        title([curSexID ': ' num2str(curCutVal) 'cm cut' ])
        set(gca,'FontSize',12)

        % set x axis and labels
        set(gca,'XTick',useSess)
        set(gca,'XTickLabels',Xname(useSess))
        xlim([useSess(1)-0.5 useSess(end)+0.5])

        % adjust for copying
        if plotNice
            axis('square')
        end

        % perform final statistics
        nUnits = 'mice';
        testCat = [curSexID ': ' num2str(curCutVal) 'cm cut' ];
        outStats(end+1,:) = [testCat anovaEffectSizeMinimal(data{1}(:,2:end),data{2}(:,2:end),nUnits)];
    end
end

% save figure
savefig('success_cutZone.fig')

% save figure as tif if indivudual plot is on to preserve transparency
if indON
    exportgraphics(gcf,'success_cutZone.tif')
end


%% Plot stop distributions

% distribution types
distTypes = {'All','Pre','Post'};
distTypesIdx = {2:11,2,8:11};
nDTypes = length(distTypes);

% define labels
fields = 'stopDist';
ttls = 'Stop distribution';
ylabs = 'Stop distribution';

% define sexes to plot
useSexes = 1:3;

% close all figures
clAll = 0;

% define x-axis to plot
xPlot = 0:5:600;

% xZones = {(20:5:475)/5,(480:5:560)/5};
xZones = {(20:5:480)/5,(510:5:560)/5};

% extract data for current field
curDataCat = cell(length(B),1);
for ii = 1:length(B)
    curData = B(ii).(fields)([B(ii).type{1}; B(ii).type{2}]);
    curDataCat{ii} = cat(1,curData{:});
end

% initialize figure
figure; tiledlayout(nSexes,nDTypes)
sgtitle(ttls)

% loop through sexes
for ss = 1:nSexes
    % current sexes
    curSexes = sexes{useSexes(ss)};
    curSexID = sexIDs{useSexes(ss)};

    for dd = 1:nDTypes
        % calculate mean distribution
        curDataMean = cellfun(@(x) mean(x(distTypesIdx{dd},:),1,'omitnan'),...
            curDataCat,'UniformOutput',false);

        % sort data by group
        data = cell(1,2);
        for g = 1:2
            data{g} = cat(1,curDataMean{intersect(groups{g},curSexes)});
        end

        % calculate p values
        [pACur,pMC] = anovaRM2W_full_BH(data{1},data{2},0);
        pAnova = pACur([1 3])'.*ones(length(pMC),2);

        % initialize figure
        nexttile(dd+(ss-1)*nDTypes); hold on

        % plot line graph
        plotErrorSig(xPlot,data{1},data{2},groupIDs,pAnova,pMC,colors)

        % set figure labels
        legend('off')
        xlabel('Track Position (cm)')
        ylabel(ylabs)
        title([distTypes{dd} ' Learning: ' curSexID])
        set(gca,'FontSize',12)

        % plot reward zone and false zone for distributions
        yPK = max(ylim);
        fill([510 510 560 560],[0 yPK yPK 0],[1 1 0.7],'LineStyle','None',...
            'FaceAlpha',0.5,'DisplayName','Reward Zone')

        % set x axis and labels
        xlim([0 600])
    end
end

% save figure
savefig('stopDistribution.fig')

% close figure
if clAll
    close
end


%% Plot rewarded stop distributions

% distribution types
distTypes = {'All','Post'};
distTypesIdx = {2:11,8:11};
nDTypes = length(distTypes);

% define labels
fields = 'yStopRew';
ttls = 'Reward Stop Distribution';
ylabs = 'Reward stop distribution';

% define sexes to plot
useSexes = 1:3;

% close all figures
clAll = 0;

% define x-axis to plot
xPlot = 510:5:560;
xLim = [500 570];

% xZones = {(20:5:475)/5,(480:5:560)/5};
xZones = {(20:5:480)/5,(510:5:560)/5};

% extract data for current field
curData = cell(length(B),1);
for ii = 1:length(B)
    tempData = B(ii).(fields)([B(ii).type{1}; B(ii).type{2}]);
    nullData = cellfun(@(a) cellfun(@(b) isempty(b),a),tempData,'UniformOutput',false);
    curData{ii} = cellfun(@(x,y) cat(1,x{~y}),tempData,nullData,'UniformOutput',false);
end

% initialize figure
figure; tiledlayout(nSexes,nDTypes)
sgtitle(ttls)

% loop through sexes
for ss = 1:nSexes
    % current sexes
    curSexes = sexes{useSexes(ss)};
    curSexID = sexIDs{useSexes(ss)};

    for dd = 1:nDTypes
        % calculate per mouse distribution
        curDataMouse = cellfun(@(x) cat(1,x{distTypesIdx{dd}}),curData,'UniformOutput',false);
        curDataDist = cellfun(@(x) ksdensity(x,xPlot,'Function','pdf'),curDataMouse,'UniformOutput',false);
        % curDataDist = cellfun(@(x) ksdensity(x,xPlot,'Function','cdf'),curDataMouse,'UniformOutput',false);

        % sort data by group
        data = cell(1,2);
        for g = 1:2
            data{g} = cat(1,curDataDist{intersect(groups{g},curSexes)});
        end

        % calculate p values
        [pACur,pMC] = anovaRM2W_full_BH(data{1},data{2},0);
        pAnova = pACur([1 3])'.*ones(length(pMC),2);

        % initialize figure
        nexttile(dd+(ss-1)*nDTypes); hold on

        % plot line graph
        plotErrorSig(xPlot,data{1},data{2},groupIDs,pAnova,pMC,colors)

        % set figure labels
        legend('off')
        xlabel('Track Position (cm)')
        ylabel(ylabs)
        title([distTypes{dd} ' Learning: ' curSexID])
        set(gca,'FontSize',12)
        axis('square')

        % set x axis and labels
        xlim(xLim)

        % perform final statistics
        nUnits = 'mice';
        testCat = [distTypes{dd} ' Learning: ' curSexID];
        outStats(end+1,:) = [testCat anovaEffectSizeMinimal(data{1},data{2},nUnits)];
    end
end

% save figure
savefig('rewardStopDistribution.fig')

% close figure
if clAll
    close
end