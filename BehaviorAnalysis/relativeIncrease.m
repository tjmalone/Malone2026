%% relativeIncrease.m
% Plots the ratio between WT and AD mouse behavior as a function of time.
% Behavior of both genotypes is normalized to WT behavior. A linear
% correlation of the AD data is used to demonstrate the change in the
% difference between WT and AD mice over time. Currently, individual AD
% values are used for correlation, but can be switched to mean values
%

clear; close all; clc


%% Set input parameters

% set correct directory
cd('D:\AD_Project\Behavior')
p = pwd;

% load mouse data
load('data\allMiceDataB4.mat')

% load genotypes {[WT],[AD]}
load('groupIDs.mat')

% manually set genotypes
% groups = {[1 2 3 5 6 7 12],[4 8 9 10 11 13 14]};

% define legend and plot colors
legs = {'WT','AD'};
colors = {[0 0 1],[1 0 0]};

% define use days
useIdx = 5:15;      % novelty days

% define x label
X = 1:length(useIdx);

% define statistics groups
dayCats = {1:length(useIdx)};
nCats = length(dayCats);

% set save folder name
svFile = [p '\Figures\Figures' char(datetime('today','format','yyyyMMdd'))];
mkdir(svFile);
mkdir([svFile '\data']);

% set save file name
svName = 'relativeIncrease';


%% Plot line graphs by session

% define fields and graph labels
fields = {'perSuccess','perSuccessTry'};
ttls = {'Percent successful ratio','Percent try ratio'};
ylabs = {'Percent','Percent'};
Xname = ['FE' string(1:10)];

% cycle through fields
for f = 1:length(fields)
    %% Process current field
    
    % intialize data
    curData = zeros(length(B),length(useIdx));
    data = cell(1,2);
    
    % extract data for current field
    for ii = 1:length(B)
        curD = B(ii).(fields{f})(useIdx);
        curData(ii,:) = curD;
    end
    
    % sort data by group
    for g = 1:2
        data{g} = curData(groups{g},:);
    end
    
    % normalize data
    normVals = mean(data{1},1);
    dataNorm = cellfun(@(a) a./normVals,data,'UniformOutput',0);
    
    % calculate p values
    pMC = nan(length(useIdx),1);
    pAnova = nan(length(useIdx),2);
    
    for ii = 1:nCats
        lenCur = length(dayCats{ii});
        if lenCur>1
            % calculate anova p values with multiple comparisons
            [pACur,pMCCur] = anovaRM2W_full_BH...
                (dataNorm{1}(:,dayCats{ii}),dataNorm{2}(:,dayCats{ii}),1);
            % extract group difference and group-time interaction
            pAnova(dayCats{ii},:) = pACur([1 3])'.*ones(lenCur,2);
        else
            % calculate single day p values
            [~,pMCCur] = ttest2(dataNorm{1}(:,dayCats{ii}),...
                dataNorm{2}(:,dayCats{ii}));
        end
        
        pMC(dayCats{ii}) = pMCCur;
    end
    
%     % calculate linear correlation for AD mean
%     meanNorm = mean(data{2},1);
%     [rCorr,pCorr] = corr((1:length(meanNorm))',meanNorm');

    % calculate linear correlation for AD data
    meanNorm = mean(data{2},1);
    [rCorr,pCorr] = corr(repmat((1:length(meanNorm))',...
        length(groups{2}),1),reshape(data{2}',[],1));
    
    
    %% Plot data
    
    F = figure; hold on
    
    % plot data
    plotErrorSig(X,dataNorm{1},dataNorm{2},legs,pAnova,pMC,colors,0)
    
    % set graph labels
    legend('Location','eastoutside')
    set(gca,'FontSize',12)
    xlabel('Session')
    xlim([X(1)-0.5 X(end)+0.5])
    ylabel(ylabs{f})
    title(ttls{f})
    ylim([0 1.5])
    set(F,'Position',[100 100 800 400])

    set(gca,'XTick',1:length(Xname))
    set(gca,'XTickLabels',Xname)

    % set title with correlation info
    title([ttls{f} ' r=' num2str(rCorr,2) ', p=' num2str(pCorr,2)])
    
    % save figure
    savefig([svFile '\' svName '_' fields{f} '.fig'])
    close
    
    % save underlying data
    save([svFile '\data\' svName '_' fields{f} '.mat'],'dataNorm')

    
end

