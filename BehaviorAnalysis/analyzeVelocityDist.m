%%  analyzeVelocityDist
% Analyzes the velocity distributions calculated in behaviorAD_all.
% Requires the saved output array. Currently runs a separate code to add
% flexibity.
%

clear; close all; clc


%% Initialize parameters

cd('D:\AD_Project\Behavior')
p = pwd;

% read environment parameters
rewStart = [240 510];
rewStop = [290 560];
trackLength = 600;

% load processed behavioral data
load('data\allMiceDataB4.mat')

% load genotypes {[WT],[AD]}
load('groupIDs.mat')

% manually set genotypes
% groups = {[1 2 3 5 6 7 12],[4 8 9 10 11 13 14]};

% define legend
legs = genotypes;

%define plot colors
colors = {[0 0 1],[1 0 0]};

% set save folder name
svFile = [p '\Figures\Figures' char(datetime('today','format','yyyyMMdd'))];
mkdir(svFile);
mkdir([svFile '\data']);

cueFolder = 'D:\AnalysisCode\PostAnalysis\Cues\6mEnv2\';
load([cueFolder 'tempL.mat'])
load([cueFolder 'tempR.mat'])
cueLR = max(tempR,tempL);
cueX = (0:5:595)';

corrX = (7:91)';


%% Calculate behavioral parameters

% define fields, titles, and labels
fields = {'velDistMean','dwellDistMean','stopDist'};
ttls = {'Velocity distribution','Dwell time distribution','Stop distribution'};
xlabs = 'Track location';
ylabs = {'Velocity','Percent time spent','Percent stops'};
X = {0:5:600,0:5:600,0:5:600};

dayNames = [string(1:10)];
daySet = {8:10};

% close all figures
clAll = 0;

for f = 1:length(fields)
    %% Calculate distribution

    meanDistAll = zeros(length(B),length(X{f}));

    for m = 1:length(B)

        dataUse = B(m).(fields{f});
        learnDays = B(m).typeLearn;

        % calculate average distribution or distribution difference
        for ii = 1:length(daySet)
            curDist = mean(cat(1,dataUse{learnDays(daySet{ii})}),1,'omitnan');

            if ii==1
                meanDist = curDist;
            else
                % only makes sense for two daySets
                meanDist = meanDist-curDist;
            end
        end

        % store distribution
        meanDistAll(m,:) = meanDist;
    end


    %% Group data by category

    % initialize data matrix
    data = cell(1,2);

    % sort data by group
    for g = 1:2
        data{g} = meanDistAll(groups{g},:);
        % data{g} = data{g}(:,1:length(X{f}));
    end


    %% Calculate statistics

    % calculate p values
    [pACur,pMC] = anovaRM2W_full_BH(data{1},data{2},0);
    pAnova = pACur([1 3])'.*ones(length(X{f}),2);
   
    % plot correlation with cue template
    corrData = cell(1,2);
    for ii = 1:2
         corrData{ii} = corr(cueLR(corrX),data{ii}(:,corrX)');
    end
    figure; hold on
    [h,pBar] = barGroup([],corrData,colors,[1 2]);


    %% Plot distributions

    % initialize figure
    figure; hold on

    % plot line graph
    plotErrorSig(X{f},data{1},data{2},legs,pAnova,pMC,colors)

    % plot reward zone
    ymin = min(ylim);
    ymax =  max(ylim);
    fill([510 510 560 560],[ymin ymax ymax ymin],[1 1 0.7],'LineStyle','None',...
         'FaceAlpha',0.5,'DisplayName','Reward Zone')
    legend('Location','eastoutside')

    plot(cueX,tempL*(ymax-ymin)-ymin,'Color',[0.5 0.5 0.5])
    plot(cueX,tempR*(ymax-ymin)-ymin,'Color',[0 0 0])

    % set figure labels
    set(gca,'FontSize',12)
    xlabel(xlabs)
    ylabel(ylabs{f})

    % define title
    if length(daySet)>1
        ttlType = ' Difference';
    else
        ttlType = '';
    end

    dayLabels = [];
    for ii = 1:length(daySet)
        if ii>1
            dayLabels = [dayLabels ' vs. '];
        end
        dayLabels = [dayLabels strcat({'Days '},dayNames(daySet{ii}(1)),'-',dayNames(daySet{ii}(end)))];
    end
    dayLabels = cat(2,dayLabels{:});

    % set title
    ttlUse = [ttls{f} ttlType ': ' dayLabels];
    title(ttlUse)
    

    % set x axis
    xlim([X{f}(1) X{f}(end)])

    % save figure
    savefig([svFile '\' fields{f} '.fig'])

    % save data
    curX = X{f};
    save([svFile '\data\' fields{f} '.mat'],'curX','data')

    % close figure
    if clAll
        close
    end

end

