%%  behaviorAD_compareThresh
% A cut down verison of behaviorAD_all that compares percent success under
% multiple threshold conditions. Other days are calculated, but only FE day
% 0 and NE day 1-10 are plotted.
%
% Required paramters:
%   Genotype groups (currently loaded from groupIDs.mat)
%   FE and NE reward locations
%   Use laps (currently 2-20 to exclude first lap where mice are pushed)
%   Threshold parameters
%

clear; close all; clc


%% Initialize parameters

cd('D:\AD_Project\Behavior')
p = pwd;

% identify mouse folders
d = dir('TM*');
mID = {d(:).name}';

% read environment parameters
rewStart = [240 510];
rewStop = [290 560];

% use laps
startLap = 2;
stopLap = 20;

% set thresholds
speedThresh = [1 1 1 1 1];    % speed threshold
timeThresh = [60 45 30 15 6];   % stop threshold for smoothing (sec*60)
winBuff = [0 0 0 0 0];       % zone window buffer


%% Initialize behavior calculations

% behavior data struct (clear so fields don't need to be initialized)
clear B

% define day category names
typeFields = {'famLearn','novLearn','famRecall','novRecall','famPre'};


%% Calculate behavioral parameters

for m = 1:length(mID)
    %% Process file info
    
    % move to current mouse directory
    cd(mID{m})
    
    % identify all virmen files
    d = dir('*T*.txt');
    fNames = {d(:).name}';
    
    % load file classification information
    load('fileType.mat','envType','type','typeLearn','typeAll');
    
    % define which track type each file corresponds to (files with no track
    % type are extra NE days, so are set to 2)
    envTrack = envType;
    envTrack(envTrack==0) = 2;
    
    
    %% Initialize data fields
    
    % initialize struct
    curB = struct();
    
    % set struct fields
    curB.type = type;
    curB.typeLearn = typeLearn;
    curB.typeAll = typeAll;


    %% Run behavior analysis for each session
    
    % cycle through all files
    for ff = 1:length(fNames)
        %% Run analysis
        
        % load data (a .mat file will be created after the first load,
        % which will be loaded directly on future tuns)
        if ~isfile([fNames{ff}(1:end-4) '.mat'])
            logData = readLog(fNames{ff}, 9);
            params = logParams(logData);
            save([fNames{ff}(1:end-4) '.mat'],'params')
        else
            load([fNames{ff}(1:end-4) '.mat'])
        end
        
        % set environment parameters
        envParams = struct();
        envParams.rewStart = rewStart(envTrack(ff));
        envParams.rewStop = rewStop(envTrack(ff));
        envParams.speedThresh = speedThresh;
        envParams.timeThresh = timeThresh;
        envParams.winBuff = winBuff;

        % set use laps
        if envTrack(ff)==2
            % set use laps for novel environment days
            useLaps = startLap:stopLap;
        else
            % set use laps for familiar environment days
            useLaps = startLap:1000;
        end

        % calculate behavior
        behavior = extractPerCorrect(params,envParams,useLaps);
        curB.perSuccess(ff,:) = behavior;

        % extract mouse and day
        curB.day{ff} = fNames{ff}(3:8);
        curB.mouse = split(mID{m},'\');

        % print file completion
        fprintf([num2str(m) '-' num2str(ff) '\n'])

        
    end
    
    %% Save results
    
    % save current mouse data
    B(m) = curB;
    
    % change to parent directory
    cd(p)
    
end

% save combined data (match to load in next section)
save('data\allMiceThresh.mat','B')


%% Set plotting parameters
% To skip recalculation, start running from here. All previously calculated
% data is loaded in allMiceThresh.mat. This file name must match the current
% save output above.
%

% reclear all variable.
clear; close all; clc

% reset directory
cd('D:\AD_Project\Behavior')
p = pwd;

% load processed behavioral data
load('data\allMiceThresh.mat')

% load genotypes {[WT],[AD]}
load('groupIDs.mat')

% manually set genotypes
% groups = {[1 2 3 5 6 7 12],[4 8 9 10 11 13 14]};

% define legend
legs = groupIDs;

%define plot colors
colors = {[0 0 1],[1 0 0]};

% set day categories (FE,NE,FER,NER)
dayCats = {1:4,5:14,15,16};
nCats = length(dayCats);

% set save folder name
svFile = [p '\Figures\Figures' char(datetime('today','format','yyyyMMdd'))];
mkdir(svFile);
mkdir([svFile '\data']);


%% Plot line graphs by session

% define fields, titles, and labels
nameThresh = {'1','0-75','0-5','0-25','0-1'};

% calculate maximum number of plotting sessions
maxSessions = max(cellfun(@length,{B(:).typeAll}));

% close all figures
clAll = 1;

% define x-axis to plot
xPlot = 4:14;
nPlot = length(xPlot);
X = xPlot-4;
Xname = ['FE' string(1:nPlot-1)];

% cycle through all fields
for f = 1:length(nameThresh)
    %% Plot current field

    % initialize data matrices
    curData = zeros(length(B),maxSessions);
    data = cell(1,2);
    
    % extract data for current field (pad data with nan)
    for ii = 1:length(B)
        curD = B(ii).perSuccess(B(ii).typeAll,f);
        curD(end+1:maxSessions) = nan;
        curData(ii,:) = curD;
    end
    
    % sort data by group
    for g = 1:2
        data{g} = curData(groups{g},xPlot);
    end
    
    % calculate p values
    pMC = nan(nPlot,1);
    pAnova = nan(nPlot,2);
    
    for ii = 1:nCats
        % find intersect of plotting x values with day category 
        curCat = find(ismember(xPlot,dayCats{ii}));
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
    figure; hold on
    
    % plot line graph
    plotErrorSig(X,data{1},data{2},legs,pAnova,pMC,colors)
    
    % set figure labels
    legend('Location','eastoutside')
    xlabel('Session')
    ylabel('Percent')
    title(['Percent successful run (' replace(nameThresh{f},'-','.') 's window)'])
    set(gca,'FontSize',12)
    
    % set x axis and labels
    set(gca,'XTick',X)
    set(gca,'XTickLabels',Xname)
    xlim([X(1)-0.5 X(end)+0.5])
    ylim([0 120])

    % save figure
    savefig([svFile '\perSuccess_' nameThresh{f} '.fig'])
    
    % save data
    save([svFile '\data\perSuccess_' nameThresh{f} '.mat'],'xPlot','data')
    
    % close figure
    if clAll
        close
    end
end



