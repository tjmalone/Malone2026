%% envExperience.m
% Plots the difference in experience between WT and AD mice. The average
% number of days, laps, and rewards in each envirnment type, as well as the
% age of mice at the start of imaging, are plotted. Experience data is
% loaded from experienceData.mat, which is contains manually generated
% matrices for each category.
%

clear; close all; clc


%% Set input parameters

% set correct directory
cd('D:\AD_Project\Behavior')
p = pwd;

% load experience data
load('experienceData.mat')

% load genotypes {[WT],[AD]}
load('groupIDs.mat')

% manually set genotypes
% groups = {[1 2 3 5 6 7 12],[4 8 9 10 11 13 14]};

% define environment types
envTypes = {'1m Passive','TM_6m_Env2_passive','TM_6m_Env2_active','TM_6m_Env1_active'};
nTypes = length(envTypes);

% set save folder name
svFile = [pwd '/Figures/Figures' char(datetime('today','format','yyyyMMdd')) '/experience'];
if ~isfolder(svFile)
    mkdir(svFile);
    mkdir([svFile '/data']);
end

% set save file name
svName = 'experience';


%% Process expereince data

% initialize experience by environment type matrices
dayUse = cell(1,nTypes);
lapUse = cell(1,nTypes);
rewUse = cell(1,nTypes);

% convert expereince to matrix
lapsN = cell2mat(laps);
rewN = cell2mat(rewards);

% cycle through environment types
for ee = 1:nTypes
    
    % find all environment matches
    curMatch = cellfun(@(x) contains(x,envTypes{ee}),environments);
    
    % remove extra days from novel environment
    maxDay = 11;
    if ee==4
        for ii = 1:size(curMatch,2)
            curTrue = find(curMatch(:,ii)==1);
            if length(curTrue)>=maxDay
                curMatch(curTrue(maxDay:end),ii) = 0;
            end
        end
    end

    % find environment days
    dayUse{ee} = curMatch;
    
    % find environment laps
    lapUse{ee} = nan(size(curMatch));
    lapUse{ee}(curMatch) = lapsN(curMatch==1);
    
    % find environment rewards
    rewUse{ee} = nan(size(curMatch));
    rewUse{ee}(curMatch) = rewN(curMatch==1);
    
    % apply maximum lap filter (not currently used)
%     maxLaps = 200;
%     lapUse{ee}(lapUse{ee}>maxLaps) = maxLaps;
end

% calculate total experience of each mouse by environment type
dayTots = cellfun(@(x) sum(x,1),dayUse,'UniformOutput',0);
lapTots = cellfun(@(x) nansum(x,1),lapUse,'UniformOutput',0);
rewTots = cellfun(@(x) nansum(x,1),rewUse,'UniformOutput',0);


%% Plot results

% set enviornment labels
X = {'1m','FE_{passive}','FE_{active}','NE_{active}'};

% set plot colors
colors = {[0 0 1],[1 0 0]};

% set experience type lapes
ylabs = {'Training Days (#)','Training Laps (#)','Training Rewards (#)','Imaging Age (months)'};

% set save labels
svType = {'days','laps','rewards','age'};

% set input data
useData = {dayTots,lapTots,rewTots,age};

% initialize stats parameters and array
outStats = {};
nUnits = 'mice';
testName = 'two-tailed unpaired Students t-test';
testPair = 0;
testMC = 0;
testLimitP = 0;

for f = 1:length(useData)
    figure; hold on
    
    if f~=4
        % plot experierence by environment type
        
        % split data by genotype
        data1 = cellfun(@(x) x(groups{1})',useData{f},'UniformOutput',0);
        data2 = cellfun(@(x) x(groups{2})',useData{f},'UniformOutput',0);
        data = [data1;data2]';
        
        % set significance pairs
        pShow = [(1:4)',(5:8)'];
        
        % plot data
        [h,pBar] = barGroup(X,data,'violin',colors,pShow);
        
        % set legend
        legend(h,groupIDs,'Location','NorthWest')
        legend boxoff
    else
        % plot mouse age
        
        % split data by genotype
        data = {useData{f}(groups{1}),useData{f}(groups{2})};
        
        % set significance pairs
        pShow = [1 2];
        
        % plot data
        [~,pBar] = barGroup(groupIDs,data,'violin',colors,pShow);
    end
    
    % set figure labels/settings
    xlabel('Environment')
    ylabel(ylabs{f})
    set(gca,'FontSize',12)
    
    % save figure
    savefig([svFile '/' svName '_' svType{f} '.fig'])
    close
    
    % save underlying data
    save([svFile '/data/' svName '_' svType{f} '.mat'],'data')

    % calculate full statistics
    if f~=4
        statData = data(3:4,:)';
    else
        statData = data';
    end
    testCat = svType{f};
    outStats(end+1,:) = [testCat ttestEffectSize(statData(1,:),statData(2,:),testName,nUnits,testPair,testMC,testLimitP)];

end

