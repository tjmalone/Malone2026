%%

clear; close all; clc

% set folder
p1 = 'D:\AD_Project\imagingData';
cd(p1)

% define and load data
useFolder = 'Figures/Figures_timecourse/timecourse/data/';
useFile = 'decodeIORaw_allType-cell';
load([useFolder useFile],'mapData')

% define test info
nUnits = 'FOV';

% define analysis parameters
sexTypes = {'allType','female','male'};
morphTypes = {'allType','ste','pyr'};
genoTypes = {'WT','PS19'};
timeTypes = {'All','Pre','Post','Diff'};
useSexes = 3:-1:2;
useMorphs = 1;
useTimes = [1];

% define time categories for anova and correlation
timeCats = {2:11,2,8:11,[]};
corrDays = 2:11;
tDays = 2:11;

% timeCats = {2:10,2,8:10,[]};
% corrDays = 2:10;
% tDays = 2:10;

% test names
testNameCorr = 'two-tailed Pearson linear correlation';
testNameTtime = 'two-tailed unpaired Students t-test, Bonferroni-Holm correction for mutiple comparsions';
testNameTavg = 'two-tailed unpaired Students t-test';

% t-test settings
pair = 0;
MCtime = 1;
MCavg = 0;
limitP = 0;

% define pairwise ttest settings
testNamePairwise = 'two-tailed unpaired Students t-test, Bonferroni-Holm correction for mutiple comparsions';
MCPairwsise = 1;
idxsPairwise = [1 2; 3 4; 2 3; 2 4];
labsPairwise = {'Male';'Female';'PS19 Male vs. WT Female';'PS19'};
nMulti = length(idxsPairwise);
useMulti = 1;


%% Perform statistical tests with printing

outStats = {};

for mm = 1:length(morphTypes)
    if ~ismember(mm,useMorphs); continue; end

    for tt = 1:length(timeTypes)
        if ~ismember(tt,useTimes); continue; end

        % define sub average data
        barData = squeeze(mapData.(timeTypes{tt})(useSexes,mm,:))';

        % perform t-test
        BCat = [morphTypes{mm} '-' timeTypes{tt}];
        outStats(end+1,:) = [BCat ttestEffectSize(...
            barData(1,:),barData(2,:),testNameTavg,nUnits,pair,MCavg,limitP)];

        % perfrom multi-t-test
        if useMulti==1
            multiData = reshape(barData,[],1)';
            multiCat = strcat(repmat({[BCat ' multi: ']},nMulti,1),labsPairwise);
            outStats(end+1:end+nMulti,:) = [multiCat ttestAllPairsEffectSize(...
                multiData,testNamePairwise,nUnits,pair,MCPairwsise,idxsPairwise)];
        end
    end

    for ss = 1:length(sexTypes)
        if ~ismember(ss,useSexes); continue; end

        catCur = [sexTypes{ss} '-' morphTypes{mm}];

        % define full time course
        timeAll = squeeze(mapData.timecourse(ss,mm,:));
        timeAll = cat(1,timeAll{:});

        % print correlation info
        for gg = 1:length(genoTypes)
            % get data means
            corrData = cellfun(@(x) mean(x,'omitnan'), timeAll(gg,corrDays));

            % perform test
            corrCat = [catCur '-' genoTypes{gg}];
            outStats(end+1,:) = [corrCat corrEffectSize(corrData,[],testNameCorr,'days')];
        end

        % perform t-test
        outStats(end+1,:) = [catCur ttestEffectSize(timeAll(1,tDays),timeAll(2,tDays),testNameTtime,nUnits,pair,MCtime,limitP)];

        for tt = 1:length(timeTypes)
            if ~ismember(tt,useTimes); continue; end

            if tt==4; continue; end

            % define time label
            TCat = [catCur '-' timeTypes{tt}];

            % define sub time course data
            timeData = timeAll(:,timeCats{tt});

            % perform ANOVA
            tData1 = cat(2,timeData{1,:});
            tData2 = cat(2,timeData{2,:});
            if ~isscalar(timeCats{tt})
                outStats(end+1,:) = [TCat anovaEffectSizeMinimal(tData1,tData2,nUnits)];
            else
                outStats(end+1,:) = [TCat ttestEffectSize(tData1,tData2,testNameTavg,nUnits,pair,MCavg,limitP)];
            end
        end
    end
end


