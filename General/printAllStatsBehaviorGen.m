%%

clear; close all; clc

% set folder
p1 = 'D:\AD_Project\Behavior';
cd(p1)

% define test info
nUnits = 'mice';

% define data
useFolder = 'Figures/manuscriptFigures';
subFolds = {'/all','/bySex','/bySex'};
subLabels = {'all','male','female'};
useFile = 'avgVelMovingRaw_';
nSex = length(subLabels);

% define groups
genoTypes = {'WT','PS19'};
nGeno = length(genoTypes);

% load data
dataAll = cell(nSex,nGeno);
for ss = 1:nSex
    load([useFolder subFolds{ss} '/data/' useFile subLabels{ss}],'data')
    dataAll(ss,:) = data;
end

% define analysis parameters
timeTypes = {'FE','NE'};
timeCats = {1,2:11};

% test names
testNameCorr = 'two-tailed Pearson linear correlation';
testNameTTimes = {'two-tailed unpaired Students t-test',...
    'two-tailed unpaired Students t-test, Bonferroni-Holm correction for mutiple comparsions'};
MCtimes = [0 1];

% t-test settings
pair = 0;
limitP = 1;


%% Perform statistical tests with printing

outStats = {};

for ss = 1:nSex

    for tt = 1:length(timeTypes)

        % define time label
        catLabel = [subLabels{ss} ': ' timeTypes{tt}];

        % define sub time course data
        timeData = cellfun(@(x) x(:,timeCats{tt}),dataAll(ss,:),'UniformOutput',false);

        % perform ANOVA
        if ~isscalar(timeCats{tt})
            outStats(end+1,:) = [catLabel anovaEffectSizeMinimal(timeData{1},timeData{2},nUnits)];
        end

        % perform t-test
        outStats(end+1,:) = [catLabel ttestEffectSize(timeData(1),timeData(2),testNameTTimes{tt},nUnits,pair,MCtimes(tt),limitP)];

    end
end



