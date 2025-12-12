%%

clear; close all; clc

% set folder
p1 = 'D:\AD_Project\Behavior';
cd(p1)

% define test info
nUnits = 'mice';

% define data
useFolder = 'Figures/manuscriptFigures/speedMatch';
useFile = 'velocityMatchbySession.mat';

% define groups
genoTypes = {'WT','PS19'};
nGeno = length(genoTypes);

% load data
load([useFolder '/data/' useFile],'storeData')

% define sexes
sexIDs = fieldnames(storeData);
nSex = length(sexIDs);

% define analysis parameters
timeTypes = {'FE','NE'};
timeCats = {1,2:11};

% test names
testNameA = 'ordinary two-way ANOVA';
testNameTTimes = {'two-tailed unpaired Students t-test',...
    'two-tailed unpaired Students t-test, Bonferroni-Holm correction for mutiple comparsions'};
MCtimes = [0 1];

% t-test settings
pair = 0;
limitP = 1;


%% Perform statistical tests with printing

outStats = {};

for ss = 1:nSex
    % define current sex data
    sexData = storeData.(sexIDs{ss});
    fields = fieldnames(sexData);

    for ff = 1:length(fields)

        for tt = 1:length(timeTypes)

            % define time label
            catLabel = [fields{ff} '-' sexIDs{ss} ': ' timeTypes{tt}];

            % define sub time course data

            % get current field data
            curData = sexData.(fields{ff})(timeCats{tt},:);

            % perform ANOVA
            if ~isscalar(timeCats{tt})
                outStats(end+1,:) = [catLabel anovaEffectSize(...
                    curData(:,1)',curData(:,2)',testNameA,nUnits,[])];
            end

            % perform t-test
            outStats(end+1,:) = [[catLabel ': mult'] ttestEffectSize(...
                curData(:,1),curData(:,2),testNameTTimes{tt},nUnits,pair,MCtimes(tt),limitP)];

        end
    end

end

