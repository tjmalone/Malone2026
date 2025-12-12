%% plotCellNumber_forFigs
% Plot active cell number and proportion of verous functional cell types
% for use in figures
%


%% Normalize cell numbers

clear; close all; clc

p1 = 'D:\AD_Project\imagingData';
cd(p1)

% load alignment data
load('foldersLearning.mat')
load('data/globalCellTypes_PreComb.mat', 'globalCellTypeIdx','commonTypes')
nFOV = length(foldersLearning);
useDay = 1;

% set save folder name
svFile = [p1 '/Figures/Figures' char(datetime('today','format','yyyyMMdd')) '/cellNumbers'];
mkdir(svFile);
mkdir([svFile '/data']);


%% Calculate cell type information

% % define final groups (full version)
% groupFields = {4,1,2,3,1:3,7,5,6,5:6,8,9,10};
% groupNames = {'grid','cueL','cueR','cueAll','cueCombInd','cueCombComb',...
%     'speedPos','speedNeg','speedCombInd','speedCombComb','othersInd','othersComb'};

% define final groups
groupFields = {4,7,8,10};
groupNames = {'grid','cue','speed','uncategorized'};
nGroups = length(groupFields);

perGlobal = zeros(nFOV,nGroups);

for ff = 1:nFOV
    % get common cell number
    nCommon = size(alignsLearning{ff}(:,useDay),1);

    for gg = 1:nGroups
        % get current sub-fields
        curFields = groupFields{gg};
        nFields = length(curFields);

        % get current group indices
        idxCur = [];
        for sf = 1:nFields
            idxSub = globalCellTypeIdx.(commonTypes{curFields(sf)}){ff}(:,useDay);
            idxCur = cat(1,idxCur,idxSub);
        end

        % store cell type percentage
        perGlobal(ff,gg) = length(idxCur)/nCommon*100;
    end
end


%% Plot cell numbers

% load genotypes {[WT],[AD]}
load('groupIDs.mat')
nGeno = length(groups);
useSexes = 1:3;

% define plot colors
colors = {[0 0 1],[1 0 0]};
pShow = generatePShow(nGroups,[1 2]);

% define loop inputs
ylabs = 'Percent of all cells';
svName = 'typeTrue';

% close all figures
clAll = 0;

% initialize stats parameters and array
% initialize stats parameters and array
outStats = {};
testUnits = {'FOV','mice'};
testName = 'two-tailed unpaired Students t-test';
testPair = 0;
testMC = 0;
testLimitP = 0;
nUnits = length(testUnits);

% cycle through all inputs
for ss = useSexes
    %% Plot current field

    % sort data by group
    plotData = cell(nGroups,nGeno);
    for gn = 1:nGeno
        for gg = 1:nGroups
            plotData{gg,gn} = perGlobal(intersect(groups{gn},sexes{ss}),gg);
        end
    end

    for uu = 1:nUnits
        if strcmp(testUnits{uu},'mice')
            plotDataUse = cellfun(@(x) mean([x(1:2:end),x(2:2:end)],2,'omitnan'),...
                plotData,'UniformOutput',false);
        elseif strcmp(testUnits{uu},'FOV')
            plotDataUse = plotData;
        else
            error('Invalid unit type')
        end

        % plot violin graph
        figure; hold on
        h = barGroup(groupNames,plotDataUse,'violin',colors,pShow);

        % set figure labels
        legend(h,groupIDs,'Location','best')

        % legend boxoff
        ylabel(ylabs)
        curTtl = [sexIDs{ss} ' (for ' testUnits{uu} ')'];
        title(curTtl)
        set(gca,'FontSize',12)
        ylim([0 100])

        % calculate full statistics
        outStats(end+1,:) = [curTtl ttestEffectSize(...
            plotDataUse(:,1),plotDataUse(:,2),testName,testUnits{uu},testPair,testMC,testLimitP)];

        % save figure
        savefig([svFile '/' svName '_' sexIDs{ss} '_' testUnits{uu} '.fig'])

        % save data
        save([svFile '/data/' svName '_' sexIDs{ss} '_' testUnits{uu} '.mat'],'groupNames','plotDataUse')

        % close figure
        if clAll
            close
        end
    end
end