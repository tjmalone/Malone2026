%% trajectoryPlot_mouseCorr_getAll.m
% Calculate neuronal activity distribution as a function of track location
% per mouse for all days. Caluclate and analyze intermouse map correlation
%
% Distribution Types:
%   Learning types - learn
%   Run types  - all runs
%   Distribution types - to others RBR
%


%% Load data

clear; clc; close all

p1 = 'D:\AD_Project\imagingData';
cd(p1)

% load distribution data
load('data/distActivity.mat','distActivity')

% load cell selections
load('data/cellSelect.mat','cellSelect')

% set save folder name
svFile = [p1 '/Figures/Figures' char(datetime('today','format','yyyyMMdd')) '/trajectory'];
if ~isfolder(svFile)
    mkdir(svFile);
    mkdir([svFile '/data']);
end


%% Set input parameters

% define settings
learningTypes = 'learn';
distTypes = 'dfof';
sexIDsAct = 'allSex';
cellTypes = {'common','grid','nongrid'};
morphTypes = 'allMorph';
genotypes = {'WT','AD'};
nGeno = length(genotypes);

% define averaging parameters
fovPerMouse = 2;

% define day category and labels
useDays = 1:11;
dayNames = ['FE' string(1:10)];

% define colors
colors2 = {[0 0 1],[1 0 0]};
colors4 = {[102 105 255]/255,[255 102 102]/255;[0 0 175]/255,[175 0 0]/255};

% close all figures
clAll = 1;


%% Process data

for ct = 1:length(cellTypes)
    % collect use data and cells
    distBase = distActivity.(learningTypes).allRuns.(distTypes).all;
    cellsBase = cellSelect.(learningTypes).(sexIDsAct).(cellTypes{ct}).(morphTypes);
    cellsComb = cellfun(@(x,y) cat(1,x,y),cellsBase.WT,cellsBase.AD,'UniformOutput',false);

    % initialize quantification struct
    nFOV = size(cellsComb,1);
    nDays = size(cellsComb,2);
    dataInterMouse = cell(nFOV/2,nDays);
    dataInterMouseMean = zeros(nFOV/2,nDays);

    % correct empty cells
    maxIdx = cellfun(@max,cellsComb,'UniformOutput',false);
    emptyCellSel = cellfun(@isempty,maxIdx);
    emptyData = cellfun(@isempty,distBase);
    fixIdx = find(~emptyCellSel & emptyData);
    refIdx = find(~emptyData,1);
    szXCur = size(distBase{refIdx},2);
    for ii = 1:length(fixIdx)
        distBase{fixIdx(ii)} = nan(maxIdx{fixIdx(ii)},szXCur);
    end

    % define current data
    useData = cellfun(@(x,y) x(y,:),distBase,cellsComb,'UniformOutput',false);

    % fill in empty cells
    emptyAll = cellfun(@isempty,useData);
    for ff = 1:nFOV
        % get all empty and first full index
        curEmpty = find(emptyAll(ff,:));
        curFull = find(~emptyAll(ff,:),1);

        % skip full or empty rows
        if isempty(curEmpty) || isempty(curFull)
            continue
        end

        % fill empty cells
        for ii = curEmpty
            useData{ff,ii} = nan(size(useData{ff,curFull}));
        end
    end

    % fill in remaining empty cells
    curEmpty = find(cellfun(@(x) size(x,2),useData)==0);
    for ii = curEmpty'
        useData{ii} = zeros(0,szXCur);
    end

    % take mouse mean
    dataMouse = cell(nFOV/fovPerMouse,nDays);
    for ff = 1:fovPerMouse:nFOV
        curFOV = ff:ff+fovPerMouse-1;

        for dd = 1:nDays
            dataMouse{(ff+fovPerMouse-1)/fovPerMouse,dd}...
                = cat(1,useData{curFOV,dd});
        end
    end

    % take mouse mean
    dataMean = cellfun(@(x) mean(x,1,'omitnan'),...
        dataMouse,'UniformOutput',false);

    % concatenate across mice and remove nan mice
    dataCat = zeros(nFOV/fovPerMouse,nDays,120);
    for ff = 1:nFOV/fovPerMouse
        for dd = 1:nDays
            dataCat(ff,dd,:) = dataMean{ff,dd};
        end
    end

    % loop across days
    for dd = 1:nDays
        % get current data
        curDataCat = squeeze(mean(dataCat(:,dd,:),2,'omitnan'));

        % calculate correlations
        corrAll = corr(curDataCat',curDataCat');
        corrAll(eye(size(corrAll,1))==1) = NaN;

        % store data by mouse
        for mm = 1:nFOV/2
            dataInterMouse{mm,dd} = corrAll(mm,:);
            dataInterMouseMean(mm,dd) = mean(corrAll(mm,:),2,'omitnan');
        end

    end

    save(['data/interMouseData_' cellTypes{ct} '.mat'],'dataInterMouse','dataInterMouseMean')


    %% Plot inter-mouse cosistency

    % load behavior groups
    load('D:\AD_Project/Behavior/groupIDs.mat')
    nSexes = length(sexIDs);
    nGroups = length(genotypes);

    F = figure; hold on

    for ss = 2:nSexes
        for gg = 1:nGroups
            % get current group data
            curIdx = intersect(sexes{ss},groups{gg});
            plotData = dataInterMouseMean(curIdx,:);

            % get mean and SEM
            plotM = mean(plotData,1,'omitnan');
            plotS = nansem(plotData,1);

            % plot group data
            errorbar(useDays,plotM,plotS,'Color',colors4{ss-1,gg})
        end
    end

    % set labels
    title('Inter-mouse Consistency')
    xlim([0.5 11.5])
    xticks(1:length(dayNames));
    xticklabels(dayNames)
    xlabel('Learning day')
    ylabel('Inter-mouse correlation')
    set(gca,'FontSize',18)

    % save figure
    savefig([svFile '/interMouseTimecourse_' cellTypes{ct}])

    % save data
    save([svFile '/data/interMouseCorr_' cellTypes{ct} '.mat'],'dataInterMouseMean')

end