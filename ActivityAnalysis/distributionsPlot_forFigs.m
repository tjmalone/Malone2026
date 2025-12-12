%% distributionsPlot_forFigs.m
% Plot neuronal activity distributions as a function of track location. A
% simplified version of distributionsPlot_master. No longer plots
% success/fail runs or distribution quantification
%
% Distribution Types:
%   Learning types - learn, recall
%   Run types  - all runs, true succes/fail, conditional success/fail
%   Distribution types - dfof, field, to mean RBR, to next RBR, to others
%       RBR
%
% Data Separations:
%   Sex - allSex, female, male
%   Cell type - allType, common, grid, cue
%   Morphology type - allMorph, ste, pyr
%   Genotype - WT, AD
%


%% Load data

clear; clc; close all

p1 = 'D:\AD_Project\imagingData';
cd(p1)

% load distribution data
load('data\distActivity.mat','distActivity')

% load cell selections
load('data\cellSelect.mat','cellSelect')

% set save folder name
svFile = [p1 '\Figures\Figures' char(datetime('today','format','yyyyMMdd')) '\distributionsDfof'];
if ~isfolder(svFile)
    mkdir(svFile);
    mkdir([svFile '\data']);
end


%% Set input parameters

% distribution key: learning type (lt), run type (rp), distribution type,
% (dt), success/fail (ds)
learningTypes = {'learn'};
runType = 'allRuns';
runSubType = 'all';
distTypes = {'dfof'};
% distTypes = {'dfof','field'};
distTypes = {'fieldAmp'};

% cell selection key: sex (ss), cell type (ct), morphology type (mt),
% genotype (gg)
sexIDs = {'allSex','female','male'};
cellType = 'common';
morphTypes = {'allMorph','ste','pyr'};
morphTypes = {'allMorph'};
genotypes = {'WT','AD'};
nSexes = length(sexIDs);
nMorphs = length(morphTypes);
nGeno = length(genotypes);

% define averaging methods (by cell/by FOV)
avgMethod = 'cell';

% set base plot labels
ylabBase = {'\DeltaF/F','Percent of cells with field','RBR correlation',...
    'RBR correlation','RBR correlation','EMD'};
ttlBase = {'\DeltaF/F Distribution','Field Distribution',...
    'mRBR Consistiency Distribution','nRBR Consistiency Distribution',...
    'oRBR Consistiency Distribution','Earth Movers Distance Distribution'};
xBase = {linspace(2.5,597.5,120),linspace(2.5,597.5,120),...
    linspace(12.5,587.5,116),linspace(12.5,587.5,116),...
    linspace(12.5,587.5,116),linspace(12.5,587.5,116)};

% initialize day category matrix
daysDef = struct();

% define learning day categories and labels
daysDef.learn.dayCats = {(8:11)};
daysDef.learn.ttlsCats = {'All NE'};
daysDef.learn.svLabelCats = {'NE-All'};
daysDef.learn.dayNames = ['FE' string(1:10)];

% define recall day categories and labels
daysDef.recall.dayCats = {1,2};
daysDef.recall.ttlsCats = {'FE-Ref','FE-Recall'};
daysDef.recall.svLabelCats = {'FE-Ref','FE-Recall'};
daysDef.recall.dayNames = ["FE-ref", "FE-recall"];

% define colors
colors6 = {[0 0 1],[1 0 0];[102 105 255]/255,[255 102 102]/255;[0 0 175]/255,[175 0 0]/255};
cueTypes = {'Left','Right'};

% whether to replot figures\
replotOn = 1;

% close all figures
clAll = 0;


%% Define cue templates

% define cue/reward info
colorsCue = [0 0 0;0.5 0.5 0.5];
colorsRew = [0.75 1 1];
cueX = 2.5:5:597.5;

% define cue template
cueData = struct();
for ii = 1:2
    if ii==1
        cueFolder = 'D:\AnalysisCode\PostAnalysis\Cues\6mEnv1\';
        cueData(ii).rewLoc = [240 290];
        cueData(ii).subIdxs = 1;
    else
        cueFolder = 'D:\AnalysisCode\PostAnalysis\Cues\6mEnv2\';
        cueData(ii).rewLoc = [510 560];
        cueData(ii).subIdxs = 2:11;
    end

    % load cue templates
    load([cueFolder 'tempL.mat'])
    load([cueFolder 'tempR.mat'])
    cueData(ii).cueTemp = [tempL,tempR];
end


%% Plot distributions

% loop through learning types
for lt = learningTypes
    % set current day categories and labels
    curDayCats = daysDef.(lt{:}).dayCats;
    curTtlsCats = daysDef.(lt{:}).ttlsCats;
    curSvLabelCats = daysDef.(lt{:}).svLabelCats;
    curDayNames = daysDef.(lt{:}).dayNames;
    nDayCats = length(curDayCats);

    % find day category titles
    dayLabels = cell(1,nDayCats);
    for dc = 1:nDayCats
        curCat = curDayCats{dc};
        dayLabels{dc} = strcat({'Days '},curDayNames(curCat(1)),'-',curDayNames(curCat(end)));
    end

    % loop through distribution types
    for dt = 1:length(distTypes)
        % get current data
        curDist = distActivity.(lt{:}).(runType).(distTypes{dt}).(runSubType);

        % loop through sexes
        for ss = 1:nSexes

            % loop through morphology types
            for mt = 1:nMorphs

                % get morphology types
                curTypeBase = cellSelect.(lt{:}).(sexIDs{ss}).(cellType).(morphTypes{mt});

                %% Plot indivivudal figure

                % check for completion
                ttlCur = [lt{:} '-' runType '-' distTypes{dt}...
                    '-' sexIDs{ss} '-' cellType '-' morphTypes{mt} '-' avgMethod];
                if replotOn==0 && isfile([svFile '\' ttlCur '.fig'])
                    continue
                end

                % initialize tiled figure
                figure
                t = tiledlayout(1, nDayCats);
                sgtitle(ttlCur)

                % initialize figure axes
                for ii = 1:nDayCats
                    nexttile();
                    title([curTtlsCats{ii} ': ' char(dayLabels{ii})])
                end

                % loop through genotypes
                for gg = 1:nGeno
                    %% Process current data

                    % get current cell selection
                    curCellSel = curTypeBase.(genotypes{gg});

                    % correct empty cells
                    maxIdx = cellfun(@max,curCellSel,'UniformOutput',false);
                    emptyCellSel = cellfun(@isempty,maxIdx);
                    emptyData = cellfun(@isempty,curDist);
                    fixIdx = find(~emptyCellSel & emptyData);
                    refIdx = find(~emptyData,1);
                    szXCur = size(curDist{refIdx},2);
                    for ii = 1:length(fixIdx)
                        curDist{fixIdx(ii)} = nan(maxIdx{fixIdx(ii)},szXCur);
                    end

                    % define current data
                    curData = cellfun(@(x,y) x(y,:),curDist,...
                        curCellSel,'UniformOutput',false);
                    nFOV = size(curData,1);
                    nDays = size(curData,2);

                    % take FOV mean
                    if strcmp(avgMethod,'fov')
                        curDataMean = cellfun(@(x) mean(x,1,'omitnan'),...
                            curData,'UniformOutput',false);
                    elseif strcmp(avgMethod,'cell')
                        curDataMean = curData;
                    else
                        error('Invalid averaging method')
                    end

                    % fill in empty cells
                    emptyAll = cellfun(@isempty,curDataMean);
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
                            curDataMean{ff,ii} = nan(size(curDataMean{ff,curFull}));
                        end
                    end

                    % fill in remaining empty cells
                    curEmpty = find(cellfun(@(x) size(x,2),curDataMean)==0);
                    for ii = curEmpty'
                        curDataMean{ii} = zeros(0,szXCur);
                    end

                    % concatenate data
                    curDataCat = cell(1,nDays);
                    for dd = 1:nDays
                        curDataCat{dd} = cat(1,curDataMean{:,dd});
                    end


                    %% Plot distributions

                    % initialize distribution data
                    distData = cell(1,nDayCats);

                    % loop through day categories
                    for dc = 1:nDayCats
                        % set current axis
                        nexttile(dc); hold on

                        % combine days of the same day category
                        distData{dc} = mean(cat(3,curDataCat{curDayCats{dc}}),3,'omitnan');

                        % normalize data
                        plotMean = mean(distData{dc},1,'omitnan');
                        plotNorm = (distData{dc}-min(plotMean))/range(plotMean);

                        % plot distribution
                        semshade(plotNorm,0.3,colors6{ss,gg},xBase{dt});
                    end
                end


                %% Complete distribution plots

                % ylimits
                ylimsCur = zeros(nDayCats,2);

                for dc = 1:nDayCats
                    % store y limits
                    nexttile(dc); hold on
                    ylimsCur(dc,:) = ylim;
                end

                % find y limit extremes
                ylimsMin = min(ylimsCur(:,1),[],1);
                ylimsMax = max(ylimsCur(:,2),[],1);

                % loop through day categories
                for dc = 1:nDayCats
                    % set current axis
                    nexttile(dc); hold on

                    % set y limits
                    ylim([ylimsMin ylimsMax])

                    if strcmp(lt,'learn') && curDayCats{dc}(1)>1
                        cueTemp = cueData(2).cueTemp;
                        rewLoc = cueData(2).rewLoc;
                    else
                        cueTemp = cueData(1).cueTemp;
                        rewLoc = cueData(1).rewLoc;
                    end

                    % plot cues
                    for ii = 1:size(cueTemp,2)
                        plotCues(cueX,cueTemp(:,ii),max(ylim),colorsCue(ii,:),min(ylim));
                    end

                    % plot rewards
                    pos = [rewLoc(1) min(ylim) diff(rewLoc) diff(ylim)];
                    rectangle('Position',pos,'FaceColor',colorsRew,...
                        'EdgeColor','none')

                    % flip children order to plot reward behind distributions
                    axCur = gca;
                    axCur.Children = flipud(axCur.Children);
                end

                %% Save current figure

                % save figure
                savefig([svFile '\' ttlCur '.fig'])

                % % save data
                save([svFile '\data\' ttlCur '.mat'],'curDataCat')

                % close figure
                if clAll
                    close
                end

            end
        end
    end
end

