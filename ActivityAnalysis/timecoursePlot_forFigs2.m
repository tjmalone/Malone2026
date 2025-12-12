%% timecoursePlot_master.m
% Plot neuronal activity variables as a function of day and generate heat
% maps of the relative differences by subgroup
%
% Activity Types:
%   Learning types - learn, recall (not implemented)
%   Activity variables (implemented) - intra-day consistency, inter-day
%       consistency
%   Activity variables (not implemented) - activity correlation to cue
%       template, spatial decoding, cue decoding, activity variation
%   Run types  - all runs, true succes/fail, conditional success/fail
%   Distribution types - dfof, field, to mean RBR, to next RBR, to others
%       RBR
%
% Data Separations:
%   Sex - allSex, female, male
%   Cell type - allType (not implemented), common, grid (not implemented),
%       cue (not implemented)
%   Morphology type - allMorph, ste, pyr
%   Genotype - WT, AD
%
% Heatmap Types: pre-learning, post-learning, learning difference, all days
%


%% Load data

clear; clc; close all

p1 = 'D:\AD_Project\imagingData';
cd(p1)

% initialize data struct
globalData = struct();
globalFields = {'RBR','inter','dfofRdgBkg','FWidthSum','spatialSelectivity',...
    'decodeIORaw','decodeI','decodeO','decodeAll','dfofSig','speedScoreAbs',...
    'dfofInField','dfofNonField','dfofI','dfofO','dfofIO','speedScoreRaw',...
    'speedScorePos','speedScoreNeg'};
globalTitles = {'Intra-day consistency (to others)','Inter-day consistency',...
    'Ridge/Background Ratio of \DeltaF/F','Total Field Coverage','Spatial Selectivity',...
    'Decoding IO ratio (raw)','In-Cue Decoding','Out-Cue Decoding','Overall Decoding',...
    'Mean \DeltaF/F_{sig} (while moving)','Speed Score','In-Field \DeltaF/F',...
    'Non-Field \DeltaF/F','\DeltaF/F In-Cue','\DeltaF/F Out-Cue',...
    'In/Out Ratio of \DeltaF/F','Raw Speed Score','Positive Speed Score',...
    'Negative Speed Score'};
globalYLabs = {'RBR correlation','inter-day correlation','rdg/bkg ratio',...
    'Field coverage (cm)','Spatial selectivity','IO ratio','% correct bins',...
    '% correct bins', '% correct bins','Mean \DeltaF/F','speed score',...
    '\DeltaF/F (%)','\DeltaF/F (%)','\DeltaF/F (%)','\DeltaF/F (%)',...
    '\DeltaF/F (%)','speed score','speed score','speed score'};
nGlobalFields = length(globalFields);

% define which heatmaps to generate for each global field
globalCats = {[1 3 4], [1 3 4], 1, 1, 1, [1 4], 1, 1, [3], [1 2 3], 1, 1, 1, 1, 1, [1 4], 1, 1, 1};

% define which global fields use common cell referencing or are
% pre-calculated by FOV. 0 = original data by all cells. 1 = original data
% by common cel. -1 = original data by FOV.
globalCommonRef = [0, 1, -1, 0, 0, -1, -1, -1, -1, 0, 0, 0, 0, -1, -1, -1, 0, 0, 0];

% define whether a high value is expected to align with good or poor
% behavior. 1 = low AD implies deficit, -1 = high AD implies decifit
globalHeatSign = [1, 1, -1, -1, 1, -1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, 1, 1, 1];

% define whether to calculate correlations for time course
corrCats = {2:11,2:10};
globalPlotCorr = [1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
corrDays = (2:11)';

% define anova classes
anovaCats = {{8:11},{2:10},{2:11},{2:5 8:11},{2 8:11}};
globalAnovaCat = [5, 2, 3, 3, 3, 3, 4, 4, 3, 5, 3, 3, 3, 3, 3, 3, 3, 3, 3];

% load intra-day consistency data
load('data/intraDayData.mat','dataActIntra')
globalData.RBR = dataActIntra.full.RBR_sigOthers;
globalData.RBR = cellfun(@(x) x',globalData.RBR,'UniformOutput',false);

% load inter-day consistency data
load('data/interDaySubset.mat','dataOut')
globalData.inter = dataOut;

% % load correlation to cue template data
load('data/corrInfoGlobal.mat','dataActivity')
globalData.dfofRdgBkg = dataActivity.dfofRdgBkg;
globalData.dfofI = dataActivity.dfofI;
globalData.dfofO = dataActivity.dfofO;
globalData.dfofIO = dataActivity.dfofIO;

% load by cell properties
load('data/dataByCell.mat','dataAll')
globalData.speedScoreRaw = dataAll.speedScore;
globalData.speedScoreAbs = cellfun(@(x) abs(x),dataAll.speedScore,'UniformOutput',false);
globalData.speedScorePos = dataAll.speedScore;
globalData.speedScoreNeg = dataAll.speedScore;
for ii = 1:numel(dataAll.speedScore)
    globalData.speedScorePos{ii}(dataAll.speedScore{ii}<=0) = NaN;
    globalData.speedScoreNeg{ii}(dataAll.speedScore{ii}>=0) = NaN;
end

globalData.spatialSelectivity = dataAll.spatialSelectivity;
globalData.FWidthSum = dataAll.FWidthSum;
globalData.dfofInField = dataAll.dfofInField;
globalData.dfofNonField = dataAll.dfofNonField;

% load decoding properties
load('data/decodeInfoGlobal_newDecodeIO.mat','decodingActivity')
globalData.decodeIORaw = decodingActivity.decodeIORaw;
globalData.decodeI = decodingActivity.decodeI;
globalData.decodeO = decodingActivity.decodeO;
globalData.decodeAll = decodingActivity.decodeAll;

% load dfof
load('data/dfofData.mat','dataDfof')
globalData.dfofSig = dataDfof;
globalData.dfofSig = cellfun(@(x) x',globalData.dfofSig,'UniformOutput',false);

% load cell selections
load('data/cellSelect.mat','cellSelect')
cellSelect = cellSelect.learn;

% load alignments
load('foldersLearning.mat','alignsLearning','trueDays')

% set save folder name
svFile = [pwd '/Figures/Figures' char(datetime('today','format','yyyyMMdd')) '/timecourse'];
if ~isfolder(svFile)
    mkdir(svFile);
    mkdir([svFile '/data']);
end


%% Set input parameters

% cell selection key: sex (ss), cell type (ct), morphology type (mt),
% genotype (gg)
sexIDs = {'allSex','female','male'};
nSex = length(sexIDs);
cellTypes = {'common','grid','nongrid'};
% cellTypes = {'common'};

morphTypes = {'allMorph','ste','pyr'};
nMorph = length(morphTypes);
genotypes = {'WT','AD'};
nGeno = length(genotypes);

% define averaging method (by cell/by FOV)
avgMethods = 'cell';

% define learning day categories and labels
catDays = {2:11,2,8:11};
catTtls = {'All NE','NE Pre-Learning','NE Post-Learning','Learning Difference'};
catSvLabels = {'NE','NE-Pre','NE-Post','NE-Diff'};
catFieldLabels = {'All','Pre','Post','Diff'};
dayNames = ['FE' string(1:10)];
nDayCats = length(catDays);

% define colors and patterns
colors2 = {[0 0 1],[1 0 0]};
colors4 = {[102 105 255]/255,[255 102 102]/255;[0 0 175]/255,[175 0 0]/255};
patterns = {'-','--'};

% whether to replot figures
replotOn = 1;

% close all figures
clAll = 1;


%% Plot activity variables

% define indices with forced all cells
allOvrIdx = 6:9;

% loop through data types
useFields = 1:nGlobalFields;
useFields = nGlobalFields-1:nGlobalFields;
useFields = [10];
% cellTypes = {'common'};

for gf = useFields
    % get current data
    fieldData = globalData.(globalFields{gf});

    % loop through cell types
    for ct = cellTypes
        %% Plot individual figure

        % overwrite decoding cell type
        if ismember(gf,allOvrIdx) && strcmp(ct,'common')
            ct = {'allType'};
        end

        % check for completion
        ttlCur = [globalTitles{gf} ': ' ct{:} ' cells (by ' avgMethods ')'];
        svLabelCur = ['/' globalFields{gf} '_' ct{:} '-' avgMethods];
        if replotOn==0 && isfile([svFile '/' svLabelCur '_timecourse.fig'])
            continue
        end

        % initialize tiled figure
        figure
        tiledlayout(3, 3);
        sgtitle(ttlCur)

        % initialize figure axes
        ax = cell(nSex,nMorph);
        mapData = struct();
        mapData.timecourse = cell(nSex,nMorph,nGeno);
        for flc = 1:length(catFieldLabels)
            mapData.(catFieldLabels{flc}) = cell(3,3,2);
        end

        % initialize correlation info
        mapData.corrR = zeros(nSex,nMorph,nGeno);
        mapData.corrP = zeros(nSex,nMorph,nGeno);

        % loop through sexes
        for ss = 1:nSex

            % loop through morphology types
            for mt = 1:nMorph

                % get current cell selections (skip invalid combinations)
                try
                    curCellSub = cellSelect.(sexIDs{ss}).(ct{:}).(morphTypes{mt});
                    if ismember(gf,allOvrIdx) && mt~=1
                        error('Skip')
                    end
                catch
                    disp('Skipped comparison)')
                    continue
                end

                % set current axis
                if ss==1
                    ax{ss,mt} = nexttile([1 1]); hold on
                    title([sexIDs{ss} '-' morphTypes{mt}])
                elseif ss==2
                    ax{2,mt} = nexttile([2 1]); hold on
                    title([morphTypes{mt} ' by sex'])
                else
                    axes(ax{2,mt}); hold on
                end

                % loop through genotypes
                for gg = 1:nGeno
                    %% Process current data

                    % get current cell selection and make data copy
                    curCellSel = curCellSub.(genotypes{gg});
                    fieldDataSub = fieldData;

                    if globalCommonRef(gf)==-1
                        curDataMean = fieldDataSub.(ct{:}).(sexIDs{ss}).(morphTypes{mt}).(genotypes{gg});
                    else
                        % convert cell selections to logical indices of common cells if necessary
                        if globalCommonRef(gf)==1
                            curCellSel = idx2common(curCellSel,alignsLearning);
                        end

                        % trim cell selections if necessary
                        if size(fieldDataSub,2)<size(curCellSel,2)
                            curCellSel = curCellSel(:,1:size(fieldDataSub,2));
                        end

                        % correct empty cells
                        maxIdx = cellfun(@max,curCellSel,'UniformOutput',false);
                        maxIdx2 = cellfun(@length,curCellSel,'UniformOutput',false);

                        emptyCellSel = cellfun(@isempty,maxIdx);
                        emptyData = cellfun(@isempty,fieldDataSub);
                        fixIdx = find(~emptyCellSel & emptyData);
                        refIdx = find(~emptyData,1);
                        for ii = 1:length(fixIdx)
                            fieldDataSub{fixIdx(ii)} = nan(maxIdx2{fixIdx(ii)},1);
                        end

                        % define current data
                        curData = cellfun(@(x,y) x(y,:),fieldDataSub,curCellSel,'UniformOutput',false);

                        % take FOV mean
                        if strcmp(avgMethods,'fov')
                            curDataMean = cellfun(@(x) mean(x,1,'omitnan'),curData,'UniformOutput',false);
                        elseif strcmp(avgMethods,'cell')
                            curDataMean = curData;
                        else
                            error('Invalid averaging method')
                        end
                    end

                    % get FOV size
                    nFOV = size(curDataMean,1);
                    nDays = size(curDataMean,2);

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
                        curDataMean{ii} = zeros(0,1);
                    end

                    % concatenate data
                    curDataCat = cell(1,nDays);
                    for dd = 1:nDays
                        curDataCat{dd} = cat(1,curDataMean{:,dd});
                    end


                    %% Plot time courses

                    % define color
                    if ss==1
                        useColor = colors2{gg};
                        usePattern = '-';
                    else
                        useColor = colors4{ss-1,gg};
                        usePattern = patterns(ss-1);
                    end

                    % calculate mean and standard error
                    catMean = cellfun(@(x) mean(x,'omitnan'),curDataCat);
                    catSEM = cellfun(@(x) nansem(x,1),curDataCat);

                    % plot line graph
                    errorbar(catMean,catSEM,'Color',useColor)


                    %% Store data

                    for flc = 1:length(catFieldLabels)
                        % take average for current day groups
                        if flc<4
                            % crop category indices
                            idxUse = catDays{flc};
                            idxUse(idxUse>size(curDataCat,2)) = [];

                            % take category mean
                            storeData = mean(cat(2,curDataCat{idxUse}),2,'omitnan');
                        else
                            % crop category indices
                            idxUse = cell(1,2);
                            for ii = 1:2
                                idxUse{ii} = catDays{ii+1};
                                idxUse{ii}(idxUse{ii}>size(curDataCat,2)) = [];
                            end

                            % take difference of category means
                            storeData = mean(cat(2,curDataCat{idxUse{2}}),2,'omitnan')...
                                - mean(cat(2,curDataCat{idxUse{1}}),2,'omitnan');
                        end

                        % store day group average for current map type
                        mapData.(catFieldLabels{flc}){ss,mt,gg} = storeData;
                    end

                    % calculate correlation
                    [corrR,corrP] = corr(corrCats{globalPlotCorr(gf)}',...
                        catMean(corrCats{globalPlotCorr(gf)})');
                    mapData.corrR(ss,mt,gg) = corrR;
                    mapData.corrP(ss,mt,gg) = corrP;

                    % store timecourse
                    mapData.timecourse{ss,mt,gg} = curDataCat;

                end
            end
        end


        %% Complete time course plots

        % perform ANOVA
        for mt = 1:nMorph
            for ss = 1:nSex
                % set statistic data
                statData = squeeze(mapData.timecourse(ss,mt,:));
                curAnovaIdx = anovaCats{globalAnovaCat(gf)};

                for aa = 1:length(curAnovaIdx)
                    if all(cellfun(@isempty,statData)); continue; end
                    statDataCat = cellfun(@(x) cat(2,x{curAnovaIdx{aa}}),...
                        statData,'UniformOutput',false);

                    for gg = 1:nGeno
                        statDataCat{gg}(all(isnan(statDataCat{gg}),2),:) = [];
                    end

                    % perform anova
                    if ~any(cellfun(@isempty,statDataCat))
                        if length(curAnovaIdx{aa})>1
                            [pACur,~] = anovaRM2W_full_BH(statDataCat{1},statDataCat{2});
                            pAnova = pACur(1);
                        else
                            [~,pAnova] = ttest2(statDataCat{1},statDataCat{2});
                        end
                    else
                        pAnova = NaN;
                    end

                    % disp([globalTitles{gf} ', ' sexIDs{ss} ', ' morphTypes{mt}...
                    %     ', ' ct{:} ': p = ' num2str(pAnova,2)])

                    % set current axis
                    if ss~=3 && ~isempty(ax{ss,mt})
                        axes(ax{ss,mt}); hold on
                    end

                    % add significance bars
                    plotLineP(curAnovaIdx{aa},pAnova,[],[],[],[],1,sexIDs{ss})
                end
            end
        end

        % store ylimits
        ylimsCur = nan(nSex,nMorph,2);
        for mt = 1:nMorph
            for ss = 1:nSex
                % set current axis
                if ss~=3 && ~isempty(ax{ss,mt})
                    axes(ax{ss,mt}); hold on
                    ylimsCur(ss,mt,:) = ylim;
                end
            end
        end

        % find y limit extremes
        outerIdx = [1 2 4 7];
        innerIdx = [5 8];
        ylimsMinOuter = min(ylimsCur(outerIdx),[],'all','omitnan');
        ylimsMaxOuter = max(ylimsCur(outerIdx+9),[],'all','omitnan');
        ylimsMinInner = min(ylimsCur(innerIdx),[],'all','omitnan');
        ylimsMaxInner = max(ylimsCur(innerIdx+9),[],'all','omitnan');

        % loop through axes
        for mt = 1:nMorph
            for ss = 1:nSex
                % set current axis
                if ss==3 || isempty(ax{ss,mt})
                    continue
                end
                axes(ax{ss,mt}); hold on

                % set y limits
                if ss==2 && mt~=1
                    ylim([ylimsMinInner ylimsMaxInner])
                else
                    ylim([ylimsMinOuter ylimsMaxOuter])
                end

                % set plot labels
                ax{ss,mt}.XTick = 1:length(dayNames);
                ax{ss,mt}.XTickLabelRotation = 0;
                ax{ss,mt}.XTickLabel = dayNames;
                xlabel('learning day')
                ylabel(globalYLabs{gf})
                legend(genotypes,'Location','best')
            end
        end

        % save figure
        savefig([svFile svLabelCur '_timecourse'])

        % save data
        save([svFile '/data' svLabelCur '.mat'],'mapData')


        %% Generate additional plots

        for flc = 1:length(catFieldLabels)
            % skip extra heatmaps
            if ~ismember(flc,globalCats{gf}); continue; end

            curMap = mapData.(catFieldLabels{flc});

            %% Generate bar plot

            % initialize figure
            figure;
            t = tiledlayout(3,3);
            sgtitle([ttlCur ': ' catTtls{flc}])
            ax = cell(nSex,nMorph);

            % loop through axes
            for ss = 1:nSex
                for mt = 1:nMorph
                    % set current axis
                    if ss==1
                        % define axis
                        ax{ss,mt} = nexttile([1 1]); hold on
                        title([sexIDs{ss} '-' morphTypes{mt}])

                        % select data
                        curBarData = squeeze(curMap(ss,mt,:))';

                        % define labels
                        curLabels = genotypes;
                        pShow = [1 2];
                    elseif ss==2
                        % define axis
                        ax{2,mt} = nexttile([2 1]); hold on
                        title([morphTypes{mt} ' by sex'])

                        % select data
                        curBarData = squeeze(curMap(2:3,mt,:));

                        % define labels
                        curLabels = sexIDs(2:3);
                        pShow = [1 2; 1 3; 2 4; 3 4; 2 3];
                    else
                        continue
                    end

                    if ~all(cellfun(@isempty,curBarData))
                        barGroup(curLabels,curBarData,'violin',colors2,pShow)
                    end

                end
            end

            % save figure
            savefig([svFile svLabelCur '_bar_' catSvLabels{flc}])


            %% Plot heat maps

            if any(cellfun(@isempty,curMap),'all')
                continue
            end

            % make heat map
            figure;
            h = heatmapGridSplit(curMap,['Sex',sexIDs],['Morphology',morphTypes],...
                globalHeatSign(gf),{2:3,2:3},{1,2:3});
            sgtitle([ttlCur ': ' catTtls{flc}],'FontSize',18)

            % save figure
            savefig([svFile svLabelCur '_' catSvLabels{flc}])
        end

        if clAll
            close all
        end
    end

end

