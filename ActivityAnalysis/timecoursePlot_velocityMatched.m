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
globalFields = {'RBR','inter','dfofRdgBkg','ampDiffCueNL','ampDiffAntiNL','ampDiffPerNL',...
    'FWidthNum','FWidthSum','FWidthMean','spatialSelectivity','ampDiffIOraw','ampDiffIOnorm',...
    'decodeIORaw','decodeIONorm','dfofSig','dfofIO','speedScore'};
globalTitles = {'Intra-day consistency (to others)','Inter-day consistency',...
    'Ridge/Background Ratio of \DeltaF/F','Cue Amplitude Difference Mean',...
    'Anti-Cue Amplitude Difference Mean','Cue Amplitude Difference Percentile (of anti-cue)',...
    'Number of Fields','Total Field Coverage','Mean Field Width','Spatial Selectivity',...
    'Amplitude Difference IO ratio (raw)','Amplitude Difference IO ratio (norm)',...
    'Decoding IO ratio (raw)','Decoding IO ratio (norm)','Mean \DeltaF/F_{sig} (while moving)',...
    'In/Out Ratio of \DeltaF/F','Speed Score'};
globalYLabs = {'RBR correlation','inter-day correlation','rdg/bkg ratio',...
    'amplitude difference','amplitude difference','percentile',...
    '# of Fields','Field coverage (cm)','Field width (cm)','Spatial selectivity',...
    'IO ratio','IO ratio','IO ratio','IO ratio','Mean \DeltaF/F','IO ratio','speed score'};
nGlobalFields = length(globalFields);

% define which heatmaps to generate for each global field
globalCats = {[1 3 4], [1 3 4], 1, 1, 1, [1 4], 1, 1, 1, 1, 1, 1, [], [], 1, 1, 1};

% define which global fields use common cell referencing or are
% pre-calculated by FOV. 0 = original data by all cells. 1 = original data
% by common cel. -1 = original data by FOV.
globalCommonRef = [0, 1, -1, 1, 1, 1, 0, 0, 0, 0, 1, 1, -1, -1, 0, -1, 0];

% define whether a high value is expected to align with good or poor
% behavior. 1 = low AD implies deficit, -1 = high AD implies decifit
globalHeatSign = [1, 1, -1, 1, 1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, 1];

% define whether to plot correlations for time course
globalPlotCorr = [0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1];
corrDays = (2:11)';

% define whether a bar graph should be made for results (1/0 = yes/no)
globalPlotBar = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];

% define anova classes
anovaCats = {8:11,2:10};
globalAnovaCat = [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];

% load intra-day consistency data
load('data\intraDayData.mat','dataActIntra')
globalData.RBR = dataActIntra.full.RBR_sigOthers;
globalData.RBR = cellfun(@(x) x',globalData.RBR,'UniformOutput',false);

% load inter-day consistency data
load('data\interDaySubset.mat','dataOut')
globalData.inter = dataOut;

% % load correlation to cue template data
load('data\corrInfoGlobal.mat','dataActivity')
globalData.dfofRdgBkg = dataActivity.dfofRdgBkg;
globalData.dfofIO = dataActivity.dfofIO;

% load amplitude difference
load('data\cueIdentitySubset.mat','dataOut')
globalData.ampDiffCueNL = dataOut.cue;
globalData.ampDiffAntiNL = dataOut.anti;
globalData.ampDiffPerNL = dataOut.cuePer;
globalData.ampDiffIOraw = dataOut.IOraw;
globalData.ampDiffIOnorm = dataOut.IOnorm;

% load by cell properties
load('data\dataByCell.mat','dataAll')
globalData.speedScore = dataAll.speedScore;
globalData.spatialSelectivity = dataAll.spatialSelectivity;
globalData.FWidthNum = dataAll.FWidthNum;
globalData.FWidthSum = dataAll.FWidthSum;
globalData.FWidthMean = dataAll.FWidthMean;

% load decoding properties
load('data\decodeInfoGlobal.mat','decodingActivity')
globalData.decodeIORaw = decodingActivity.decodeIORaw;
globalData.decodeIONorm = decodingActivity.decodeIONorm;

% load dfof
load('data\dfofData.mat','dataDfof')
globalData.dfofSig = dataDfof;
globalData.dfofSig = cellfun(@(x) x',globalData.dfofSig,'UniformOutput',false);

% load cell selections
load('data\cellSelect.mat','cellSelect')
cellSelect = cellSelect.learn;

% load alignments
load('foldersLearning.mat','alignsLearning','trueDays')

% load velocity match info
load('groupIDs.mat','groups','sexes')
genoIdx = groups;
load('D:\AD_Project\Behavior\data\velocityMatch.mat','validPairsAll')

% set save folder name
svFile = [p1 '\Figures\Figures' char(datetime('today','format','yyyyMMdd')) '\timecourse_velMatch'];
if ~isfolder(svFile)
    mkdir(svFile);
    mkdir([svFile '\data']);
end


%% Set input parameters

% cell selection key: sex (ss), cell type (ct), morphology type (mt),
% genotype (gg)
sexIDs = {'allSex','female','male'};
nSex = length(sexIDs);
cellTypes = {'common','grid','nongrid'};
cellTypes = {'common'};
morphTypes = {'allMorph','ste','pyr'};
nMorph = length(morphTypes);
genotypes = {'WT','AD'};
nGeno = length(genotypes);

% define averaging method (by cell/by FOV)
avgMethods = 'cell';

% define learning day categories and labels
dayNames = ['FE' string(1:10)];

% define colors and patterns
colors2 = {[0 0 1],[1 0 0]};
colors4 = {[102 105 255]/255,[255 102 102]/255;[0 0 175]/255,[175 0 0]/255};
patterns = {'-','--'};

% whether to replot figures
replotOn = 1;

% close all figures
clAll = 0;


%% Plot activity variables

% loop through data types
useFields = 1:nGlobalFields;
useFields = [1,2,3,8,10,13,15];
% useFields = 17;
for gf = useFields
    % get current data
    fieldData = globalData.(globalFields{gf});

    % loop through cell types
    for ct = cellTypes
        %% Plot individual figure

        % overwrite decoding cell type
        if ismember(gf,[13 14]) && strcmp(ct,'common')
            ct = {'allType'};
        end

        % check for completion
        ttlCur = [globalTitles{gf} ': ' ct{:} ' cells (by ' avgMethods ')'];
        svLabelCur = ['\' globalFields{gf} '_' ct{:} '-' avgMethods];
        if replotOn==0 && isfile([svFile '\' svLabelCur '_timecourse.fig'])
            continue
        end

        % initialize tiled figure
        figure
        t = tiledlayout(3, 3);
        sgtitle(ttlCur)

        % initialize figure axes
        ax = cell(nSex,nMorph);
        mapData = struct();
        mapData.timecourse = cell(nSex,nMorph,nGeno);

        % initialize correlation info
        if globalPlotCorr(gf)==1
            mapData.corrR = zeros(nSex,nMorph,nGeno);
            mapData.corrP = zeros(nSex,nMorph,nGeno);
        end

        % loop through sexes
        for ss = 1:nSex

            % loop through morphology types
            for mt = 1:nMorph

                % get current cell selections (skip invalid combinations)
                try
                    curCellSub = cellSelect.(sexIDs{ss}).(ct{:}).(morphTypes{mt});
                    if ismember(gf,[13 14]) && mt~=1
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
                        curGroupFull = intersect(genoIdx{gg},sexes{ss});
                        curMacthMice = validPairsAll{ss,dd}(:,gg);
                        curMacthFOV = [curMacthMice*2-1;curMacthMice*2];
                        curUseFOV = curGroupFull(curMacthFOV);
                        curDataCat{dd} = cat(1,curDataMean{curUseFOV,dd});
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

                    % calculate correlation
                    if globalPlotCorr(gf)==1
                        [corrR,corrP] = corr(corrDays,catMean(corrDays)');
                        mapData.corrR(ss,mt,gg) = corrR;
                        mapData.corrP(ss,mt,gg) = corrP;
                    end

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
                if all(cellfun(@isempty,statData)); continue; end
                
                statDataSub = cellfun(@(x) x(anovaCats{globalAnovaCat(gf)}),...
                    statData,'UniformOutput',false);
                statDataCat = cat(1,statDataSub{:});

                % perform anova
                if ~any(cellfun(@isempty,statDataCat))
                    [pACur,~] = anovaO2W_BH(statDataCat);
                    pAnova = pACur(1);
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
                plotLineP(anovaCats{globalAnovaCat(gf)},pAnova,[],[],[],[],1,sexIDs{ss})
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
        save([svFile '\data' svLabelCur '.mat'],'mapData')

        if clAll
            close all
        end
    end

end

