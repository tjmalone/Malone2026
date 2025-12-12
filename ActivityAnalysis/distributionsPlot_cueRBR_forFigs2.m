%% distributionsPlot_cueRBR.m
% Plot neuronal activity distributions as a function of track location and
% quantify in-cue versus out-cue ratios.
%
% Distribution Types:
%   Learning types - learn
%   Run types  - all runs
%   Distribution types - to others RBR
%
% Data Separations:
%   Sex - allSex, female, male
%   Cell type - common
%   Morphology type - allMorph, ste, pyr
%   Genotype - WT, AD
%


%% Load data

clear; clc; close all

p1 = '/MATLAB Drive/FY2025/imagingData';
cd(p1)

% load learning alignments
load('foldersLearning.mat','alignsLearning')

% load distribution data
load('data/distActivityIntra_dfof.mat','distActivityIntra')
% load('data/distActivityInter.mat','distActivityInter')

% load cell selections
load('data/cellSelect.mat','cellSelect')

% set save folder name
svFile = [p1 '/Figures/Figures' char(datetime('today','format','yyyyMMdd')) '/distributions_dil-cut'];
if ~isfolder(svFile)
    mkdir(svFile);
    mkdir([svFile '/data']);
end


%% Set input parameters

% cell selection key: sex (ss), cell type (ct), morphology type (mt),
% genotype (gg)
sexIDs = {'allSex','female','male'};
% cellTypes = {'allType','common'};
% cellTypes = {'common','grid','cue','other'};
cellTypes = {'common'};
morphTypes = {'allMorph','ste','pyr'};
genotypes = {'WT','AD'};
nGeno = length(genotypes);
nSexes = length(sexIDs);
nMorphs = length(morphTypes);

% quantification field names
quantFields  = {'CueIn','CueOut','RdgBkg'};

% define averaging methods (by cell/by FOV)
% avgMethods = {'cell','fov'};
avgMethods = {'cell'};

% set base plot labels
distActivity = struct();
distActivity.oRBR = distActivityIntra;
% distActivity.inter = distActivityInter;

distTypes = {'oRBR','inter'};
distTypes = {'oRBR'};
% distTypes = {'inter'};
ylabBase = {'RBR correlation','Inter-day correlation'};
ttlBase = {'oRBR Consistiency Distribution','Inter-day Consistency Distribution'};
xBase = linspace(12.5,587.5,116);

% define which global fields use common cell referencing or are
% pre-calculated by FOV. 0 = original data by all cells. 1 = original data
% by common cell. -1 = original data by FOV.
globalCommonRef = [0, 1];

% define which heatmaps to generate for each global field
globalCats = {[1 2 3 4], [2 4]};

% define learning day categories and labels
dayCats = {1,2,8:11,2:11};
ttlsCats = {'FE','NE Pre-Learning','NE Post-Learning','All NE'};
dayNames = ['FE' string(1:10)];
nDayCats = length(dayCats);

% define heatmap days
heatDays = {2:11,2,8:11};
heatTtls = {'All NE','NE Pre-Learning','NE Post-Learning','Learning Difference'};
heatFieldLabels = {'All','Pre','Post','Learning'};
heatUseDays = 2:11;

% define cue calculation parameters
% cueDil = 2;
% cueMaxBin = 100;

cueDil = 2;
cueMaxBin = [];

% define colors
colors2 = {[0 0 1],[1 0 0]};
colors4 = {[102 105 255]/255,[255 102 102]/255;[0 0 175]/255,[175 0 0]/255};
cueTypes = {'Left','Right'};

% whether to replot figures/
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
        cueFolder = '/MATLAB Drive/FY2025/AnalysisCode/PostAnalysis/Cues/6mEnv1/';
        cueData(ii).rewLoc = [240 290];
        cueData(ii).subIdxs = 1;
    else
        cueFolder = '/MATLAB Drive/FY2025/AnalysisCode/PostAnalysis/Cues/6mEnv2/';
        cueData(ii).rewLoc = [510 560];
        cueData(ii).subIdxs = 2:11;
    end

    % load cue templates
    load([cueFolder 'tempL.mat'])
    load([cueFolder 'tempR.mat'])
    cueData(ii).cueTemp = [tempL,tempR];
end


%% Plot distributions

% find day category titles
dayLabels = cell(1,nDayCats);
for dc = 1:nDayCats
    curCat = dayCats{dc};
    dayLabels{dc} = strcat({'Days '},dayNames(curCat(1)),'-',dayNames(curCat(end)));
end

% loop through distribution types
for dt = 1:length(distTypes)
    % get current data
    curDist = distActivity.(distTypes{dt});
    curColors = colors2;
    curPShow = [1 2];

    % loop through averaging method (by cell or by FOV)
    for am = avgMethods

        % loop through cell types
        for ct = cellTypes

            % initialize quantification struct
            dataQuant = struct;
            for ii = quantFields
                dataQuant.(ii{:}) = cell(nSexes,nMorphs,nGeno);
            end

            % loop through sexes
            for ss = 1:nSexes

                % get morpholgoy types
                curTypeBase = cellSelect.learn.(sexIDs{ss}).(ct{:});

                % loop through morphology types
                for mt = 1:nMorphs
                    % get current cell selections
                    curCellSub = curTypeBase.(morphTypes{mt});

                    %% Plot indivivudal figure

                    % check for completion
                    ttlCur = [distTypes{dt} '-' sexIDs{ss} '-'...
                        ct{:} '-' morphTypes{mt} '-' am{:}];
                    if replotOn==0 && isfile([svFile '/' ttlCur '.fig'])
                        continue
                    end

                    % initialize tiled figure
                    figure
                    t = tiledlayout(1,nDayCats);
                    sgtitle(ttlCur)

                    % initialize figure axes
                    ax = {};
                    for ii = 1:nDayCats
                        ax{ii} = nexttile();
                        title([ttlsCats{ii} ': ' char(dayLabels{ii})])
                    end

                    % loop through genotypes
                    for gg = 1:nGeno
                        %% Process current data

                        % get current cell selection
                        curCellSel = curCellSub.(genotypes{gg});
                        if globalCommonRef(dt)==1
                            curCellSel = idx2common(curCellSel,alignsLearning);
                        end

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
                        if strcmp(am,'fov')
                            curDataMean = cellfun(@(x) mean(x,1,'omitnan'),...
                                curData,'UniformOutput',false);
                        elseif strcmp(am,'cell')
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

                        % quantify distribution
                        curDataQuantSub = cell(1,2);
                        for ii = 1:2
                            cueTemp = sum(cueData(ii).cueTemp,2);
                            szData = size(curDataCat{1},2);
                            newSt = round((length(cueTemp)-szData)/2)+1;
                            cueTemp = cueTemp(newSt:newSt+szData-1);
                            rewLoc = round(cueData(ii).rewLoc/5)-newSt;
                            subData = curDataCat(:,cueData(ii).subIdxs);
                            curDataQuantSub{ii} = distrQuant(subData,cueTemp,rewLoc,cueDil,cueMaxBin);
                        end

                        % combine quantifications
                        for ii = quantFields
                            % concatenate run type data
                            dataQuant.(ii{:}){ss,mt,gg} = cat(2,...
                                curDataQuantSub{1}.(ii{:}),curDataQuantSub{2}.(ii{:}));
                        end


                        %% Plot distributions

                        % initialize distribution data
                        distData = cell(1,nDayCats);

                        % loop through day categories
                        for dc = 1:nDayCats
                            % set current axis
                            axes(ax{dc}); hold on

                            % combine days of the same day category
                            if strcmp(ct,'allType') && strcmp(am,'cell')
                                tempCat = cellfun(@(x) mean(x,1,'omitnan'),curDataCat,'UniformOutput',false);
                                distData{dc} = mean(cat(1,tempCat{dayCats{dc}}),1,'omitnan');
                            else
                                distData{dc} = mean(cat(3,curDataCat{dayCats{dc}}),3,'omitnan');
                            end

                            % plot distribution
                            semshade(distData{dc},0.3,colors2{gg},xBase);
                        end
                    end

                    %% Complete distribution plots

                    % ylimits
                    ylimsCur = zeros(nDayCats,2);

                    for dc = 1:nDayCats
                        % store y limits
                        axes(ax{dc}); hold on
                        ylimsCur(dc,:) = ylim;
                    end

                    % find y limit extremes
                    ylimsMin = min(ylimsCur(:,1),[],1);
                    ylimsMax = max(ylimsCur(:,2),[],1);

                    % loop through day categories
                    for dc = 1:nDayCats
                        % set current axis
                        axes(ax{dc}); hold on

                        if dc<=nDayCats
                            ylim([ylimsMin ylimsMax])
                        else

                        end

                        if dc>1
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
                    savefig([svFile '/' ttlCur])

                    % % save data
                    % save([svFile '/distributions/' ttlCur '.mat'],'curData')

                    % close figure
                    if clAll
                        close
                    end
                end
            end

            close all


            %% Calculate time course comparisons

            quantFieldsHeat  = {'CueIn','CueOut'};

            % quantFieldsHeat  = {'RdgBkg','RdgBkg'};

            % define in/out ratio
            % dataQuant.IO = cellfun(@(x,y) cellfun(@(xx,yy) (xx-yy)./(xx+yy),x,y,...
            %     'UniformOutput',false),dataQuant.CueIn,dataQuant.CueOut,'UniformOutput',false);
            % quantFieldsHeat = {'IO','IO'};

            % % define raw in/out ratio
            % CueIn = cellfun(@(x) cellfun(@(xx) xx,x,'UniformOutput',false),...
            %     dataQuant.CueIn,'UniformOutput',false);
            % quantFieldsHeat = {'IO','IO'};
            % dataQuant.IO = cellfun(@(x,y) cellfun(@(xx,yy) log(xx./yy),x,y,...
            %     'UniformOutput',false),dataQuant.CueIn,dataQuant.CueOut,'UniformOutput',false);
            % quantFieldsHeat = {'IO','IO'};

            heatY = {'Correlation','Correlation'};
            nQuantHeat = length(quantFieldsHeat);

            % initialize data struct
            mapData = struct();

            for qq = 1:nQuantHeat
                curQuant = dataQuant.(quantFieldsHeat{qq});

                % initialize data struct
                mapData.(quantFieldsHeat{qq}).timecourse = curQuant;
                for flc = 1:length(heatFieldLabels)
                    mapData.(quantFieldsHeat{qq}).(heatFieldLabels{flc}) = cell(nSexes,nMorphs,nGeno);
                end

                % loop through sexes
                for ss = 1:nSexes

                    % loop through morphology types
                    for mt = 1:nMorphs

                        % loop through genotypes
                        for gg = 1:nGeno
                            %% Plot time courses

                            % define current time course
                            dataQuantSub = curQuant{ss,mt,gg};

                            %% Store data

                            for flc = 1:length(heatFieldLabels)
                                if flc~=4
                                    % crop category indices
                                    curDays = heatDays{flc};

                                    % take day averages
                                    storeDMean = mean(cat(3,dataQuantSub{curDays}),3,'omitnan');
                                else
                                    storeA = mean(cat(3,dataQuantSub{heatDays{2}}),3,'omitnan');
                                    storeB = mean(cat(3,dataQuantSub{heatDays{3}}),3,'omitnan');
                                    storeDMean = (storeB-storeA)/mean(storeA,'omitnan')*100;
                                end
                                    % store averages for current map type
                                    mapData.(quantFieldsHeat{qq}).(heatFieldLabels{flc}){ss,mt,gg} = storeDMean;

                            end
                        end
                    end
                end
            end

            for flc = 1:length(heatFieldLabels)
                for ss = 1:nSexes
                    for mt = 1:nMorphs
                        for gg = 1:nGeno
                            data1 = mapData.(quantFieldsHeat{1}).(heatFieldLabels{flc}){ss,mt,gg};
                            data2 = mapData.(quantFieldsHeat{2}).(heatFieldLabels{flc}){ss,mt,gg};

                            dataRatio = data1./data2;
                            dataRatio(dataRatio<=0) = NaN;
                            IO = log(dataRatio);
                            % calculate and store log ratio
                            mapData.IO.(heatFieldLabels{flc}){ss,mt,gg} = IO;
                        end
                    end
                end
            end


            %% Plot bar graphs

            % define sex/morph categories to plot
            useMorphs = {'allMorphs','ste','pyr'};
            % useMorphs = {'allMorphs'};
            useSexes = {'All','Female','Male'};
            nUMorphs = length(useMorphs);
            nUSexes = length(useSexes);
            smSubs = {1,2:3};

            for mt = 1:nUMorphs

                for flc = 1:length(heatFieldLabels)
                    % skip extra heatmaps
                    if ~ismember(flc,globalCats{dt}); continue; end

                    % get current map
                    curMap = mapData.IO.(heatFieldLabels{flc});
                    dataMapSM = squeeze(curMap(:,mt,:));

                    % intitialize figure
                    figure;
                    tiledlayout(1,length(smSubs))

                    % plot bar graphs
                    for ss = 1:length(smSubs)
                        nexttile(); hold on

                        % generate pshow
                        if length(smSubs{ss})==1
                            pShow = [1 2];
                            curColor = colors2;
                        else
                            pShow = [1 2; 3 4; 1 3; 2 4];
                            curColor = colors4;
                        end

                        curPlotData = dataMapSM(smSubs{ss},:);
                        Mlog = cellfun(@(x) mean(x,'omitnan'),dataMapSM(smSubs{ss},:));
                        SEMlog = cellfun(@(x) nansem(x,1),dataMapSM(smSubs{ss},:));
                        Morig = exp(Mlog);
                        SEMorig = Morig.*SEMlog;
                        X = useSexes(smSubs{ss});
                        
                        barGroup(X,curPlotData,'bar',curColor,pShow,[],[],[],Morig,SEMorig);

                        % set labels
                        axis('square')
                        ylabel(ylabBase{dt})
                    end

                    % set labels
                    ttlCur = [distTypes{dt} '-' ct{:} '-' am{:}];
                    sgtitle([ttlCur ': ' heatTtls{flc} ' (' useMorphs{mt} ')'],'FontSize',18)

                    % save figure
                    savefig([svFile '/cueBars_' ttlCur '_' useMorphs{mt} '_' heatFieldLabels{flc}])
                end
            end

            if clAll
                close all
            end
        end
    end
end
