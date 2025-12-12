%% distributionsPlot_master.m
% Plot neuronal activity distributions as a function of track location and
% quantification of distribution characteristics.
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
svFile = [p1 '\Figures\Figures' char(datetime('today','format','yyyyMMdd')) '\distributionsMaster'];
if ~isfolder(svFile)
    mkdir(svFile);
    mkdir([svFile '\data']);
end


%% Set input parameters

% distribution key: learning type (lt), run type (rp), distribution type,
% (dt), success/fail (ds)
% learningTypes = {'learn','recall'};
learningTypes = {'learn'};
% runTypes = {'allRuns','trueSF','condSF'};
runTypes = {'allRuns'};
% distTypes = {'dfof','field','mRBR','nRBR','oRBR','oEMD'};
% distTypes = {'dfof','field','nRBR','oRBR'};
distTypes = {'dfof','field','fieldAmp'};
sfTypes = {'Success','Fail'};

distTypes = {'fieldAmp'};

% cell selection key: sex (ss), cell type (ct), morphology type (mt),
% genotype (gg)
sexIDs = {'allSex','female','male'};
% cellTypes = {'allType','common','grid','cue'};
cellTypes = {'common'};
% cellTypes = {'common','grid','cue','other'};
% morphTypes = {'allMorph';'ste';'pyr'};
morphTypes = {'allMorph'};
genotypes = {'WT','AD'};
nGeno = length(genotypes);

% quantification field names
quantFields  = {'CueIn','CueOut','RewPre','RewIn','RewPost','RdgBkg','CorrNoLag','LagToPeakCorr'};
quantGroups = {1:2,3:5,6:8};
nQGroups = length(quantGroups);
nQSubs = cellfun(@length,quantGroups);
qGNames = {'Cue Quantification','Reward Quantification','Cue Template Correlation'};

% define averaging methods (by cell/by FOV)
% avgMethods = {'cell','fov'};
avgMethods = {'cell'};

% set base plot labels
ylabBase = {'\DeltaF/F','Percent of cells with field','RBR correlation',...
    'RBR correlation','RBR correlation','EMD'};
ttlBase = {'\DeltaF/F Distribution','Field Distribution',...
    'mRBR Consistiency Distribution','nRBR Consistiency Distribution',...
    'oRBR Consistiency Distribution','Earth Movers Distance Distribution'};
xBase = {linspace(2.5,597.5,120),linspace(2.5,597.5,120),...
    linspace(12.5,587.5,116),linspace(12.5,587.5,116),...
    linspace(12.5,587.5,116),linspace(12.5,587.5,116)};

xBase = {linspace(2.5,597.5,120)};

% initialize day category matrix
daysDef = struct();

% define learning day categories and labels
daysDef.learn.dayCats = {1,2,8:11};
daysDef.learn.ttlsCats = {'FE','NE Pre-Learning','NE Post-Learning'};
daysDef.learn.svLabelCats = {'FE','NE-Pre','NE-Post'};
daysDef.learn.dayNames = ['FE' string(1:10)];

% define recall day categories and labels
daysDef.recall.dayCats = {1,2};
daysDef.recall.ttlsCats = {'FE-Ref','FE-Recall'};
daysDef.recall.svLabelCats = {'FE-Ref','FE-Recall'};
daysDef.recall.dayNames = ["FE-ref", "FE-recall"];

% define colors
colors2 = {[0 0 1],[1 0 0]};
colors4 = {[0 0 1],[1 0 0];[0.5 0.5 1],[1 0.5 0.5]};
cueTypes = {'Left','Right'};

% whether to replot figures\
replotOn = 1;

% close all figures
clAll = 1;


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

    % loop through run types
    for rt = runTypes

        % loop through distribution types
        for dt = 1:length(distTypes)
            % get current data
            curDistBase = distActivity.(lt{:}).(rt{:}).(distTypes{dt});
            curDistFields = fieldnames(curDistBase);
            nCurFields = length(curDistFields);

            % set bar colors
            if nCurFields==1
                curColors = colors2;
                curPShow = [1 2];
            elseif nCurFields==2
                curColors = colors4;
                curPShow = [1 2;1 3;2 4;3 4];
            end

            % loop through sexes
            for ss = sexIDs

                % loop through cell types
                for ct = cellTypes
                    % get morpholgoy types
                    curTypeBase = cellSelect.(lt{:}).(ss{:}).(ct{:});
                    morphFields = fieldnames(curTypeBase)';

                    % loop through morphology types
                    for mt = intersect(morphFields,morphTypes)
                        % get current cell selections
                        curCellSub = curTypeBase.(mt{:});

                        % loop through averaging method (by cell or by FOV)
                        for am = avgMethods
                            %% Plot indivivudal figure

                            % check for completion
                            ttlCur = [lt{:} '-' rt{:} '-' distTypes{dt}...
                                '-' ss{:} '-' ct{:} '-' mt{:} '-' am{:}];
                            if replotOn==0 && isfile([svFile '\' ttlCur '.fig'])
                                continue
                            end
                            
                            % initialize tiled figure
                            figure
                            t = tiledlayout(3, 12);
                            sgtitle(ttlCur)
                            legs = cell(nCurFields,nGeno);

                            % initialize figure axes
                            ax = {};
                            for ii = 1:nDayCats
                                ax{ii} = nexttile([1, 12/nDayCats]);
                                title([curTtlsCats{ii} ': ' char(dayLabels{ii})])
                            end
                            endTop = ii;

                            for ii = 1:8
                                ax{endTop+ii} = nexttile([1, 3]);
                                if ii==1
                                    title('NE Post-Pre Learning')
                                end
                            end

                            % initialize quantification struct
                            distQuant = cell(1,nQGroups*2);

                            % loop through distribution fields
                            for ds = 1:nCurFields
                                % get current distribution
                                curDistSub = curDistBase.(curDistFields{ds});

                                % loop through genotypes
                                for gg = 1:nGeno
                                    %% Process current data

                                    % set legend
                                    legs{ds,gg} = [genotypes{gg} '-' curDistFields{ds}];

                                    % get current cell selection
                                    curCellSel = curCellSub.(genotypes{gg});
                                    
                                    % correct empty cells
                                    maxIdx = cellfun(@max,curCellSel,'UniformOutput',false);
                                    emptyCellSel = cellfun(@isempty,maxIdx);
                                    emptyData = cellfun(@isempty,curDistSub);
                                    fixIdx = find(~emptyCellSel & emptyData);
                                    refIdx = find(~emptyData,1);
                                    szXCur = size(curDistSub{refIdx},2);
                                    for ii = 1:length(fixIdx)
                                        curDistSub{fixIdx(ii)} = nan(maxIdx{fixIdx(ii)},szXCur);
                                    end

                                    % define current data
                                    curData = cellfun(@(x,y) x(y,:),curDistSub,...
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
                                    if strcmp(lt,'learn')
                                        % calculate each environment separately
                                        curDataQuantSub = cell(1,2);
                                        for ii = 1:2
                                            cueTemp = sum(cueData(ii).cueTemp,2);
                                            szData = size(curDataCat{1},2);
                                            newSt = round((length(cueTemp)-szData)/2)+1;
                                            cueTemp = cueTemp(newSt:newSt+szData-1);
                                            rewLoc = round(cueData(ii).rewLoc/5);
                                            subData = curDataCat(:,cueData(ii).subIdxs);
                                            curDataQuantSub{ii} = distrQuant(subData,cueTemp,rewLoc);
                                        end

                                        % combine quantifications
                                        curDataQuant = struct();
                                        for ii = quantFields
                                            curDataQuant.(ii{:}) = cat(2,...
                                                curDataQuantSub{1}.(ii{:}),curDataQuantSub{2}.(ii{:}));
                                        end
                                    elseif strcmp(lt,'recall')
                                        cueTemp = sum(cueData(1).cueTemp,2)';
                                        szData = size(curDataCat{1},2);
                                        newSt = round((length(cueTemp)-szData)/2)+1;
                                        cueTemp = cueTemp(newSt:newSt+szData-1);
                                        rewLoc = round(cueData(1).rewLoc/5);
                                        curDataQuant = distrQuant(curDataCat,cueTemp',rewLoc);
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
                                            distData{dc} = mean(cat(1,tempCat{curDayCats{dc}}),1,'omitnan');
                                        else
                                            distData{dc} = mean(cat(3,curDataCat{curDayCats{dc}}),3,'omitnan');
                                        end

                                        % plot distribution
                                        semshade(distData{dc},0.3,colors4{ds,gg},xBase{dt});
                                    end

                                    % set current axis
                                    axes(ax{endTop+1}); hold on

                                    % calculate pre versus post differnce
                                    distDiffData = distData{end}-distData{end-1};

                                    % plot difference distribution
                                    semshade(distDiffData,0.3,colors4{ds,gg},xBase{dt});


                                    %% Store quantification data

                                    for ii = 1:nQGroups
                                        for jj = 1:nQSubs(ii)
                                            curSub = quantFields{quantGroups{ii}(jj)};
                                            curSubData = curDataQuant.(curSub);
                                            for kk = 1:2

                                                % combine days of the same day category
                                                if strcmp(ct,'allType') && strcmp(am,'cell')
                                                    tempCat = cellfun(@(x) mean(x,1,'omitnan'),curSubData);
                                                    curSubMean = mean(tempCat(curDayCats{end-2+kk}),'omitnan');
                                                else
                                                    curSubMean = mean(cat(2,curSubData{curDayCats{end-2+kk}}),2,'omitnan');
                                                end

                                                curSubMean(curSubMean==inf) = NaN;

                                                % find store indices
                                                curIdx1 = nQGroups*(kk-1)+ii;
                                                curIdx2 = nCurFields*(gg-1)+ds;

                                                % store quantification data
                                                distQuant{curIdx1}{jj,curIdx2} = curSubMean;
                                            end
                                        end
                                    end

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
                            for dc = 1:nDayCats+1
                                % set current axis
                                axes(ax{dc}); hold on

                                if dc<=nDayCats
                                    ylim([ylimsMin ylimsMax])
                                else

                                end

                                if strcmp(lt,'learn') && dc>1
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


                            %% Plot bar graphs

                            % ylimits
                            ylimsCur = zeros(nQGroups,2,2);

                            for ii = 1:nQGroups*2
                                % define current axis
                                qMod = mod(ii-1,nQGroups)+1;
                                qRem = ceil(ii/nQGroups);
                                qAx = endTop + 4*(qRem-1)+qMod + 1;
                                axes(ax{qAx}); hold on

                                % set title
                                title([qGNames{qMod} ': ' curTtlsCats{end-2+qRem}])

                                % plot current quantification
                                curDistQuant = distQuant{ii};
                                xBar = quantFields(quantGroups{qMod});
                                pShow = generatePShow(nQSubs(qMod),curPShow);
                                barGroup(xBar,curDistQuant,'violin',curColors,pShow)

                                % store y limits
                                ylimsCur(qMod,qRem,:) = ylim;
                            end

                            % find y limit extremes
                            ylimsMin = min(ylimsCur(:,:,1),[],2);
                            ylimsMax = max(ylimsCur(:,:,2),[],2);

                            for ii = 1:nQGroups*2
                                % define current axis
                                qMod = mod(ii-1,nQGroups)+1;
                                qRem = ceil(ii/nQGroups);
                                qAx = endTop + 4*(qRem-1)+qMod + 1;
                                axes(ax{qAx}); hold on

                               ylim([ylimsMin(qMod),ylimsMax(qMod)])
                            end

                            % set legend
                            axes(ax{endTop+nQGroups+2}); hold on
                            for ii = 1:nCurFields*nGeno
                                plot(0,0,'Color',curColors{ii})
                            end
                            axis('off')
                            legend(reshape(legs,nCurFields*nGeno,1))


                            %% Save current figure

                            % save figure
                            savefig([svFile '\' ttlCur])

                            % % save data
                            % save([svFile '\distributions\' ttlCur '.mat'],'curData')

                            % close figure
                            if clAll
                                close
                            end

                        end
                    end
                end
            end
        end
    end
end

