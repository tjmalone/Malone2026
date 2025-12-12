%% trajectoryPlot_mouseCorr.m
% Calculate neuronal activity distribution as a function of track location
% per mouse for a given set of days. Caluclate and analyze intermouse map
% correlation split by sex and genotype
%
% Distribution Types:
%   Learning types - learn
%   Run types  - all runs
%   Distribution types - to others RBR
%
% Data Separations:
%   Sex - allSex, female, male
%   Cell type - all, common
%   Morphology type - allMorph, ste, pyr
%   Genotype - WT, AD
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

% distribution key: learning type (lt), run type (rp), distribution type,
% (dt), success/fail (ds)
learningTypes = {'learn'};
distTypes = {'dfof'};

% cell selection key: sex (ss), cell type (ct), morphology type (mt),
% genotype (gg)
sexIDs = {'female','male'};
% cellTypes = {'allType','common'};
% cellTypes = {'common'};
cellTypes = {'common','grid','nongrid'};
morphTypes = {'allMorph','ste','pyr'};
% morphTypes = {'pyr'};

genotypes = {'WT','AD'};
nGeno = length(genotypes);
nSexes = length(sexIDs);
nMorphs = length(morphTypes);

% define averaging parameters
fovPerMouse = 2;

% define day category and labels
dayCatIdxs = {2,8:11};
dayCats = {'Pre','Post'};
dayNames = ['FE' string(1:10)];

% define colors
colorsBase = [0.85, 0.35, 0.1;0.0, 0.45, 0.55];
col1 = [1 1 1];
colorSpread = [colorsBase(1,:); mean([colorsBase(1,:);col1],1); col1; mean([colorsBase(2,:);col1],1); colorsBase(2,:)];
colorSpreadH = colorSpread(1:3,:);

% close all figures
clAll = 1;


%% Plot distributions

% initialize example data
exData = cell(nGeno,nSexes);

% initialize output statistics
outStats = {};
testName = 'two-tailed unpaired Students t-test';
nUnits = 'mouse pairs';
pair = 0;
MCavg = 0;
limitP = 0;

% define pairwise ttest settings
testNamePairwise = 'two-tailed unpaired Students t-test, Bonferroni-Holm correction for mutiple comparsions';
MCPairwsise = 1;
idxsPairwise = [1 2; 3 4; 2 3; 2 4];
labsPairwise = {'Male';'Female';'PS19 Male vs. WT Female';'PS19'};
nMulti = length(idxsPairwise);
useMulti = 1;

% loop through learning types
for lt = learningTypes

    % loop day categories
    for dc = 1:length(dayCats)

        dayLabels = strcat({'Days '},dayNames(dayCatIdxs{dc}(1)),'-',dayNames(dayCatIdxs{dc}(end)));

        % loop through distribution types
        for dt = 1
            % get current data
            curDistBase = distActivity.(lt{:}).allRuns.(distTypes{dt});

            % loop through cell types
            for ct = cellTypes

                % loop through morphology types
                for mt = 1:nMorphs
                    % skip invalid combination
                    if strcmp(ct,'allType') && ~strcmp(morphTypes{mt},'allMorph')
                        continue
                    end

                    % initialize quantification struct
                    dataQuant = cell(nSexes,nGeno);

                    % loop through sexes
                    for ss = 1:nSexes

                        % get morpholgoy types
                        curTypeBase = cellSelect.(lt{:}).(sexIDs{ss}).(ct{:});

                        % get current cell selections
                        curCellSub = curTypeBase.(morphTypes{mt});

                        %% Plot indivivudal figure

                        % get current distribution
                        curDistSub = curDistBase.all;

                        % loop through genotypes
                        for gg = 1:nGeno
                            %% Process current data

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

                            % fill in empty cells
                            emptyAll = cellfun(@isempty,curData);
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
                                    curData{ff,ii} = nan(size(curData{ff,curFull}));
                                end
                            end

                            % fill in remaining empty cells
                            curEmpty = find(cellfun(@(x) size(x,2),curData)==0);
                            for ii = curEmpty'
                                curData{ii} = zeros(0,szXCur);
                            end

                            % take mouse mean
                            curDataMouse = cell(nFOV/fovPerMouse,nDays);
                            for ff = 1:fovPerMouse:nFOV
                                curFOV = ff:ff+fovPerMouse-1;

                                for dd = 1:nDays
                                    curDataMouse{(ff+fovPerMouse-1)/fovPerMouse,dd}...
                                        = cat(1,curData{curFOV,dd});
                                end
                            end

                            % take mouse mean
                            curDataMean = cellfun(@(x) mean(x,1,'omitnan'),...
                                curDataMouse,'UniformOutput',false);

                            % concatenate across mice and remove nan mice
                            curDataCat = zeros(nFOV/fovPerMouse,nDays,120);
                            for ff = 1:nFOV/fovPerMouse
                                for dd = 1:nDays
                                    curDataCat(ff,dd,:) = curDataMean{ff,dd};
                                end
                            end
                            nanMice = all(isnan(curDataCat),[2 3]);
                            curDataCatUse = curDataCat(~nanMice,:,:);

                            if strcmp(ct,'common') && mt==1
                                exData{gg,ss} = curDataMouse;
                            end


                            %% Plot distributions

                            % take mean across days
                            distData = squeeze(mean(curDataCatUse(:,dayCatIdxs{dc},:),2,'omitnan'));
                            dataQuant{ss,gg} = distData;
                        end
                    end

                    %% Calculate inter-mouse correlation

                    % reshape data
                    dataByGroup = reshape(dataQuant,[],1);
                    dataComb = cat(1,dataByGroup{:})';
                    nMouseGroups = cellfun(@(x) size(x,1),dataByGroup);
                    groupIdxs =  [0; cumsum(nMouseGroups)];
                    nCorrGroups = length(nMouseGroups);

                    % calculate labels
                    xVals = mean([groupIdxs(1:end-1),groupIdxs(2:end)],2)+0.5;
                    xLabs = {'WT-F','WT-M','AD-F','AD-M'};

                    % calculate correlations
                    corrAll = corr(dataComb,dataComb);
                    corrAll(eye(size(corrAll,1))==1) = 0;
                    corrAllN = corrAll;
                    corrAllN(eye(size(corrAllN,1))==1) = NaN;

                    figure;
                    set(gca,'FontSize',18)
                    t = tiledlayout(1,2);
                    ax = cell(1,2);
                    sgtitle(['Inter-mouse dfof correlation: ' ct{:} '-'...
                        morphTypes{mt} ' (' char(dayLabels) ')'])

                    % calculate group means
                    corrMeans = zeros(nCorrGroups,nCorrGroups);
                    for ii = 1:nCorrGroups
                        for jj = 1:nCorrGroups
                            
                            curCorrs = corrAllN(groupIdxs(ii)+1:groupIdxs(ii+1),groupIdxs(jj)+1:groupIdxs(jj+1));
                            corrMeans(ii,jj) = mean(curCorrs,'all','omitnan');
                        end
                    end

                    for kk = 1:2
                        % plot current heatmap
                        ax{kk} = nexttile(); hold on

                        if kk==1
                            % plot by mouse correlations
                            imAlpha = ones(size(corrAllN));
                            imAlpha(isnan(corrAllN))=0;
                            imagesc(corrAllN,'AlphaData',imAlpha,'YData',[size(corrAllN,1) 1]);

                            % % plot border lines
                            xline(groupIdxs(1)+0.5,'k');
                            yline(groupIdxs(end)-groupIdxs(1)+0.5,'k');
                            for ii = 1:nCorrGroups
                                xline(groupIdxs(ii+1)+0.5,'k');
                                yline(groupIdxs(end)-groupIdxs(ii+1)+0.5,'k');
                            end
                        
                            % set specific labels
                            xlim([0.5 groupIdxs(end)+0.5])
                            ylim([0.5 groupIdxs(end)+0.5])
                            xticks(xVals)
                            yticks(xVals)
                            curCMap = customcolormap([0 0.25 0.5 0.75 1],colorSpread);
                            clim([-1 1])
                        else
                            % plot by group correlations
                            imagesc(corrMeans,'YData',[size(corrMeans,1) 1]);

                            % set specific labels
                            xlim([0.5 nCorrGroups+0.5])
                            ylim([0.5 nCorrGroups+0.5])
                            xticks(1:nCorrGroups)
                            yticks(1:nCorrGroups)

                            % set cLimits
                            curCMap = customcolormap([0 0.25 0.5 0.75 1],colorSpread);
                            cLimit = ceil(10*max(0.5,max(corrMeans(:))))/10;
                            clim([-cLimit cLimit])
                        end

                        % set general labels
                        colorbar
                        ax{kk}.XAxisLocation = 'top';
                        ax{kk}.Colormap = curCMap;
                        axis('square')
                        xticklabels(xLabs)
                        yticklabels(fliplr(xLabs))
                        set(gca,'FontSize',16)
                    end


                    %% Perform t-test

                    combs = combvec(1:4, 1:4)';
                    combs(combs(:,1)>combs(:,2),:) = [];
                    nCombs = size(combs,1);

                    pVals = zeros(nCombs);
                    for ii = 1:nCombs
                        set1A = groupIdxs(combs(ii,1))+1:groupIdxs(combs(ii,1)+1);
                        set1B = groupIdxs(combs(ii,2))+1:groupIdxs(combs(ii,2)+1);
                        corrs1 = corrAllN(set1A,set1B);

                        for jj = 1:nCombs
                            set2A = groupIdxs(combs(jj,1))+1:groupIdxs(combs(jj,1)+1);
                            set2B = groupIdxs(combs(jj,2))+1:groupIdxs(combs(jj,2)+1);
                            corrs2 = corrAllN(set2A,set2B);

                            [~,curP] = ttest2(corrs1(:),corrs2(:));

                            pVals(ii,jj) = curP;
                        end
                    end

                    % define significance sets
                    sets = {[2 4 5],[7 8 9];[1 3 6],10};
                    pGroup = zeros(1,size(sets,1));

                    % calculate significance for each set
                    for st = 1:size(sets,1)
                        corrsGroup = cell(1,2);
                        for ii = 1:nCombs
                            setA = groupIdxs(combs(ii,1))+1:groupIdxs(combs(ii,1)+1);
                            setB = groupIdxs(combs(ii,2))+1:groupIdxs(combs(ii,2)+1);
                            corrs = corrAllN(setA,setB);

                            % correct for duplicates
                            if combs(ii,1)==combs(ii,2)
                                corrsTril = corrs(tril(true(size(corrs)),-1));
                            else
                                corrsTril = corrs(:);
                            end

                            for jj= 1:2
                                if ismember(ii,sets{st,jj})
                                    corrsGroup{jj} = [corrsGroup{jj};corrsTril];
                                end
                            end
                        end
                        [~,pGroup(st)] = ttest2(corrsGroup{1},corrsGroup{2});
                    end

                    % set figure title
                    sgtitle(['Inter-mouse dfof correlation: ' ct{:} '-'...
                        morphTypes{mt} ' (' char(dayLabels) '), p_{inter-group} = ' num2str(pGroup(1),2)...
                        ', p_{intra-group}' num2str(pGroup(2),2)],'FontSize',16)

                    % save figure
                    svLabelCur = [ct{:} '-' morphTypes{mt} '_' dayCats{dc}];
                    savefig([svFile '/interMouseCorr_' svLabelCur])

                    % save data
                    save([svFile '/data/interMouseCorr_' svLabelCur '.mat'],'corrAll','groupIdxs')


                    %% Generate between/within group bar graph

                    sets = struct();
                    labels = struct();
                    pShow = struct();
                    colors = struct();
                    plotTypes = {'betweenMP','betweenFP','betweenMW','betweenFW','within'};
                    nPlotGroups = [1 1 1 1 2];
                    nTypes = length(plotTypes);
                    colors2 = {[0.5 0.5 0.5],[175 0 0]/255};
                    colors4 = {[0 0 175]/255,[175 0 0]/255,[102 105 255]/255,[255 102 102]/255};

                    % define between groups
                    sets.betweenMP = {[2 4 5],[7 8 9]};
                    labels.betweenMP = {'Other to Other','PS19 Male to Other'};
                    pShow.betweenMP = [1 2];
                    colors.betweenMP = colors2;

                    % define between groups
                    sets.betweenFP = {[2 7 8],[4 5 9]};
                    labels.betweenFP = {'Other to Other','PS19 Female to Other'};
                    pShow.betweenFP = [1 2];
                    colors.betweenFP = colors2;

                    % define between groups
                    sets.betweenMW = {[4 7 9],[2 5 8]};
                    labels.betweenMW = {'Other to Other','WT Male to Other'};
                    pShow.betweenMW = [1 2];
                    colors.betweenMW = colors2;

                    % define between groups
                    sets.betweenFW = {[5 8 9],[2 4 7]};
                    labels.betweenFW = {'Other to Other','WT Female to Other'};
                    pShow.betweenFW = [1 2];
                    colors.betweenFW = colors2;

                    % define within groups
                    sets.within = {3, 10, 1, 6};
                    labels.within = {'WT Male','PS19 Male','WT Female','PS19 Female'};
                    pShow.within = [1 2;1 3;2 4; 3 4];
                    colors.within = colors4;

                    figure;
                    tiledlayout(1,nTypes)

                    for tt = 1:nTypes
                        % define current data
                        curSets = sets.(plotTypes{tt});

                        % initialize plot data array
                        nBars = length(curSets);
                        plotData = cell(1,nBars);

                        % set bar data
                        for bb = 1:nBars
                            corrsGroup = [];
                            for ii = 1:nCombs
                                setA = groupIdxs(combs(ii,1))+1:groupIdxs(combs(ii,1)+1);
                                setB = groupIdxs(combs(ii,2))+1:groupIdxs(combs(ii,2)+1);
                                corrs = corrAllN(setA,setB);

                                % correct for duplicates
                                if combs(ii,1)==combs(ii,2)
                                    corrsTril = corrs(tril(true(size(corrs)),-1));
                                else
                                    corrsTril = corrs(:);
                                end

                                if ismember(ii,curSets{bb})
                                    corrsGroup = [corrsGroup;corrsTril];
                                end
                            end

                            plotData{bb} = corrsGroup;
                        end

                        % plot data
                        nexttile(tt); hold on
                        barGroup(labels.(plotTypes{tt}),plotData,'violin',...
                            colors.(plotTypes{tt}),pShow.(plotTypes{tt}))

                        % add full statistics output
                        curCat = [ct{:} '-' morphTypes{mt} ': ' plotTypes{tt} ' (' char(dayLabels) ')'];
                        if nPlotGroups(tt)>1
                            statData = reshape(plotData,[],2);
                        else
                            statData = plotData';
                        end

                        outStats(end+1,:) = [curCat ttestEffectSize(...
                            statData(1,:),statData(2,:),testName,nUnits,pair,MCavg,limitP)];


                        if tt==5
                            multiCat = strcat(repmat({[curCat ' multi: ']},nMulti,1),labsPairwise);
                            outStats(end+1:end+nMulti,:) = [multiCat ttestAllPairsEffectSize(...
                                plotData,testNamePairwise,nUnits,pair,MCPairwsise,idxsPairwise)];
                        end
                    end

                    sgtitle(['Inter-mouse dfof correlation: ' ct{:} '-'...
                        morphTypes{mt} ' (' char(dayLabels) ')'],'FontSize',16)
                    savefig([svFile '/interMouseCorrViolin_' svLabelCur])

                end
            end

            close all
        end
    end
end


%% Plot example distribtuions

close all

% define example FOV
exMouse = 14;
% exMouse = 17;
exDays = 8:11;
nMice = length(exMouse);
nDays = length(exDays);

% load and process example data
exDataUse = exData{2,1}(exMouse,exDays);
exDataMean = cellfun(@(x) mean(x,1,'omitnan'),exDataUse,'UniformOutput',false);

% find limits
valMax = max(cellfun(@max,exDataMean),[],'all');
valMin = min(cellfun(@min,exDataMean),[],'all');

% plot activity
figure; hold on
tiledlayout(nMice,nDays)

for mm = 1:nMice
    for dd = 1:nDays
        nexttile(); hold on
        plot(exDataMean{mm,dd})
        ylim([valMin-0.01 valMax+0.01])
        axis('square')
    end
end

% plot examples
A = cat(1,exDataMean{:});
AA = mean(A,1,'omitnan');
figure; hold on
plot(AA/mean(AA))

BB = dataComb(:,[9 13]);
BB = BB./mean(BB,1);
% BB = dataComb(:,14);
plot(BB)

corr(dataComb(:,[9 13]))


%% Plot full group overlays

% calculate all mouse means
exDays = 8:11;

% define cue/reward info
colorsCue = [0 0 0;0.5 0.5 0.5];
colorsRew = [0.75 1 1];
cueX = 2.5:5:597.5;

% define cue template
cueData = struct();
cueFolder = 'D:/AnalysisCode/PostAnalysis/Cues/6mEnv2/';
rewLoc = [510 560];

% load cue templates
load([cueFolder 'tempL.mat'])
load([cueFolder 'tempR.mat'])
cueTemp = [tempL,tempR];

colors4 = {[102 105 255]/255,[255 102 102]/255;[0 0 175]/255,[175 0 0]/255};

figure; hold on
tiledlayout(2,2)

exDataMean = cell(nSexes,nGeno);
for gg = 1:nGeno
    for ss = 1:nSexes
        % define current data
        curData = exData{gg,ss}(:,exDays);

        % calculate averages
        for ff = 1:nFOV/2
            curDataCat = cat(3,curData{ff,:});
            curDataMean = mean(curDataCat,3,'omitnan');
            curDataMouse = mean(curDataMean,1,'omitnan');

            if ~all(isnan(curDataMouse))
                exDataMean{ss,gg}(end+1,:) = curDataMouse;
            end
        end

        
        % curPlot = exDataMean{ss,gg}./mean(exDataMean{ss,gg},2,'omitnan');

        % curPlot = normalize(exDataMean{ss,gg},2);

        valMin = min(exDataMean{ss,gg},[],2);
        valMax = max(exDataMean{ss,gg},[],2);
        curPlot = (exDataMean{ss,gg}-valMin)./(valMax-valMin);

        nexttile(); hold on

        % plot mean
        meanPlotNorm = mean(curPlot,1,'omitnan');
        % meanPlotNorm = (meanPlot-min(meanPlot))/(max(meanPlot)-min(meanPlot));
        plot(cueX,meanPlotNorm,'Color',[0 0 0],'LineWidth',1.5)

        % plot overlays
        plot(cueX,curPlot','Color',colors4{ss,gg})
        % 
        % % plot cues and reward
        % for ii = 1:size(cueTemp,2)
        %     plotCues(cueX,cueTemp(:,ii),max(ylim),colorsCue(ii,:),min(ylim));
        % end
        % 
        % % plot rewards
        % pos = [rewLoc(1) min(ylim) diff(rewLoc) diff(ylim)];
        % rectangle('Position',pos,'FaceColor',colorsRew,...
        %     'EdgeColor','none')

        % flip children order to plot reward behind distributions
        axCur = gca;
        axCur.Children = flipud(axCur.Children);
    end
end


%% Plot mean group overlay with SEM

% calculate all mouse means
exDays = 8:11;

% define cue/reward info
colorsCue = [0 0 0;0.5 0.5 0.5];
colorsRew = [0.75 1 1];
cueX = 2.5:5:597.5;

% define cue template
cueData = struct();
cueFolder = 'D:/AnalysisCode/PostAnalysis/Cues/6mEnv2/';
rewLoc = [510 560];

% load cue templates
load([cueFolder 'tempL.mat'])
load([cueFolder 'tempR.mat'])
cueTemp = [tempL,tempR];

colors4 = {[102 105 255]/255,[255 102 102]/255;[0 0 175]/255,[175 0 0]/255};

figure; hold on
tiledlayout(2,2)

exDataMean = cell(nSexes,nGeno);
for gg = 1:nGeno
    for ss = 1:nSexes
        % define current data
        curData = exData{gg,ss}(:,exDays);

        % calculate averages
        for ff = 1:nFOV/2
            curDataCat = cat(3,curData{ff,:});
            curDataMean = mean(curDataCat,3,'omitnan');
            curDataMouse = mean(curDataMean,1,'omitnan');

            if ~all(isnan(curDataMouse))
                exDataMean{ss,gg}(end+1,:) = curDataMouse;
            end
        end

        valMin = min(exDataMean{ss,gg},[],2);
        valMax = max(exDataMean{ss,gg},[],2);
        curPlot = (exDataMean{ss,gg}-valMin)./(valMax-valMin);

        nexttile(); hold on
        semshade(curPlot,0.3,colors4{ss,gg},cueX);
    end
end

