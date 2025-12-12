%% Compare behavior and activity


%% Initialize figure parameters

clear; close all; clc

p1 = 'D:\AD_Project\Behavior';
cd(p1)

% load common cell activity
load('data/model_vB/activityDataVB.mat','activityData')
fullData.common = activityData;

% load grid/non-grid cell activity
load('data/model_vB/activityData_GnG_vB.mat','activityData')
fullData.grid = activityData.grid;
fullData.nongrid = activityData.nongrid;

% load stellate/pyramidal cell activity
load('data/model_vB/activityData_PS_vB.mat','activityData')
fullData.ste = activityData.ste;
fullData.pyr = activityData.pyr;

% define cell types
cellTypes = fieldnames(fullData);
cellTypesN = length(cellTypes);

% define fields
baseFields = {'speedScoreAbs'};
cmpFields = {'dfofIORaw','decodeIORaw'};
baseN = length(baseFields);
cmpN = length(cmpFields);

% define skipped comparisons (1 = skip)
cmpSkip = {[], [], [], [2], [2]};

% define axes limits
baseX = [0 0.15; 0 0.2; 0 0.2; 0 0.3; 0 0.3];
cmpY = {[0.6 1.8],[0.9 1.6];...
    [0 3.5],[0.8 1.8];...
    [0 3.5],[0.8 1.8];...
    [0 10],[];...
    [0 10],[]};

% day selections
useDays = 2:11;

% load genotypes {[WT],[AD]}
load('D:\AD_Project\imagingData\groupIDs.mat')
nGroups = length(groups);
nSexes = length(sexes);
sexIDs = {'allSex','female','male'};

% define legend
legs = [groupIDs 'All'];

% define plot colors
colors = {[0 0 1],[1 0 0]};
colors6 = {[0 0 1],[1 0 0];[102 105 255]/255,[255 102 102]/255;[0 0 175]/255,[175 0 0]/255};

% set save folder name
svFile = [p1 '/Figures/Figures' char(datetime('today','format','yyyyMMdd')) '/activityCorr_vB'];
if ~isfile(svFile)
    mkdir(svFile);
    mkdir([svFile '/data']);
end


%% Plot correlation pairs

% set figure close
clAll = 0;

% initialize statistics
outStats = {};
testNameCorr = 'two-tailed Pearson linear correlation';

% loop through cell types
for tt = 1:cellTypesN

    % loop through base fields
    for bb = 1:baseN
        %% Prepare data

        % extract base data
        baseData = {fullData.(cellTypes{tt})(:).(baseFields{bb})};

        % take base data mean
        baseDataMean = cellfun(@(x)...
            cellfun(@(y) mean(y,'omitnan'),x)...
            ,baseData,'UniformOutput',false);
        baseDataCat = cat(1,baseDataMean{:});


        %% Plot correlations

        % loop through comparison fields
        for cc = 1:cmpN
            % skip comparisons
            if ismember(cc,cmpSkip{tt})
                continue
            end

            % extract comparison data
            cmpData = {fullData.(cellTypes{tt})(:).(cmpFields{cc})};

            % take comparison data mean
            cmpDataMean = cellfun(@(x)...
                cellfun(@(y) mean(y,'omitnan'),x)...
                ,cmpData,'UniformOutput',false);
            cmpDataCat = cat(1,cmpDataMean{:});

            % define final correlation data
            plotData = cat(3,baseDataCat(:,useDays),cmpDataCat(:,useDays));

            % % normalize data
            % for dd = 1:2
            %     % normalize data
            %     nData = plotData(:,:,dd);
            %     nData = normalize(nData(:));
            %     plotData(:,:,dd) = reshape(nData(:),size(plotData(:,:,dd)));
            % end

            %% Generate single figure for all subtypes

            % initialize figure
            figure
            xLab = baseFields{bb};
            yLab = cmpFields{cc};
            t = tiledlayout(1,nSexes);
            sgtitle(['Correlation (' cellTypes{tt} ' cells): ' xLab ' vs. ' yLab])

            % initialize save data
            dataSave = cell(1,nSexes);

            for ss = 1:nSexes
                % separate data by genotype
                dataFinal = cell(1,nGroups);
                for ii = 1:nGroups
                    curGroup = intersect(groups{ii},sexes{ss});
                    df = plotData(curGroup,:,:);
                    curData = squeeze(reshape(df,[],1,2));
                    curData(any(isnan(curData),2),:) = [];
                    dataFinal{ii} = curData;
                end

                % calculate linear correlation
                dataComb = cat(1,dataFinal{:});
                [rCorr,pCorr] = corr(dataComb(:,1),dataComb(:,2));

                % plot scatter
                for ii = 1:nGroups
                    xCur = dataFinal{ii}(:,1);
                    yCur = dataFinal{ii}(:,2);

                    if ss~=1
                        nexttile(ss); hold on
                        scatter(xCur,yCur,25,colors6{ss,ii},'filled','LineWidth',0.75,'MarkerEdgeAlpha',1);

                        nexttile(1); hold on
                        scatter(xCur,yCur,25,colors6{ss,ii},'filled','LineWidth',0.75,'MarkerEdgeAlpha',1);
                    else
                        nexttile(1); hold on
                        scatter(xCur,yCur,1,'MarkerEdgeAlpha',0);
                    end
                end

                nexttile(ss); hold on

                % set limits
                xlim(baseX(tt,:))
                ylim(cmpY{tt,cc})

                % initialize labels
                labelLeg = cell(1,3);
                h = zeros(1,nGroups);

                % plot trend lines
                for ii = 1:nGroups+1
                    if ii<=nGroups
                        xCur = dataFinal{ii}(:,1);
                        yCur = dataFinal{ii}(:,2);
                        curColor = colors6{ss,ii};
                    else
                        xCur = dataComb(:,1);
                        yCur = dataComb(:,2);
                        curColor = 'k';
                    end

                    % calculate trendline
                    coefficients = polyfit(xCur,yCur, 1);
                    xFit = [min(xlim) max(xlim)];
                    yFit = polyval(coefficients, xFit);

                    % plot trendline
                    h(ii) = plot(xFit,yFit,'-','Color',curColor,'LineWidth',0.5);

                    % calculate correlation
                    [rCur,pCur] = corr(xCur,yCur);
                    labelLeg{ii} = [legs{ii} ': r=' num2str(rCur,2) ', p=' num2str(pCur,2)];

                    % calculate full statistics
                    corrCat = [cmpFields{cc} ' ' cellTypes{tt} '-' sexIDs{ss} ': ' legs{ii}];
                    outStats(end+1,:) = [corrCat corrEffectSize(xCur,yCur,testNameCorr,'FOV sessions')];
                end

                legend(h,labelLeg,'Location','best','AutoUpdate','off')
                xlabel(xLab)
                ylabel(yLab)
                set(gca,'FontSize',12)
                title(sexIDs{ss})
                % axis('square')

                % store data for saving
                dataSave{ss} = dataFinal;
            end

            % save figure and data
            svCur = ['/Corr_' baseFields{bb} '-'...
                cmpFields{cc} '_' cellTypes{tt} '_byFOV'];
            savefig([svFile svCur '.fig'])
            save([svFile '/data' svCur '.mat'],'dataSave')

            % close figures
            if clAll==1
                close
            end


            %% Plot mouse averaged correlation

            if tt~=1; continue; end

            % initialize figure
            figure
            xLab = baseFields{bb};
            yLab = cmpFields{cc};
            t = tiledlayout(1,nSexes);
            sgtitle(['Mean Correlation (' cellTypes{tt} ' cells): ' xLab ' vs. ' yLab])

            for ss = 1:nSexes
                % separate data by genotype
                dataMean = cell(1,nGroups);
                for ii = 1:nGroups
                    curGroup = intersect(groups{ii},sexes{ss});
                    curData = squeeze(mean(plotData(curGroup,:,:),1,'omitnan'));
                    curData(any(isnan(curData),2),:) = [];
                    dataMean{ii} = curData;
                end

                % calculate linear correlation
                dataComb = cat(1,dataMean{:});
                [rCorr,pCorr] = corr(dataComb(:,1),dataComb(:,2));

                % plot scatter
                for ii = 1:nGroups
                    xCur = dataMean{ii}(:,1);
                    yCur = dataMean{ii}(:,2);

                    if ss~=1
                        nexttile(ss); hold on
                        scatter(xCur,yCur,25,colors6{ss,ii},'filled','LineWidth',0.75,'MarkerEdgeAlpha',1);

                        nexttile(1); hold on
                        scatter(xCur,yCur,25,colors6{ss,ii},'filled','LineWidth',0.75,'MarkerEdgeAlpha',1);
                    else
                        nexttile(1); hold on
                        scatter(xCur,yCur,1,'MarkerEdgeAlpha',0);
                    end
                end

                nexttile(ss); hold on

                % set limits
                % xlim(baseX{bb})
                % ylim(cmpY{cc})

                % initialize labels
                labelLeg = cell(1,3);
                h = zeros(1,nGroups);

                % plot trend lines
                for ii = 1:nGroups+1
                    if ii<=nGroups
                        xCur = dataMean{ii}(:,1);
                        yCur = dataMean{ii}(:,2);
                        curColor = colors6{ss,ii};
                    else
                        xCur = dataComb(:,1);
                        yCur = dataComb(:,2);
                        curColor = 'k';
                    end

                    % calculate trendline
                    coefficients = polyfit(xCur,yCur, 1);
                    xFit = [min(xlim) max(xlim)];
                    yFit = polyval(coefficients, xFit);

                    % plot trendline
                    h(ii) = plot(xFit,yFit,'-','Color',curColor,'LineWidth',0.5);

                    % calculate correlation
                    [rCur,pCur] = corr(xCur,yCur);
                    labelLeg{ii} = [legs{ii} ': r=' num2str(rCur,2) ', p=' num2str(pCur,2)];
                end

                legend(h,labelLeg,'Location','best','AutoUpdate','off')
                xlabel(xLab)
                ylabel(yLab)
                set(gca,'FontSize',12)
                title(sexIDs{ss})
                % axis('square')
            end
        end


        %% Plot time course by FOV

        if tt~=1; continue; end

        figure;
        tiledlayout(1,nSexes)

        for ss = 1:nSexes
            % separate data by genotype
            dataTime = cell(1,nGroups);
            for ii = 1:nGroups
                curGroup = intersect(groups{ii},sexes{ss});
                dataTime{ii} = baseDataCat(curGroup,useDays);
            end

            nexttile(ss); hold on
            for ii = 1:nGroups

                M = mean(dataTime{ii},1,'omitnan');
                S = nansem(dataTime{ii},1);
                errorbar(M,S,'Color',colors{ii})
            end
        end
    end
end




