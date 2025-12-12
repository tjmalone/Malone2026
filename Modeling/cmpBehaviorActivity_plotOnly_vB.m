%% Compare behavior and activity


%% Initialize figure parameters

clear; close all; clc

p1 = 'D:\AD_Project\Behavior';
cd(p1)

% load data
load('data/model_vB/corrData_vB.mat')

nFOV = size(corrData,2);
nDays = size(corrData,2);

% load genotypes {[WT],[AD]}
load('groupIDs.mat')
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

close all

% set correlation parameters
params = [];
useB = [15];
useA = [19 24];  % all behavior and combined variables
nB = length(useB);
nA = length(useA);
params(1,:) = repmat(useB,1,nA);
params(2,:) = repelem(useA,1,nB);

% close all figures
clAll = 1;
% curmax = 0;
yMax = 4.5;

% whether to average across days
avgON = [1,2];
avgLabels = {' (by session)',' (by mouse)'};
avgSaves = {'session','mouse'};

% day selections
useDays = {2:11,2:11};
% useDays = {2:11,2:11};

corrType = 'Pearson';
% corrType = 'Kendall';
% corrType = 'Spearman';

% loop through correlation pairs
for ff = 1:size(params,2)
    % loop through average types
    for gg = 1:length(useDays)
        %% Generate single figure for all subtypes

        % initialize figure
        figure
        xLab = corrLabels{params(1,ff)};
        yLab = corrLabels{params(2,ff)};
        t = tiledlayout(1,nSexes);
        sgtitle(['Correlation: ' xLab ' vs. ' yLab avgLabels{avgON(gg)}])

        % initialize save data
        dataSave = cell(1,nSexes);

        for ss = 1:nSexes
            % define current data
            data = corrData(:,useDays{gg},params(:,ff));
            useDaysCur = 1:length(useDays{gg});

            % normalize data
            nData = data(:,:,2);
            nData = normalize(nData(:));
            data(:,:,2) = reshape(nData(:),size(data(:,:,2)));

            % curmax = max(curmax,max(abs(data(:,:,2)),[],'all'));

            % define final use data
            if avgON(gg)==1
                useData = data;
            else
                useData = mean(data(:,useDaysCur,:),2,'omitnan');
            end

            % separate data by genotype
            dataFinal = cell(1,nGroups);
            for ii = 1:nGroups
                curGroup = intersect(groups{ii},sexes{ss});
                df = useData(curGroup,:,:);
                curData = squeeze(reshape(df,[],1,2));
                curData(any(isnan(curData),2),:) = [];
                dataFinal{ii} = curData;
            end

            % calculate linear correlation
            dataComb = cat(1,dataFinal{:});
            [rCorr,pCorr] = corr(dataComb(:,1),dataComb(:,2),'type',corrType);

            % plot scatter
            for ii = 1:nGroups
                xCur = dataFinal{ii}(:,1);
                yCur = dataFinal{ii}(:,2);

                if ss~=1
                    nexttile(ss); hold on
                    scatter(xCur,yCur,25,colors6{ss,ii},'LineWidth',0.75,'MarkerEdgeAlpha',0.5);

                    nexttile(1); hold on
                    scatter(xCur,yCur,25,colors6{ss,ii},'filled','LineWidth',0.75,'MarkerEdgeAlpha',1);
                else
                    nexttile(1); hold on
                    scatter(xCur,yCur,1,'MarkerEdgeAlpha',0);
                end
            end

            nexttile(ss); hold on
            if ss==1
                xLims = xlim;
                yLims = ylim;
                yLims = [-yMax yMax];
                ylim(yLims)
            else
                xlim(xLims)
                ylim(yLims)
            end

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
                [rCur,pCur] = corr(xCur,yCur,'type',corrType);
                labelLeg{ii} = [legs{ii} ': r=' num2str(rCur,2) ', p=' num2str(pCur,2)];
            end

            legend(h,labelLeg,'Location','best','AutoUpdate','off')
            xlabel(xLab)
            ylabel(yLab)
            set(gca,'FontSize',12)
            title([sexIDs{ss} avgLabels{avgON(gg)}...
                ': r=' num2str(rCorr,2) ', p=' num2str(pCorr,2)])
            % axis('square')

            % store data for saving
            dataSave{ss} = dataFinal;
        end

        % save figure
        savefig([svFile '/Corr_' corrLabels{params(1,ff)} '-'...
            corrLabels{params(2,ff)} '_by' avgSaves{avgON(gg)} '.fig'])

        % save data
        save([svFile '/data/Corr_' corrLabels{params(1,ff)} '-'...
            corrLabels{params(2,ff)} '_by' avgSaves{avgON(gg)} '.mat'],'corrLabels','dataSave')

        % close figure
        if clAll
            close
        end
    end
end

