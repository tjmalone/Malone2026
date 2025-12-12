%% interDaybyLap.m
% generate pairwise correlation calculation for groupings of laps across
% days. Requires output from findImageLaps (post-Analysis).

clear; clc; close all;


%% Set input parameters

p1 = 'D:\AD_Project\imagingData';
cd(p1)

% load alignment data
load('foldersLearning.mat')

% load genotypes {[WT],[AD]}
load('groupIDs.mat')
nGroups = length(groups);

% remove mice with skip days
% groups = {[1	2	3	4	5	6	9	10	23	24	31	32	37	38	41	42],...
%     [7	8	19	20	21	22	25	26	27	28	29	30	33	34	35	36	39	40]};

% remove mouse with bad inter-day
% groups = {[1:6 9:14 23:24	31	32	37	38	41	42],...
%     [7:8 15:16 19:22 25:28	29	30	33	34	35	36	39	40]};

useDays = 1:11;
nDays = length(useDays);
nFOV = length(foldersLearning);

% plot by day
lapGroups = {{2:20};{2:5,6:9,10:12,13:16,17:20};num2cell(2:20)};
% lapGroups = {{2:20}};


%% Find common cell activity for all mice

dataActInter = struct();
dataActInter.dfof = cell(nFOV,nDays);
dataActInter.lapID = cell(nFOV,nDays);

% cycle through FOV
for ii = 1:nFOV
    disp(ii)

    % cycle through days
    for jj = 1:nDays
        cd(foldersLearning{ii}{jj})

        % load activity and lap data
        load('RunByRun_sig\dfofM_sig.mat')
        load('imageLapIdx.mat')

        % extract activity data
        dfofInterp = findDfofInterM(dfofM_sig,imageLapIdx.lapIdx);
        curDfof = dfofInterp(alignsLearning{ii}(:,jj));
        curM = cat(3,curDfof{:});

        % identify true day index
        trueIdx = trueDays(ii,jj);
        if trueIdx>jj
            % set skipped day data to nan
            if trueDays(ii,jj-1)==trueIdx-2
                dataActInter.dfof{ii,jj} = nan(size(curM));
                dataActInter.lapID{ii,jj} = NaN;
            end

            % skip days past day limits
            if trueIdx>max(useDays)
                continue
            end
        end

        % store activity data
        dataActInter.dfof{ii,trueIdx} = curM;
        dataActInter.lapID{ii,trueIdx} = imageLapIdx.lapRun;
    end
end

cd(p1)


%% Generate interday matrices

dataAllInter = struct();

% cycle through plot types
for ff = 1:length(lapGroups)
    %% Generate lap averages
    disp(ff)

    % set current lap groups
    curLapGroup = lapGroups{ff};
    nCurLapGroup = length(curLapGroup);

    % initialize mean activity matrix
    dataAllInter(ff).meanDfof = cell(nFOV,nDays);

    for ii = 1:nFOV
        for jj = 1:nDays
            % get current data
            curDfof = dataActInter.dfof{ii,jj};
            curLaps = dataActInter.lapID{ii,jj};

            % initialize current mean activity cell
            nLaps = size(curDfof,1);
            nBins = size(curDfof,2);
            nCells = size(curDfof,3);
            dataAllInter(ff).meanDfof{ii,jj} = nan(nCurLapGroup,nBins,nCells);

            for kk = 1:nCurLapGroup
                % find laps in corresponding lap group
                useLaps = ismember(curLaps,curLapGroup{kk});

                % calculate mean activity from corresponding laps
                if all(useLaps==0)
                    curMean = nan(1,nBins,nCells);
                else
                    curMean = mean(curDfof(useLaps,:,:),1,'omitnan');
                end

                % store mean activity
                dataAllInter(ff).meanDfof{ii,jj}(kk,:,:) = curMean;
            end
        end
    end


    %% Calculate correlations by group

    nCmps = nCurLapGroup*nDays;
    dataAllInter(ff).meanCombCell = cell(nGroups,nDays);
    dataAllInter(ff).meanCombAll = cell(1,nGroups);
    dataAllInter(ff).corrMMat = cell(1,nGroups);
    dataAllInter(ff).corrMFOV = cell(1,nGroups);
    dataAllInter(ff).corrMInd = cell(1,nGroups);

    for gg = 1:nGroups
        %%

        % concatenate data by lap
        meanCombCell = cell(1,nDays);
        for ii = 1:nDays
            meanCombCell{ii} = cat(3,dataAllInter(ff).meanDfof{groups{gg},ii});
        end

        % concatenate data by day
        meanCombAll = cat(1,meanCombCell{:});

        % initialize correlation matrices
        nCells = size(meanCombAll,3);
        corrMMat = nan(nCmps,nCmps);
        corrMFOV = nan(nCmps,nCmps,length(groups{gg}));
        corrMInd = nan(nCmps,nCmps,nCells);

        % calculate pairwise correlations for combined matrix
        for ii = 1:nCmps
            A = squeeze(meanCombAll(ii,:,:))';
            useA = ~isnan(A(:,1));

            for jj = ii:nCmps
                B = squeeze(meanCombAll(jj,:,:))';
                useB = ~isnan(B(:,1));

                useC = useA&useB;
                corrMMat(ii,jj) = corr2(A(useC,:),B(useC,:));
                corrMMat(jj,ii) = corrMMat(ii,jj);
            end
        end

        % calculate pairwise correlations for FOV
        nSzs = cellfun(@(x) size(x,1),alignsLearning(groups{gg}));
        cumSzs = [0;cumsum(nSzs)];

        for ii = 1:length(groups{gg})
            curCombAll = meanCombAll(:,:,cumSzs(ii)+1:cumSzs(ii+1));

            for jj = 1:nCmps
                A = squeeze(curCombAll(jj,:,:))';
                useA = ~isnan(A(:,1));

                for kk = jj:nCmps
                    B = squeeze(curCombAll(kk,:,:))';
                    useB = ~isnan(B(:,1));

                    useC = useA&useB;
                    corrMFOV(jj,kk,ii) = corr2(A(useC,:),B(useC,:));
                    corrMFOV(kk,jj,ii) = corrMFOV(jj,kk,ii);
                end
            end
        end

        % calculate pairwise correlations per cell
        for ii = 1:nCells
            curCell = meanCombAll(:,:,ii)';
            corrMInd(:,:,ii) = corr(curCell,curCell);
        end

        % store results
        dataAllInter(ff).meanCombCell(gg,:) = meanCombCell;
        dataAllInter(ff).meanCombAll{gg} = meanCombAll;
        dataAllInter(ff).corrMMat{gg} = corrMMat;
        dataAllInter(ff).corrMFOV{gg} = corrMFOV;
        dataAllInter(ff).corrMInd{gg} = corrMInd;

    end
end

save('data\interDayData.mat','dataActInter','dataAllInter','lapGroups','-v7.3')


%% Initialize figure parameters

clear; close all; clc

p1 = 'D:\AD_Project\imagingData';
cd(p1)

% load data
load('data\interDayData.mat','dataAllInter')

% load genotypes {[WT],[AD]}
load('groupIDs.mat')
nGroups = length(groups);

% manually set genotypes
% groups = {[1:6 9:14 23:24],[7:8 15:22 25:28]};
% groups = {[1:6 9:14 23:24],[7:8 15:16 19:22 25:28]}; % exclude mouse with two maps
% groups = {[1	2	3	4	5	6	9	10	11	12	13	14	23	24	31	32	37	38	41	42],...
%     [7	8	15	16	17	18	19	20	21	22	25	26	27	28	29	30	33	34	35	36	39	40]};

% % remove mice with skip days
% groupsFinal = {[1	2	3	4	5	6	9	10	23	24	31	32	37	38	41	42],...
%     [7	8	19	20	21	22	25	26	27	28	29	30	33	34	35	36	39	40]};
% gType = 'noSkipMice';

% keep mice with skip days
% groupsFinal = groups;
% gType = 'allMice';

% female mice
% groupsFinal = cell(1,2);
% for gg = 1:2
%     groupsFinal{gg} = intersect(groups{gg},sexes{2});
% end
% gType = 'female';

% male mice
groupsFinal = cell(1,2);
for gg = 1:2
    groupsFinal{gg} = intersect(groups{gg},sexes{3});
end
gType = 'male';

% remove mouse with bad inter-day
% groups = {[1:6 9:14 23:24	31	32	37	38	41	42],...
%     [7:8 15:16 19:22 25:28	29	30	33	34	35	36	39	40]};

% define legend
legs = groupIDs;

% define plot colors
colors = {[0 0 1],[1 0 0]};

% define lap type labels
lapTypeName = {'binAll','bin5','bin1'};
groupSz = [1,5,19];
nCurLapGroup = length(lapTypeName);

% set save folder name
svFile = [p1 '\Figures\Figures' char(datetime('today','format','yyyyMMdd')) '\consistency\' gType];
mkdir(svFile);
mkdir([svFile '\data']);


%% Define cell types by genotype

% select cell types to analyze
load('data\globalCellTypes.mat', 'globalCellTypeLogical')
cellTypes = {'common','grid','cue'};
fieldIDs = fieldnames(globalCellTypeLogical);
cellTypeIDs = {7:12,10,7:9};

cellTypeIdx = cell(length(cellTypes),nGroups);

% caluclate group cell type logical
for ii = 1:length(cellTypes)
    curTypeID = cellTypeIDs{ii};
    for kk = 1:length(curTypeID)
        curTypeLogical = globalCellTypeLogical.(fieldIDs{curTypeID(kk)});
        for jj = 1:nGroups
            % full current logical
            curFull = curTypeLogical(groups{jj});

            idxRemove = find(~ismember(groups{jj},groupsFinal{jj}));
            for mm = 1:length(idxRemove)
                curFull{idxRemove(mm)} = false(size(curFull{idxRemove(mm)}));
            end

            % current group cell type logical
            groupTypeLogical = cat(1,curFull{:});

            % add cell type to logical
            if kk==1
                cellTypeIdx{ii,jj} = groupTypeLogical;
            else
                cellTypeIdx{ii,jj} = cellTypeIdx{ii,jj} | groupTypeLogical;
            end
        end
    end
end


%% Plot correlations all

close all

% define cell types
useCellTypes = [1 2 3];
useCellTypes = [1];

% define use types
useLapTypes = {[1 3],[1 3],[1 3]};
% useData = 'corrMMat';
% useData = 'corrMFOV';
useData = 'corrMInd';

% close all figures
clAll = 0;

% useCellTypes = 1;
% useLapTypes = {[1]};

for gg = 1:length(useCellTypes)
    curUseTypeLaps = useLapTypes{gg};

    for ii = 1:length(curUseTypeLaps)
        curDataRaw = dataAllInter(curUseTypeLaps(ii)).(useData);
        if strcmp(useData,'corrMInd')
            curDataRawSub = cellfun(@(x,y) x(:,:,y),curDataRaw,cellTypeIdx(useCellTypes(gg),:),'UniformOutput',false);
        else
            if useCellTypes(gg)>1; continue; end
            curDataRawSub = cellfun(@(x,y,z) x(:,:,ismember(y,z)),curDataRaw,groups,groupsFinal,'UniformOutput',false);
        end
        curDataMean = cellfun(@(x) mean(x,3,'omitnan'),curDataRawSub,'UniformOutput',false);
        curDataMax = max(cellfun(@(x) max(x(~eye(size(x))),[],'all'),curDataMean));

        for jj = 1:3
            %%
            if jj<=2
                % exract current data
                curData = curDataMean{jj};

                % set data plotting parameters
                curCMap = 'parula';
                curCLim = [0 curDataMax];
                svName = groupIDs{jj};
            else
                % extract data difference paramters
                curData = curDataMean{1}-curDataMean{2};
                diffMax = max(abs(curData(~eye(size(curData)))),[],'all');

                % set difference plotting
                curCMap = customcolormap([0 0.25 0.5 0.75 1],[0 0 1; 0.5 0.5 1; 1 1 1; 1 0.5 0.5; 1 0 0]);
                curCLim = [-diffMax diffMax];
                svName = 'diff';
            end

            nCmps = size(curData,1);

            figure

            % set diagonal to nan
            dataDiag = curData;
            dataDiag(eye(size(dataDiag))==1) = nan;

            % plot correlations
            imAlpha = ones(size(dataDiag));
            imAlpha(isnan(dataDiag))=0;
            imagesc(dataDiag,'AlphaData',imAlpha,'YData',[size(dataDiag,1) 1]);
            set(gca,'color',0*[1 1 1]);
            clim(curCLim)
            colormap(curCMap)
            colorbar
            hold on

            % plot grid
            [X,Y]=meshgrid(0.5:groupSz(curUseTypeLaps(ii)):nCmps+0.5);
            plot(X,Y,'k');
            plot(Y,X,'k');

            % turn axis off
            h = get(gca);
            axis('square')
            h.XAxis.Visible='off';
            h.YAxis.Visible='off';
            h.XAxis.Label.Visible='on';
            h.YAxis.Label.Visible='on';

            % set figure labels
            RIGHT  = char(8594);
            xlabel(['Lap group ' RIGHT])
            ylabel(['Lap group ' RIGHT])
            set(gca,'FontSize',12)
            title(['Interday consistency (' cellTypes{useCellTypes(gg)} ' cells): ' svName])

            % save figure
            savefig([svFile '\interDayMat-' useData '-' lapTypeName{curUseTypeLaps(ii)}...
                '_' svName '-' cellTypes{useCellTypes(gg)}])

            % save data
            save([svFile '\data\interDayMat-' useData '-' lapTypeName{curUseTypeLaps(ii)}...
                '_' svName '-' cellTypes{useCellTypes(gg)} '.mat'],'curData')

            % close figure
            if clAll
                close
            end
        end
    end
end


%% Plot adjacent lap curve

useLapTypes = 3;
useData = 'corrMInd';
% useData = 'corrMFOV';

% close all figures
clAll = 0;

% cycle through lap groupings
for ii = 1:length(useLapTypes)
    %% Plot figure

    curDataRaw = dataAllInter(useLapTypes(ii)).(useData);
    if strcmp(useData,'corrMInd')
        curDataRawSub = cellfun(@(x,y) x(:,:,y),curDataRaw,cellTypeIdx(1,:),'UniformOutput',false);
    else
        curDataRawSub = cellfun(@(x,y,z) x(:,:,ismember(y,z)),curDataRaw,groups,groupsFinal,'UniformOutput',false);
    end
    data = cell(1,2);

    % cycle through genotypes
    for jj = 1:2
        % exract current data
        curData = curDataRawSub{jj};
        nCmps = size(curData,1);
        curDaySz = groupSz(useLapTypes(ii));

        adjDiag = zeros(size(curData,3),size(curData,1)-1);
        for kk = 1:size(curData,3)
            adjDiag(kk,:) = diag(curData(:,:,kk),1);
        end

        adjDiagUse = adjDiag(:,curDaySz+1:end);
        adjDiagUse(:,mod(1:size(adjDiagUse,2),curDaySz)==0) = [];
        % adjDiagSmooth = zeros(size(adjDiagUse));
        % for kk = 1:size(curData,3)
        %     adjDiagSmooth(kk,:) = smooth(adjDiagUse(kk,:),round(curDaySz/2));
        % end
        % adjDiagSmooth = adjDiagUse;

        data{jj} = adjDiagUse;
    end

    % perform statistics
    lenCur = size(data{1},2);
    [pACur,pMC] = anovaRM2W_full_BH(data{1},data{2},1);
    pAnova = pACur([1 3])'.*ones(lenCur,2);

    % initialize figure
    figure; hold on
    X = 1:lenCur;

    % plot line graph
    plotErrorSig(X,data{1},data{2},legs,pAnova,pMC,colors)

    % plot day boundaries
    [X,~]=meshgrid(curDaySz-0.5:curDaySz-1:nCmps);
    plot(X(1:2,1:end-1),repmat([0;1],1,10),'color',[0 .5 0.5 0.5],'HandleVisibility','off');

    % set figure labels
    legend('Location','eastoutside')
    xlabel('Pre lap group')
    ylabel('correlation')
    title(['Adjacent Laps: ' lapTypeName{ii}])
    set(gca,'FontSize',12)

    % save figure
    savefig([svFile '\interDayAdjLaps-' useData '-' lapTypeName{ii} '.fig'])

    % save data
    save([svFile '\data\interDayAdjLaps-' useData '-' lapTypeName{ii} '.mat'],'X','data')

    % close figure
    if clAll
        close
    end
end


%% Plot adjacent day correlation curve

% define cell types
% useCellTypes = [1 2 3];
useCellTypes = [1];

% define use types
useLapTypes = {1,1,1};
useData = 'corrMInd';
% useData = 'corrMFOV';

dayCats = {1,2:10};
nCats = length(dayCats);

% set cell type axes
yAX = [0 1;0 1.2;0 1.2];

% whether to plot individual lines
indON = 0;

% close all figures
clAll = 0;

% cycle through cell types
for gg = 1:length(useCellTypes)
    curUseTypeLaps = useLapTypes{gg};

    % cycle through lap groupings
    for ii = 1:length(curUseTypeLaps)
        %% Plot figure

        curDataRaw = dataAllInter(curUseTypeLaps(ii)).(useData);

        if strcmp(useData,'corrMInd')
            curDataRawSub = cellfun(@(x,y) x(:,:,y),curDataRaw,cellTypeIdx(useCellTypes(gg),:),'UniformOutput',false);
        else
            if useCellTypes(gg)>1; continue; end
            curDataRawSub = cellfun(@(x,y,z) x(:,:,ismember(y,z)),curDataRaw,groups,groupsFinal,'UniformOutput',false);
        end

        data = cell(1,2);

        % cycle through genotypes
        for jj = 1:2
            % exract current data
            curData = curDataRawSub{jj};

            adjDiag = zeros(size(curData,3),size(curData,1)-1);
            for kk = 1:size(curData,3)
                adjDiag(kk,:) = diag(curData(:,:,kk),1);

            end

            data{jj} = adjDiag;
        end

        nCmps = size(data{1},2);

        Xname = cell(1,nCmps);
        X = 1:nCmps;
        for kk = 1:nCmps
            Xname{kk} = [num2str(kk-1) '-' num2str(kk)];
        end

        % initialize statistics
        pMC = nan(nCmps,1);
        pAnova = nan(nCmps,2);

        % perform statistics
        for jj = 1:nCats
            % find intersect of plotting x values with day category
            curCat = find(ismember(X,dayCats{jj}));
            lenCur = length(curCat);
            if lenCur>1
                % calculate anova p values with multiple comparisons
                [pACur,pMCCur] = anovaRM2W_full_BH...
                    (data{1}(:,curCat),data{2}(:,curCat),1);
                pAnova(curCat,:) = pACur([1 3])'.*ones(lenCur,2);
            else
                % calculat single day p values
                [~,pMCCur] = ttest2(data{1}(:,curCat),data{2}(:,curCat));
            end

            pMC(curCat) = pMCCur;
        end

        % check individual line condition
        if indON==1 && ~strcmp(useData,'corrMFOV')
            disp('Individual values can only be ploted by FOV')
            return
        end

        % initialize figure
        figure; hold on

        % plot line graph
        h = plotErrorSig(X,data{1},data{2},legs,pAnova,pMC,colors,indON);

        delete(h)
        h = zeros(1,2);
        for jj = 1:2
            dataMean = mean(data{jj},1,'omitnan');
            dataSEM = nansem(data{jj},1);
            h(jj) = errorbar(X(1),dataMean(1),dataSEM(1),'color',colors{jj},'LineWidth',1);
            errorbar(X(2:end),dataMean(2:end),dataSEM(2:end),'color',colors{jj},'LineWidth',1);
        end

        % set figure labels
        legend(h,legs,'Location','southeast')
        xlabel('Correlation day pair')
        ylabel('Mean activity correlation')
        title(['Adjacent Day Correlation (' cellTypes{useCellTypes(gg)} ' cells)'])
        set(gca,'FontSize',12)

        % set x axis and labels
        set(gca,'XTick',X)
        set(gca,'XTickLabels',Xname)
        xlim([X(1)-0.5 X(end)+0.5])

        % set and break y axis
        ylim(yAX(gg,:))
        % breakyaxis([0.2 0.45]);

        % save figure
        savefig([svFile '\interDayAdjDays-' useData '-' lapTypeName{ii} '-' cellTypes{useCellTypes(gg)} '.fig'])

        % save data
        save([svFile '\data\interDayAdjDays-' useData '-' lapTypeName{ii} '-' cellTypes{useCellTypes(gg)} '.mat'],'X','data')

        % close figure
        if clAll
            close
        end
    end
end


%% Plot inter-day memory correlation curve

useLapTypes = 1;
useData = 'corrMInd';
% useData = 'corrMFOV';

useDays = 2:11;     % which days to use for comparison
meanDays = 9:11;    % which days to use as reference
intersectDays = ismember(useDays,meanDays);

% close all figures
clAll = 0;

% cycle through lap groupings
for ii = 1:length(useLapTypes)
    %% Plot figure

    curDataRaw = dataAllInter(useLapTypes(ii)).(useData);
    if strcmp(useData,'corrMInd')
        curDataRawSub = cellfun(@(x,y) x(:,:,y),curDataRaw,cellTypeIdx(1,:),'UniformOutput',false);
    else
        curDataRawSub = cellfun(@(x,y,z) x(:,:,ismember(y,z)),curDataRaw,groups,groupsFinal,'UniformOutput',false);
    end

    data = cell(1,2);

    % cycle through genotypes
    for jj = 1:2
        % exract current data
        curData = curDataRawSub{jj}(useDays,useDays,:);
        curData(:,~intersectDays,:) = nan;

        % calculate mean of all diagonals
        adjDiag = zeros(size(curData,3),size(curData,1)-1);
        for kk = 1:size(curData,3)
            for mm = 1:size(curData,1)-1
                curDiag = diag(curData(:,:,kk),mm);
                adjDiag(kk,mm) = mean(curDiag,'omitnan');
            end
        end

        data{jj} = adjDiag;
    end

    % define x axis
    nCmps = size(data{1},2);
    X = 1:nCmps;
    Xname = string(X);

    % calculate anova p values with multiple comparisons
    [pACur,pMC] = anovaRM2W_full_BH(data{1},data{2},1);
    pAnova = pACur([1 3])'.*ones(nCmps,2);

    % initialize figure
    figure; hold on

    % plot line graph
    plotErrorSig(X,data{1},data{2},legs,pAnova,pMC,colors);

    % set figure labels
    legend(legs,'Location','southwest')
    xlabel('Day pair distance (days)')
    ylabel('Mean activity correlation')
    title('Day Correlation Memory')
    set(gca,'FontSize',12)

    % set x axis and labels
    set(gca,'XTick',X)
    set(gca,'XTickLabels',Xname)
    xlim([X(1)-0.5 X(end)+0.5])

    % set y axis
    ylim([0.3 0.9])

    % save figure
    savefig([svFile '\interDayMemDays-' useData '-' lapTypeName{ii} '.fig'])

    % save data
    save([svFile '\data\interDayMemDays-' useData '-' lapTypeName{ii} '.mat'],'X','data')

    % close figure
    if clAll
        close
    end

    % prep for Prism
    dataMean = cellfun(@(x) nanmean(x,1)',data,'UniformOutput',0);
    dataMean = cat(2,dataMean{:});

end


%% Plot intra-day memory correlation surface

useLapTypes = 3;
useData = 'corrMInd';
% % useData = 'corrMFOV';

useDays = 2:11;
nDays = length(useDays);

% close all figures
clAll = 0;

xRange = 1:10;
yRange = 1:15;
nX = length(xRange);
nY = length(yRange);


% cycle through lap groupings
for ii = 1:length(useLapTypes)
    %% Calculate memory surfaces

    curDataRaw = dataAllInter(useLapTypes(ii)).(useData);
    if strcmp(useData,'corrMInd')
        curDataRawSub = cellfun(@(x,y) x(:,:,y),curDataRaw,cellTypeIdx(1,:),'UniformOutput',false);
    else
        curDataRawSub = cellfun(@(x,y,z) x(:,:,ismember(y,z)),curDataRaw,groups,groupsFinal,'UniformOutput',false);
    end
    
    data = cell(1,2);

    % define use laps
    curGroupSz = groupSz(useLapTypes(ii));
    startLaps = (useDays-1)*curGroupSz+1;
    stopLaps = useDays*curGroupSz;
    useLaps = [startLaps', stopLaps'];

    % cycle through genotypes
    for jj = 1:2
        data{jj} = zeros(curGroupSz-1,nDays,size(curDataRawSub{jj},3));

        % cycle through days
        for gg = 1:nDays
            % exract current data
            useLapsExp = useLaps(gg,1):useLaps(gg,2);
            curData = curDataRawSub{jj}(useLapsExp,useLapsExp,:);

            % calculate mean of all diagonals
            adjDiag = zeros(size(curData,3),size(curData,1)-1);
            for kk = 1:size(curData,3)
                for mm = 1:size(curData,1)-1
                    curDiag = diag(curData(:,:,kk),mm);
                    adjDiag(kk,mm) = mean(curDiag,'omitnan');
                end
            end

            data{jj}(:,gg,:) = adjDiag';
        end
    end


    %% Plot figure

    figure; hold on
    dataMean = cellfun(@(x) mean(x,3,'omitnan'),data,'UniformOutput',false);
    curDataMin = min(cellfun(@(x) min(x(~eye(size(x))),[],'all'),dataMean));
    curDataMax = max(cellfun(@(x) max(x(~eye(size(x))),[],'all'),dataMean));

    for jj = 1:2
        subplot(1,2,jj);

        % plot surface
        curDataMean = dataMean{jj};
        h = imagesc(flipud(curDataMean(yRange,xRange)));

        % define view
        col = colorbar;
        clim([curDataMin curDataMax])
        curCMap = customcolormap([0 0.5 1],[0 0 0; 0.5 0.5 0.5; 1 1 1]);
        colormap(curCMap)
        colormap('jet')

        % set plot labels
        title(groupIDs{jj})
        xlabel('Day')
        ylabel('Lap distance')
        ylabel(col,'RBR correlation')
        set(gca,'FontSize',18)
        set(col,'FontSize',18)

        % set axes limits
        axis('equal')
        xlim([xRange(1)-0.5 xRange(end)+0.5])
        ylim([yRange(1)-0.5 yRange(end)+0.5])
        set(gca,'YTick',yRange)
        set(gca,'yTickLabel',flipud(get(gca,'YTickLabel')))
    end

    % calculate p values
    pData = cellfun(@(x) x(yRange,xRange,:),data,'UniformOutput',false);
    M1 = mean(pData{1},'all','omitnan');
    M2 = mean(pData{2},'all','omitnan');
    [pAnova,pLabels] = anovaRM3W(pData,{'LapDistance','Day'});

    ttl = [];
    for jj = 1:length(pAnova)
        ttl = [ttl ', ' pLabels{jj} ': ' num2str(pAnova(jj),2)];
    end
    sgtitle(['WT mean: ' num2str(M1,3) ', AD mean: ' num2str(M2,3) ttl])

    % save figure
    savefig([svFile '\intraDaySurface-' useData '-' lapTypeName{ii} '.fig'])

    % save data
    save([svFile '\data\intraDaySurface-' useData '-' lapTypeName{ii} '.mat'],'xRange','yRange','dataMean')

    % close figure
    if clAll
        close
    end


    %% Plot p vlaue grid

    % calculate p-values
    pAll= zeros(nY,nX);
    signCur = zeros(nY,nX);
    for jj = 1:nY
        for kk = 1:nX
            [~,pAll(jj,kk)] = ttest2(pData{1}(jj,kk,:),pData{2}(jj,kk,:));
            signCur(jj,kk) = mean(pData{1}(jj,kk,:),'omitnan') > mean(pData{2}(jj,kk,:),'omitnan');
        end
    end

    % [~,pAllMC] = bonferroni_holm(pAll);
    % pAllMC2 = holm_sidak(pAll);
    pAllMC = pAll;

    pAllPlot = -log10(pAllMC).*(signCur*2-1);

    figure; hold on
    h = imagesc(pAllPlot);

    % define view
    col = colorbar;
    cRange = [-4,4];
    cBar = [0.051 0.05];
    cBarMap = abs(log10(cBar))/range(cRange);

    clim(cRange)
    curCMap = customcolormap([0 0.5-cBarMap(1) 0.5-cBarMap(2) 0.5+cBarMap(2) 0.5+cBarMap(1) 1],[0 0 1; 0.95 0.95 1; 1 1 1; 1,1,1; 1 0.95 0.95; 1 0 0]);
    colormap(curCMap)

    % set plot labels
    title('p-value table')
    xlabel('Day')
    ylabel('Lap distance')
    ylabel(col,'p-value (log10)')
    set(gca,'FontSize',18)
    set(col,'FontSize',18)

    % set axes limits
    axis('equal')
    xlim([xRange(1)-0.5 xRange(end)+0.5])
    ylim([yRange(1)-0.5 yRange(end)+0.5])
    yticks(yRange)
    yticklabels(yRange)
    xticks(xRange)
    xticklabels(xRange)


    % save figure
    savefig([svFile '\intraDaySurface-p-' useData '-' lapTypeName{ii} '.fig'])

    % save data
    save([svFile '\data\intraDaySurface-p-' useData '-' lapTypeName{ii} '.mat'],'xRange','yRange','pAllPlot')

    % close figure
    if clAll
        close
    end


end


%% Plot inter-day activity matrices

useLapTypes = 1;

% set comparison days
plotDays = [2 3;10 11];
nPlots = size(plotDays,1);
nPlotCmp = size(plotDays,2);

% set titles
ttls = {'Day 1','Day 2';'Day 9','Day 10'};

allMatrix = cell(nGroups,1);

for ii = 1:nGroups

    allMatrix{ii} = cell(nPlots,nPlotCmp);


    figure
    sgtitle(groupIDs{ii})

    for jj = 1:nPlots

        for kk = 1:nPlotCmp
            %% Normalize and sort data

            curMatrix = squeeze(dataAllInter(useLapTypes).meanCombCell{ii,plotDays(jj,kk)})';

            % get maximum row values and indices
            [dataMax,idxMax] = max(curMatrix,[],2);

            % normalize by row maximum
            curMatrix = curMatrix./dataMax;

            if kk==1
                % find sorted row order
                [~,idxSort] = sort(idxMax);
            end

            % sort by reference rows
            curMatrix = curMatrix(idxSort,:);


            %% Plot sorted matrix

            % set subplot
            subplot(nPlots,nPlotCmp,(jj-1)*nPlots+kk)

            % plot matrix
            imagesc(curMatrix,[0 1])
            axis('square')


            title(ttls{jj,kk})

            allMatrix{ii}{jj,kk} = curMatrix;
        end
    end

    % save figure
    savefig([svFile '\interDayActivityMat-' lapTypeName{ii} '_' groupIDs{ii} '.fig'])

    % close figure
    if clAll
        close
    end

end


% save data
save([svFile '\data\interDayActivityMat-' lapTypeName{ii} '.mat'],'legs','allMatrix')


%% Plot inter-day consistency bar graph

% define legend
legs4 = {'WT_{Pre-Learning}','AD_{Pre-Learning}';'WT_{Post-Learning}','AD_{Post-Learning}'};

% define plot colors
colors4 = {[0.5 0.5 1],[1 0.5 0.5];[0 0 1],[1 0 0]};

dayGroups = {2,8:10};
useLapTypes = 1;
useData = 'corrMInd';
% useData = 'corrMFOV';

% close all figures
clAll = 0;

curDataRaw = dataAllInter(useLapTypes).(useData);

if strcmp(useData,'corrMInd')
    curDataRawSub = cellfun(@(x,y) x(:,:,y),curDataRaw,cellTypeIdx(1,:),'UniformOutput',false);
else
    curDataRawSub = cellfun(@(x,y,z) x(:,:,ismember(y,z)),curDataRaw,groups,groupsFinal,'UniformOutput',false);
end

data = cell(1,2);

% cycle through genotypes
for jj = 1:2
    % exract current data
    curData = curDataRawSub{jj};

    adjDiag = zeros(size(curData,3),size(curData,1)-1);
    for kk = 1:size(curData,3)
        adjDiag(kk,:) = diag(curData(:,:,kk),1);
    end

    data{jj} = adjDiag;
end

barData = {};
% generate bar data
for jj = 1:2
    for kk = 1:length(dayGroups)
        barData{end+1} = mean(data{jj}(:,dayGroups{kk}),2);
    end
end

% define t-test pairs
pShow = [1 2;1 3;2 4;3 4];

% initialize figure
figure; hold on

% plot bar graph
[hh,pCorr] = barGroup(legs4,barData,colors4,pShow);

% set figure labels
ylim([0 1])
ylabel('Mean activity correlation')
title('Inter-day consistency correlation')
set(gca,'FontSize',12)

% save figure
savefig([svFile '\interDayBar-' useData '-' lapTypeName{ii} '.fig'])

% save data
save([svFile '\data\interDayBar-' useData '-' lapTypeName{ii} '.mat'],'legs4','barData')


%% Interpolate function

function dfofMInterpM_sig = findDfofInterM(dfofM,lapIdx)

% initialize output array
dfofMInterpM_sig = cell(size(dfofM));

for n = 1:length(dfofM)
    % initialize cell matrix
    dfofMInterpM_sig{n} = NaN(length(lapIdx),size(dfofM{1},2));

    % load only use laps
    thisdfofM = dfofM{n};

    if ~isempty(thisdfofM)
        % interpolate by row
        for m = 1:length(lapIdx)
            dfofMInterpM_sig{n}(m,:) = naninterp(thisdfofM(lapIdx(m),:));
        end
    end
end

end

