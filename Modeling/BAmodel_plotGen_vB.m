%% BAmodel_plotGen
% Makes summary plots using a final model determined elsewhere. This
% generic version can plot based on training and test data, using lme of
% glm, and using LOOCV or LTOCV. BAmodel_postProcess should be run first to
% process data.
%


%% Initialize plotting

clear; close all;

% set folder
p1 = 'D:\AD_Project\Behavior';
cd(p1)

% select data
modelType = 'glm';
useGroup = 1;
groupSaveName = {'All-LOO','WT-LOO','PS19-LOO','All-LTO'};
svSuff = '_SG';

% whether to show equation
showEQ = 1;

% load processed data
load(['data/model_vB/LOOCVplot_' modelType '_' groupSaveName{useGroup} svSuff '_vB.mat'],...
    'dataProcessed','LOOData','nTrials','tblNaN')

% load groups
load('groupIDs.mat','groups','groupIDs','sexes','sexIDs')
nGroups = length(groups);
nSexes = length(sexes);
groupIDsExt = [groupIDs, 'All'];

% define parameter names
paramNames = dataProcessed.labsFinal';
% paramNames = {'Corr_{interLap}';'Corr_{SF}';'\DeltaF/F_{Sig}';'RdgBkg'};
nFinal = length(paramNames);

% set save folder name
svFile = [p1 '/Figures/Figures' char(datetime('today','format','yyyyMMdd')) '/BAmodel_vB'];
if ~isfolder(svFile)
    mkdir(svFile);
    mkdir([svFile '/data'])
end
svLabel = ['_' modelType '_' groupSaveName{useGroup} svSuff];

% define plot colors
colors2 = {[0 0 1],[1 0 0]};
colors3 = {[0 0 1],[1 0 0],[0 0 0]};
colors4 = {[102 105 255]/255,[255 102 102]/255;[0 0 175]/255,[175 0 0]/255};
colors6 = {[0 0 1],[1 0 0];[102 105 255]/255,[255 102 102]/255;[0 0 175]/255,[175 0 0]/255};

% define labels and data
szText = 16;
dataTypes = {'domTrain','domTest'};
dataLabs = {'train','test'};
SEVar = 'R2';


%% Plot data

% initialize stats parameters and array
outStats = {};
testNameCorr = 'two-tailed Pearson linear correlation';
nUnits = 'days';

for dt = 1:length(dataTypes)
    %% Plot contribution waterfall

    % get contribution data
    waterData = dataProcessed.(dataTypes{dt}).R2.domTot'*100;

    % sort data
    if dt==1
        [waterDataSort,sortIdx] = sort(waterData,'descend');
        labelsSort = paramNames(sortIdx);

        % get coefficients
        coeff = dataProcessed.coeff([1;sortIdx+1]);

        % generate equation
        EQ = ['y_{beh} = ' num2str(coeff(1),2)];
        for xx = 1:nFinal
            EQ = [EQ ' + ' num2str(coeff(xx+1),2) '*X_{' labelsSort{xx} '}'];
        end
    else
        waterDataSort = waterData(sortIdx);
    end

    % add unexplained variance
    waterDataFull = [waterDataSort; 100-sum(waterData)];
    labelsFull = [labelsSort;'unexplained'];

    % generate stack input data
    waterPlot = [[0;cumsum(waterDataFull(1:end-1))],waterDataFull];

    % plot waterfall bar graph
    figure; hold on
    h = bar(labelsFull,waterPlot,'stacked');

    % remove underbar
    delBar = 1;
    h(delBar).FaceColor = 'none';
    h(delBar).EdgeColor = 'none';

    % add text labels
    nWater = size(waterPlot,1);
    sumVals = sum(waterPlot,2);
    for ii = 1:nWater
        text(ii,sumVals(ii),num2str(waterDataFull(ii),'%.1f'),'FontSize',szText,...
            'VerticalAlignment','bottom','HorizontalAlignment','center')
    end

    if showEQ==1
        % add equation as text
        text(0.5,100,EQ,'VerticalAlignment','middle','HorizontalAlignment','left',...
            'FontSize',szText-2)
    end

    % set plot labels
    xlabel('Predictors')
    ylabel('% variance explained (R^2)')
    title (['Relative Importance of Predictors (' dataLabs{dt} ')'])
    ylim([0 110])
    set(gca,'FontSize',szText)

    % save figure
    savefig([svFile '/waterfall_' dataLabs{dt} svLabel '.fig'])

    % save data
    save([svFile '/data/waterfall_' dataLabs{dt} svLabel '.mat'],'labelsFull','waterPlot')


    %% Plot contribution pie chart

    if dt==1
        % define pie order and colors
        CP = orderedcolors('gem');
        pieOrder = [1 2 3 6 7 5 4 8 9];
        pieColors = [CP([1 1 4 4 4 2 5 5],:); [0.5 0.5 0.5]];

        % apply pie order
        pieData = waterDataFull(pieOrder);
        pieLabels = labelsFull(pieOrder);

        % plot pie chart
        figure
        piechart(pieData,pieLabels)
        colororder(pieColors)

        % set plot labels
        title (['Relative Importance of Predictors (' dataLabs{dt} ')'])
        set(gca,'FontSize',szText)

        % save figure
        savefig([svFile '/pie_' dataLabs{dt} svLabel '.fig'])

        % save data
        save([svFile '/data/pie_' dataLabs{dt} svLabel '.mat'],'labelsFull','waterPlot')
    end

    %% Plot sub/expansion model factor contribution

    % define sub labels, p, and r
    pSub = dataProcessed.(dataLabs{dt}).(SEVar).pRem';
    rSub = dataProcessed.(dataLabs{dt}).(SEVar).Rem';
    nSub = length(paramNames);

    % sort sub factors
    [rSubSort,sortSub] = sort(rSub,'ascend');
    pSubSort = pSub(sortSub);
    labsSubSort = paramNames(sortSub);

    % define final labels, p, and r
    labsFinal = 'FinalModel';
    rFinal = dataProcessed.(dataLabs{dt}).(SEVar).Pred;
    pFinal = NaN;

    % define expansion labels, p, and r
    labsExp = string(dataProcessed.labsOther)';
    pExp = dataProcessed.(dataLabs{dt}).(SEVar).pExp';
    rExp = dataProcessed.(dataLabs{dt}).(SEVar).Exp';
    nExp = length(labsExp);

    % sort expansion factors
    [rExpSort,sortExp] = sort(rExp,'descend');
    pExpSort = pExp(sortExp);
    labsExpSort = labsExp(sortExp);

    % define combined variables
    plotLab = [labsSubSort; labsFinal; labsExpSort];
    plotP = [pSubSort; pFinal; pExpSort];
    plotR = [rSubSort; rFinal; rExpSort];
    plotIdx = 1:length(plotLab);

    % define colors and sizes
    colorsScatt = [0.2, 0.5, 0.9;0.85, 0.1, 0.1;0.2, 0.8, 0.3];
    szScatt = [80 100 60];

    % set colors and sizes
    plotSizes = [repmat(szScatt(1),nSub,1); szScatt(2); repmat(szScatt(3),nExp,1)];
    plotColors = [repmat(colorsScatt(1,:),nSub,1); colorsScatt(2,:); repmat(colorsScatt(3,:),nExp,1)];

    % plot scatter
    figure; hold on
    scatter(plotIdx,plotR,plotSizes,plotColors,'filled')

    % set labels
    title(['Factor Interaction Contributions (' dataLabs{dt} ')'] )
    xlim([0.5 plotIdx(end)+0.5])
    xticks(plotIdx)
    xticklabels(plotLab)
    xtickangle(90)
    ylabel('% variance explained (R^2)')
    set(gca,'FontSize',12)

    % disp(num2str(plotP,2))

    % save figure
    savefig([svFile '/SubExp_' dataLabs{dt} svLabel '.fig'])

    % save data
    save([svFile '/data/SubExp_' dataLabs{dt} svLabel '.mat'],'plotLab','plotR')


    %% Plot single-mouse examples

    % load mouse data
    idxMouse = unique(LOOData.idxMouse(:));
    tblM = LOOData.tblM;
    nMouse = length(idxMouse);
    nRow = ceil(nMouse^0.5);

    % load prediction data
    curSuff = [upper(dataLabs{dt}(1)) dataLabs{dt}(2:end)];
    dataReal = LOOData.(['y' curSuff]);
    dataPred = LOOData.(['yPred' curSuff]);
    nItr = size(dataReal,1);

    % concatenate data by mouse
    mouseReal = cell(1,nMouse);
    mousePred = cell(1,nMouse);
    for mm = 1:nMouse
        % collect all data and predictions for current mouse
        curMouseReal = cell(1,nItr);
        curMousePred = cell(1,nItr);

        for ii = 1:nItr
            % get current mouse labels
            idxUse = xor(LOOData.idxUse{ii},dt==2);
            curIdxM = tblM(idxUse)==idxMouse(mm);

            % store current mouse data
            curMouseReal{ii} = dataReal{ii}(curIdxM);
            curMousePred{ii} = dataPred{ii}(curIdxM);
        end

        % remove empty cells
        idxEmpty = cellfun(@isempty,curMouseReal);
        curMouseReal = curMouseReal(~idxEmpty);
        curMousePred = curMousePred(~idxEmpty);

        % concatenate and take mean
        mouseReal{mm} = mean(cat(2,curMouseReal{:}),2,'omitnan')/nTrials;
        mousePred{mm} = mean(cat(2,curMousePred{:}),2,'omitnan')/nTrials;
    end

    figure;
    tiledlayout(3,7)
    sgtitle(['Individual Mouse Predictions (' dataLabs{dt} ')'])
    combinedData = cell(nGroups,nSexes,2);
    for gg = 1:nGroups
        for ss = 1:nSexes
            for ii = 1:2
                combinedData{gg,ss,ii} = {};
            end
        end
    end


    idxCount = 0;
    for ss = nSexes:-1:2
        for gg = 1:nGroups
            idxCurrent = 0;

            for mm = 1:nMouse
                % skip mice not in current group
                if ~ismember(mm,intersect(groups{gg},sexes{ss}))
                    continue
                end
                idxCount = idxCount+1;
                idxCurrent = idxCurrent+1;

                % get current mouse nans
                keepIdxCur = tblNaN.keepIdx(ismember(tblNaN.mouse,num2str(mm)));

                % get current mouse data
                curReal = nan(size(keepIdxCur));
                curPred = nan(size(keepIdxCur));
                curReal(keepIdxCur) = mouseReal{mm};
                curPred(keepIdxCur) = mousePred{mm};

                % get current mouse label
                curGroup = find(cellfun(@(x) ismember(idxMouse(mm),x),groups));
                curSex = find(cellfun(@(x) ismember(idxMouse(mm),x),sexes(2:end)));

                % calculate correlation
                [r,p] = corr(curReal,curPred,'rows','pairwise');
                R2 = findR2(curReal,curPred);

                % plot current mouse
                nexttile(idxCount); hold on
                plot(curReal,'-','Color',colors4{curSex,curGroup},'LineWidth',1)
                plot(curPred,'-','Color',[0.5 0.5 0.5],'LineWidth',1.5)
                xlabel('Session')
                ylabel('Behavior')
                xlim([0.5 10.5])
                ylim([0 1])

                legend(['p = ' num2str(p,2)],'Location','best')
                title([groupIDs{curGroup} '-' sexIDs{curSex+1} ' Mouse ' num2str(idxCurrent)...
                    ': r = ' num2str(r,2) ', R^2 = ' num2str(R2,2)],'FontSize',5)

                % store data
                combinedData{curGroup,1,1}{end+1} = curReal;
                combinedData{curGroup,curSex+1,1}{end+1} = curReal;
                combinedData{curGroup,1,2}{end+1} = curPred;
                combinedData{curGroup,curSex+1,2}{end+1} = curPred;

                % calculate full statistics
                corrCat = ['Group Predictions (' dataLabs{dt} '): '...
                    groupIDs{gg} ' ' sexIDs{ss} ' ' num2str(idxCurrent)];
                outStats(end+1,:) = [corrCat corrEffectSize(curReal,curPred,testNameCorr,nUnits)];
            end
        end
    end

    % save figure
    savefig([svFile '/indMouseCurves_' dataLabs{dt} svLabel '.fig'])

    % save data
    save([svFile '/data/indMouseCurves_' dataLabs{dt} svLabel '.mat'],'mouseReal','mousePred')


    %% Plot group averaged time courses

    figure
    t = tiledlayout(nGroups,nSexes);
    sgtitle(['Group Predictions (' dataLabs{dt} ')'])

    for gg = 1:nGroups
        for ss = 1:nSexes
            % get current mouse data
            realCat = cat(2,combinedData{gg,ss,1}{:});
            predCat = cat(2,combinedData{gg,ss,2}{:});
            realMean = mean(realCat,2,'omitnan');
            predMean = mean(predCat,2,'omitnan');

            % skip if empty
            if isempty(realMean); continue; end

            % calculate correlation
            [r,p] = corr(realMean,predMean);
            R2 = findR2(realMean,predMean);

            % plot current mouse
            nexttile(tilenum(t,gg,ss)); hold on
            plot(realMean,'-','Color',colors6{ss,gg},'LineWidth',1)
            plot(predMean,'-','Color',[0.5 0.5 0.5],'LineWidth',1.5)
            xlabel('Session')
            ylabel('Behavior')
            xlim([0.5 10.5])
            ylim([0 1])

            title([groupIDs{gg} ' ' sexIDs{ss} ': r = ' num2str(r,2)...
                ', p = ' num2str(p,2) ', R^2 = ' num2str(R2,2)])

            % calculate full statistics
            corrCat = ['Group Predictions (' dataLabs{dt} '): ' groupIDs{gg} ' ' sexIDs{ss}];
            outStats(end+1,:) = [corrCat corrEffectSize(realMean,predMean,testNameCorr,nUnits)];
        end
    end

    % save figure
    savefig([svFile '/groupMouseCurves_' dataLabs{dt} svLabel '.fig'])

    % save data
    save([svFile '/data/groupMouseCurves_' dataLabs{dt} svLabel '.mat'],'combinedData')


    %% Plot full prediction validation

    figure; hold on

    h = [];
    labelLeg = {};
    for gg = 1:nGroups+1
        % collect trendline data
        if gg<=nGroups
            plotReal = cat(1,combinedData{gg,1,1}{:});
            plotPred = cat(1,combinedData{gg,1,2}{:});
        else
            plotReal = [cat(1,combinedData{1,1,1}{:}); cat(1,combinedData{2,1,1}{:})];
            plotPred = [cat(1,combinedData{1,1,2}{:}); cat(1,combinedData{2,1,2}{:})];
        end
        plotReal = plotReal(~isnan(plotReal));
        plotPred = plotPred(~isnan(plotPred));

        % skip if empty
        if isempty(plotReal); continue; end

        % calculate trendline
        coefficients = polyfit(plotReal,plotPred,1);
        xFit = [min(xlim) max(xlim)];
        yFit = polyval(coefficients, xFit);

        % plot trendline
        h(end+1) = plot(xFit,yFit,'-','Color',colors3{gg},'LineWidth',0.5);

        % calculate correlation
        [rCorr,pCorr] = corr(plotReal,plotPred);
        R2Corr = findR2(plotReal,plotPred);

        % define legend
        labelLeg{end+1} = [groupIDsExt{gg} ': r=' num2str(rCorr,2)...
            ', p=' num2str(pCorr,2) ', R^2=' num2str(R2Corr,2)];

        if gg<=nGroups
            % plot scatter by sex
            for ss = 1:2
                plotReal = cat(1,combinedData{gg,ss+1,1}{:});
                plotReal = plotReal(~isnan(plotReal));
                plotPred = cat(1,combinedData{gg,ss+1,2}{:});
                plotPred = plotPred(~isnan(plotPred));

                % plot points
                scatter(plotReal,plotPred,25,colors4{ss,gg},'filled')
            end
        end
    end

    % set plot labels
    xlabel('Behavior')
    ylabel('Model Prediction')
    title(['Behavior versus Prediction (' dataLabs{dt} ')'])
    legend(h,labelLeg,'Location','best','AutoUpdate','off')
    set(gca,'FontSize',18)
    axis('square')
    plot([0 1],[0 1],'-')

    % save figure
    savefig([svFile '/groupPredScatter_' dataLabs{dt} svLabel '.fig'])

    % save data
    save([svFile '/data/groupPredScatter_' dataLabs{dt} svLabel '.mat'],'combinedData')


    % %% Plot data distribution
    %
    % figure; hold on
    %
    % ksdensity(plotReal,'Function','pdf')
    % ksdensity(plotPred,'Function','pdf')


end


% %%
% 
% X = {'Anchoring_{decode}','Anchoring_{\DeltaF/F}','Genotype','Corr_{inter-lap}',...
%     'Mean \DeltaF/F','Spatial selectivity','Sex','Unexplained'};
% 
% xticklabels(X)
% xtickangle(0)
