%%

clear; close all; clc

folderBase = 'D:\AD_Project\Behavior\Figures\manuscriptFigures\recall\data\';
sexes = {'all','female','male'};
varNames = {'perSuccess_raw','dPrimeGlobal_raw'};
nSexes = length(sexes);
nVars = length(varNames);
useDays = 2:11;

cLims = {[-40,40],[-1.4,1.4]};



%% Loop through variables

for vv = 1:nVars
    %% Collect data

    valsAll = cell(2,3);
    for ss = 1:3
        if ss==1
            sexFold = 'all';
        else
            sexFold = 'bySex';
        end
        load([folderBase 'recall_' varNames{vv} '_' sexes{ss} '.mat'],'data')
        valsAll(:,ss) = data([2 4]);
    end

    valsPlot = cellfun(@(x) mean(x,2,'omitnan'),valsAll,'UniformOutput',false);

    % calculate behavior difference
    meansAll = cellfun(@(x) mean(x,'omitnan'),valsPlot);
    matPlot = -diff(meansAll,1);

    % calculate significance
    pAll = zeros(1,3);
    for ss = 1:3
        [~,pAll(ss)] = ttest2(valsPlot{1,ss},valsPlot{2,ss});
    end

    %% Plot heatmap

    figure; hold on
    title(varNames{vv})

    F = imagesc(matPlot([1 3 2]),'YData',[size(matPlot,1) 1]);

    col = colorbar();

    curCMap = customcolormap([0 0.25 0.5 0.75 1],[0.75 0 0.75; 1 0.5 1; 1 1 1; 0.5 1 0.5; 0 0.75 0]);

    ylabel(col,'AD deficit')

    % set colors
    clim(cLims{vv})
    colormap(curCMap);

    % set plot sizes
    nX = 3;
    nY = 1;
    curXRange = [0.5 nX+0.5];
    curYRange = [0.5 nY+0.5];
    labelX = 'Sex';

    % set plot limits
    axis('equal')
    xlim(curXRange)
    ylim(curYRange)
    set(gca,'FontSize',16)
    set(col,'FontSize',16)

    % set x axis labels
    set(gca,'XAxisLocation','top');
    set(gca,'XTick',1:nX);
    set(gca,'XTickLabel',{'All','Male','Female'});
    xlabel('Sex')

    disp(sexes)
    disp(pAll)

    % set y axis labels
    set(gca,'YTick',[]);
end

