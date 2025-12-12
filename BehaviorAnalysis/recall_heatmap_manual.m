%%

clear; close all; clc

types = {'raw','relativeGroup'};

for tt = 1:2
    % vals = cell(2,2,2); % ste/pyr, M/F, WT/PS19

    valsAll = cell(2,3);
    for ff = 1:3
        if ff==1
            load(['D:\AD_Project\Behavior\Figures\manuscriptFigures\recall\data\recall_perSuccess_' types{tt} '_all.mat'],'data')
        elseif ff==2
            load(['D:\AD_Project\Behavior\Figures\manuscriptFigures\recall\data\recall_perSuccess_' types{tt} '_female.mat'],'data')
        elseif ff==3
            load(['D:\AD_Project\Behavior\Figures\manuscriptFigures\recall\data\recall_perSuccess_' types{tt} '_male.mat'],'data')
        end

        valsAll(:,ff) = data([2 4]);
    end


    %%

    meansAll = cellfun(@(x) mean(x,'omitnan'),valsAll);
    matPlot = -diff(meansAll,1);

    pAll = zeros(1,3);
    for ff = 1:3
        [~,pAll(ff)] = ttest2(valsAll{1,ff},valsAll{2,ff});
    end
    disp(pAll)

    figure; hold on
    title(types{tt})

    F = imagesc(matPlot,'YData',[size(matPlot,1) 1]);

    col = colorbar();

    if tt==1
        curCLim = [-40,40];
    else
        curCLim = [-0.5 0.5];
    end

    curCMap = customcolormap([0 0.25 0.5 0.75 1],[0.75 0 0.75; 1 0.5 1; 1 1 1; 0.5 1 0.5; 0 0.75 0]);

    ylabel(col,'AD deficit')

    % set colors
    clim(curCLim)
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
    set(gca,'XTickLabel',{'All','Female','Male'});
    xlabel('Sex')

    % set y axis labels
    set(gca,'YTick',[]);
end