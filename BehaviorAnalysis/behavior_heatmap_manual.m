%%

clear; close all; clc

load('D:\AD_Project\Behavior\Figures\Figures20241218_fullUpdate\all\data\perSuccess_all.mat')

% vals = cell(2,2,2); % ste/pyr, M/F, WT/PS19

valsAll = cell(2,3);
for ff = 1:3
    if ff==1
        load('D:\AD_Project\Behavior\Figures\Figures20241218_fullUpdate\all\data\perSuccess_all.mat','data')
    elseif ff==2
        load('D:\AD_Project\Behavior\Figures\Figures20241218_fullUpdate\female\data\perSuccess_female.mat','data')
    elseif ff==3
        load('D:\AD_Project\Behavior\Figures\Figures20241218_fullUpdate\male\data\perSuccess_male.mat','data')
    end

    dataCur = cellfun(@(x) x(:,8:11),data,'UniformOutput',false);
    valsAll(:,ff) = dataCur;
end

valsPlot = cellfun(@(x) mean(x,2,'omitnan'),valsAll,'UniformOutput',false);


%%

meansAll = cellfun(@(x) mean(x,'omitnan'),valsPlot);
matPlot = -diff(meansAll,1);

pAll = zeros(1,3);
for ff = 1:3
    [~,pAll(ff)] = ttest2(valsPlot{1,ff},valsPlot{2,ff});
end

figure; hold on

F = imagesc(matPlot,'YData',[size(matPlot,1) 1]);

col = colorbar();

curCLim = [-50,50];

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
