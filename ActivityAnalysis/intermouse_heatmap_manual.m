%%

clear; close all; clc
% vals = cell(2,2,2); % ste/pyr, M/F, WT/PS19

valsAll = cell(4,3);
for ff = 1:3
    if ff==1
        load('D:\AD_Project\imagingData\Figures\Figures_current\trajectory\data\interMouseCorr_common-allMorph_Post.mat','corrAll','groupIdxs')
    elseif ff==2
        load('D:\AD_Project\imagingData\Figures\Figures_current\trajectory\data\interMouseCorr_common-ste_Post.mat','corrAll','groupIdxs')
    elseif ff==3
        load('D:\AD_Project\imagingData\Figures\Figures_current\trajectory\data\interMouseCorr_common-pyr_Post.mat','corrAll','groupIdxs')
    end

    corrAll(eye(size(corrAll,1))==1) = NaN;

    for ii = 1:4
        valsAll{ii,ff} = corrAll(groupIdxs(ii)+1:groupIdxs(ii+1),:);
    end
end

valsMean = cellfun(@(x) mean(x,2,'omitnan'),valsAll,'UniformOutput',false);
meansAll = cellfun(@(x) mean(x,'all','omitnan'),valsMean);

meansMap = reshape(meansAll,2,2,[]);
meansMap2 = permute(meansMap,[1 3 2]);
meansMap3 = cat(1,meansMap2(1,:,:),meansMap2);


valsMap = reshape(valsMean,2,2,[]);
valsMap = permute(valsMap,[1 3 2]);
valsMap = cat(1,valsMap(1,:,:),valsMap);

save('data\interMouseHeat.mat','valsMap')

%%

sexIDs = {'allSex','female','male'};
morphTypes = {'allMorph','ste','pyr'};

% make heat map
% h = heatmapGrid(valsMap,['Sex',sexIDs],['Morphology',morphTypes]);

[h,p] = heatmapGridSplit(valsMap,['Sex',sexIDs],['Morphology',morphTypes],...
    1,{2:3,2:3},{1,2:3});

% save figure
% savefig([svFile svLabelCur '_' catSvLabels{flc}])

% initialize stats parameters and array
outStats = {};
nUnits = 'mice';
testName = 'two-tailed unpaired Students t-test';
testPair = 0;
testMC = 0;
testLimitP = 0;

idxX = [2:3;2:3];
idxY = 2:3;
cats = {'stellate','pyramidal'};

for ii = 1:2
    curData = squeeze(valsMap(idxX(ii,:),idxY(ii),:));

    outStats(end+1,:) = [cats{ii} ttestEffectSize(...
        curData(:,1),curData(:,2),testName,nUnits,testPair,testMC,testLimitP)];
end


%%

matPlot = (meansAll(1:2,:)-meansAll(3:4,:))';
pAll = zeros(2,2);

for ii = 1:2
    for jj = 1:2
        [~,pAll(jj,ii)] = ttest2(valsAll{ii,jj}(:),valsAll{ii+2,jj}(:));
    end
end

figure; hold on

F = imagesc(matPlot,'YData',[size(matPlot,1) 1]);

col = colorbar();

curCLim = [-0.1,0.1];

curCMap = customcolormap([0 0.25 0.5 0.75 1],[0.75 0 0.75; 1 0.5 1; 1 1 1; 0.5 1 0.5; 0 0.75 0]);

ylabel(col,'AD deficit')

% set colors
clim(curCLim)
colormap(curCMap);

% set plot sizes
nX = 2;
nY = 2;
curXRange = [0.5 nX+0.5];
curYRange = [0.5 nY+0.5];
labelX = 'Sex';
labelY = 'Morph';

% set plot limits
axis('equal')
xlim(curXRange)
ylim(curYRange)
set(gca,'FontSize',16)
set(col,'FontSize',16)

% set x axis labels
set(gca,'XAxisLocation','top');
set(gca,'XTick',1:nX);
set(gca,'XTickLabel',{'Female','Male'});
xlabel('Sex')

% set y axis labels
set(gca,'YTick',1:nY);
set(gca,'YTickLabel',{'ste','pyr'});
set(gca,'YTickLabelRotation',90);
ylabel('Morph')
