%% heatmap_mean.m
% generate mean activity heatmap

clear; clc; close all

p1 = 'D:\AD_Project\imagingData';
cd(p1)

% initialize raw data struct
fullData = struct();

% load mean dfof heatmap
load('Figures/Figures_current/timecourse/data/dfofSig_common-cell.mat','mapData')
fullData.dfofPre = mapData.Pre;
fullData.dfofPost = mapData.Post;

% load spatial selectivity heatmap
load('Figures/Figures_current/timecourse/data/spatialSelectivity_common-cell.mat','mapData')
fullData.SpatialSelectivity = mapData.All;

% load absolute speed score heatmap
load('Figures/Figures_current/timecourse/data/speedScore_common-cell.mat','mapData')
fullData.SpeedScore = mapData.All;

% load within-day consistency heatmaps
load('Figures/Figures_current/timecourse/data/RBR_common-cell.mat','mapData')
fullData.InterLap_Post = mapData.Post;
fullData.InterLap_Diff = mapData.Diff;

% load inter-day consistency heatmaps
load('Figures/Figures_current/timecourse/data/inter_common-cell.mat','mapData')
fullData.InterDayAll = mapData.All;
fullData.InterDayPost = mapData.Post;

% load inter-mouse consistency heatmaps
load('data/interMouseHeat.mat','valsMap')
fullData.InterMouse_Post = valsMap;

% load dfof decoding heatmaps
load('Figures/Figures_current/timecourse/data/dfofIO_common-cell.mat','mapData')
fullData.dfofIO = mapData.All;

% define deficit directionality (1 = low is good)
mapDir = [-1, 1, -1, -1, -1, -1, -1, -1, -1, 1];

% set save folder name
svFile = [p1 '/Figures/Figures' char(datetime('today','format','yyyyMMdd')) '/heatmaps'];
if ~isfolder(svFile)
    mkdir(svFile);
    mkdir([svFile '/data']);
end


%% Normalize individual heatmaps

% set field parameters
fldNames = fieldnames(fullData);
nFields = length(fldNames);
useFields = 1:10;

% define x axis parameters
useSex = 3:-1:2;
labsSex = {'male','female'};
nSex = length(useSex);

% define y axis parameters
useMorph = 2:3;
labsMorph = {'ste','pyr'};
nMorph = length(useMorph);

% initialize normalized data struct
normData = zeros(nMorph,nSex,nFields);

% normalize each heatmap
for ff = 1:nFields
    % select current raw data
    curRaw = fullData.(fldNames{ff})(useSex,useMorph,:);
    curMean = cellfun(@(x) mean(x,'omitnan'),curRaw);

    % calculate genotype difference
    curDiff = diff(curMean,[],3)'*mapDir(ff);

    % normalize heatmap
    curNorm = curDiff/max(abs(curDiff),[],'all');

    % store normalized
    normData(:,:,ff) = curNorm;
end

% calculate averaged heatmap
normMean = mean(normData(:,:,useFields),3);


%% Plot averaged heatmap

% plot heatmap
figure; hold on
imagesc(normMean,'YData',[size(normMean,1) 1]);

% set color map
clim([-1 1])
curCMap = customcolormap([0 0.25 0.5 0.75 1],[0.75 0 0.75; 1 0.5 1; 1 1 1; 0.5 1 0.5; 0 0.75 0]);
colormap(curCMap);

% set colorbar label
col = colorbar();
ylabel(col,'AD deficit')

% set plot limits
axis('equal')
xlim([0.5 nSex+0.5])
ylim([0.5 nMorph+0.5])
set(gca,'FontSize',16)
set(col,'FontSize',16)

% set x axis labels
set(gca,'XAxisLocation','top');
xticks(1:nSex);
xticklabels(labsSex);
xlabel('Sex')

% set y axis labels
yticks(1:nMorph);
yticklabels(fliplr(labsMorph));
set(gca,'YTickLabelRotation',90);
ylabel('Morphology');

% define labels
title('Averaged Activity Heatmap')

% save figure and data
savefig([svFile '/averagedHeatmap.fig'])
save([svFile '/data/averagedHeatmap.mat'],'normData','normMean',...
    'labsSex','labsMorph','fldNames')

