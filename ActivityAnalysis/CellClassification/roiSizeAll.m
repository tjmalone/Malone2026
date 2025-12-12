%% roiSizeAll
% Performs initial processing of manually selected ROIs and saves roi
% locations and parameters. Performs initial plotting
%

clear; close all; clc

p1 = 'D:\AD_Project\imagingData\analysis_ManualSelection';
cd(p1)

% image raw image and data path
imPath = 'D:\AD_Project\imagingData\data_FOV\refImages\';
dataPath = 'FOV_data\';

numUse = 1:42;

% update 10,20
nFile = length(numUse);

pixelSize = 750/512; %for 750um FOV under  512x512 pixels

allLong = [];
allShort = [];
allArea = [];

for ii = 1:nFile
    numCur = num2str(numUse(ii));
    numCur2 = num2str(numUse(ii),'%02d');

    % define image size
    A = loadtiff([imPath 'refImage' numCur '.tif']);
    rows = size(A,1);
    columns = size(A,2);

    % extract roi info
    d = dir([dataPath 'RoiSet' numCur2 '*.zip']);
    strROIArchiveFilename = [d(1).folder '\' d(1).name];
    [ allROIsManual ] = extractImageJROI(rows,columns,strROIArchiveFilename );
    save([dataPath 'allROIsManual' numCur2 '.mat'],'allROIsManual');

    % calculate and plot ROI size info
    params = cellParams2D(allROIsManual, pixelSize);
    cellParams2D_plot(params.longAxisum,params.shortAxisum,params.areaum,pixelSize);

    % save individual size information
    saveas(gcf,[dataPath 'cluster' numCur2 '.fig']);
    save([dataPath 'params' numCur2 '.mat'],'params');
    % close all

    % store global parameters information
    allLong = [allLong; params.longAxisum];
    allShort = [allShort; params.shortAxisum];
    allArea = [allArea; params.areaum];
end

% plot combine figure
cellParams2D_plot(allLong,allShort,allArea,pixelSize);

figure;
scatter(allLong./allShort,allArea,'b.')
xlabel('long-short ratio')
ylabel('area')


%%

clear; close all; clc

p1 = 'D:\AD_Project\imagingData\analysis_ManualSelection\';
cd(p1)

% data path
dataPath = 'FOV_data\';

numUse = 1:42;
% numUse = [31 32 37:40];

% update 10,20
nFile = length(numUse);

pixelSize = 750/512; %for 750um FOV under  512x512 pixels

allLong = [];
allShort = [];
allArea = [];

for ii = 1:nFile
    numCur2 = num2str(numUse(ii),'%02d');

    % save individual size information
    load([dataPath 'params' numCur2 '.mat'],'params');

    % store global parameters information
    allLong = [allLong; params.longAxisum];
    allShort = [allShort; params.shortAxisum];
    allArea = [allArea; params.areaum];
end

save('roiSizeFinal.mat','allArea','allShort','allLong')


%%
% plot combine figure
cellParams2D_plot_tight(allLong,allShort,allArea,0.4);

% save global figure
saveas(gcf,'globalCluster_sub.fig');

%%

% plot combine figure
cellParams2D_plot_tight_noScatter(allLong,allShort,allArea,4);

% save global figure
% saveas(gcf,'globalCluster.fig');


% cellParams2D_plot(allLong/pixelSize,allShort/pixelSize,allArea/pixelSize/pixelSize,pixelSize);




