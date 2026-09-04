%%

clear; close all; clc

p1 = 'Z:\labMembers\TM\AD_Project\imagingData\data_FOV\refImages\highConfidenceFOVs\Final';
cd(p1)

% image raw image and data path
imPath = 'Z:\labMembers\TM\AD_Project\imagingData\data_FOV\refImages\highConfidenceFOVs\Final\';
dataPath = imPath;

useNums = [9,13,14,15,17,19,20,21,22,23,24,27,29,30,31,32,35,37,38,39,40];
nFile = length(useNums);

pixelSize = 750/512; %for 750um FOV under  512x512 pixels
areaWidth = [6,8,10,12,14,16];
buff = 10;

allLong = [];
allShort = [];
allArea = [];


%%

for ii = 1:nFile
    numCur = num2str(useNums(ii));
    numCur2 = num2str(useNums(ii),'%02d');

    % load image
    A = loadtiff([imPath 'refImage' numCur '.tif']);

    if ~isfile([dataPath 'params' numCur2 '.mat'])
        % define image size
        rows = size(A,1);
        columns = size(A,2);

        % extract roi info
        d = dir([dataPath 'RoiSet' numCur2 '*.zip']);
        strROIArchiveFilename = [d(end).folder '\' d(end).name];
        [ allROIsManual ] = extractImageJROI(rows,columns,strROIArchiveFilename );

        % calculate and plot ROI size info
        params = cellParams2D(allROIsManual,pixelSize);
        figure; tiledlayout(1,length(areaWidth))
        for jj = 1:length(areaWidth)
            nexttile()
            cellParams2D_plotArea(params.areaum,areaWidth(jj));
        end

        % save individual size information
        saveas(gcf,[dataPath 'cluster' numCur2 '.fig']);

        % set identification threshold
        splitThr = input('Select bimodal minimum: ');
        ksSize = input('Select ks width: ');

        % save parameters
        save([dataPath 'params' numCur2 '.mat'],'allROIsManual','params','splitThr','ksSize');
    else
        % load parameters
        load([dataPath 'params' numCur2 '.mat'])
    end

    figure; tiledlayout(1,2)
    sgtitle(['FOV ' numCur])
    nexttile()
    cellParams2D_plotArea(params.areaum,ksSize);
    xline(splitThr,'r')
    xline(splitThr-buff,'r:')
    xline(splitThr+buff,'r:')


    %% Plot Cell Types

    cutOff = splitThr;

    % get pyramidal cell boundaries
    idxPyr = params.areaum<cutOff-buff;

    % get stellate cell boundaries
    idxSte = params.areaum>cutOff+buff;

    % plot overlay
    nexttile()
    plotSelectionOverlayTransparent(A,allROIsManual,idxPyr,idxSte)

    saveas(gcf,[dataPath 'overlay' numCur2 '.fig']);

    close all; clear splitThr
end

