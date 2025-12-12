%% identifyPyrSteVisualize
% Allows visualization of a range of overlap criteria to allow selection of
% parameters that minimize false negative/positive overlaps. Parameters are
% centroid distance and overlap percentage.
%

% clear; close all; clc

% set directory
foldMan = 'D:\AD_Project\imagingData\analysis_ManualSelection';
cd(foldMan);

% load image directory
imDir = 'D:\AD_Project\imagingData\data_FOV\refImages\';

% load mutaul index information
load('infoMatch.mat','infoMatch')

% load data folders
load('D:\AD_Project\imagingData\foldersLearning.mat','foldersLearning')
useFOV = [13 23 31];
% useFOV = [13];

nFOV = length(useFOV);
useDay = 1;

% set save directories
figPath = 'overlay_figures\';

% define threshold parameters
threshDist = 2:1:7;             % threshold for cell matching
threshDist = 10;             % threshold for cell matching
nD = length(threshDist);

threshOvr = 0.3:0.1:0.5;       % threshold for overlap of cell matching
nO  = length(threshOvr);


%% Generate overlay plots

% loop over overlaps
for gg = 1:nO

    % loop through distances
    for ii = 1:nD
        %% Calculate overlap percentages

        % initialize data structures for current parameter set
        dataIdxs = struct();

        % loop through FOV
        for jj = 1:nFOV
            %% Load data and calculate matches

            % current FOV
            curFOV = useFOV(jj);

            % load active rois
            load([foldersLearning{curFOV}{useDay} '\allROIs.mat'])
            roiActive = logical(roi);

            % load manual rois
            load(['FOV_data\allROIsManual' num2str(curFOV,'%02d') '.mat'])
            roiManual = logical(allROIsManual);

            % current match information
            curMatch = infoMatch{curFOV}{useDay};

            % apply distance threshold
            useIdx_thOvr = curMatch(:,5)>=threshOvr(gg);
            useIdx_thDist = curMatch(:,4)<=threshDist(ii);
            matchIdxManual_miss  = curMatch(~useIdx_thOvr & useIdx_thDist,1);
            matchIdxManual_hit  = curMatch(useIdx_thOvr & useIdx_thDist,1);


            %% Generate roi overlay

            % find active cell borders
            borderActive = false(size(roiActive));
            for kk = 1:size(roiActive,3)
                borderActive(:,:,kk) = bwperim(roiActive(:,:,kk));
            end

            % find manual cell borders
            borderManual = false(size(roiManual));
            for kk = 1:size(roiManual,3)
                borderManual(:,:,kk) = bwperim(roiManual(:,:,kk));
            end

            % define plot matrices
            sumActive = double(max(roiActive,[],3));
            borderSumActive = max(borderActive,[],3);
            borderSumManual = max(borderManual,[],3);
            borderSumManMatched_miss = max(borderManual(:,:,matchIdxManual_miss),[],3);
            borderSumManMatched_hit = max(borderManual(:,:,matchIdxManual_hit),[],3);

            % define colors (active, manual, matched manual miss, match
            % manual hit)
            colors = {[0 0 1],[1 1 0],[.25 .8 .25],[1 0 0]};
            showType = [1 0 1 1];
            
            % define loop data
            dataPlot = {borderSumActive,borderSumManual,borderSumManMatched_miss,borderSumManMatched_hit};

            % generate background RBG plot with acvtive cells as white
            plotData = repmat(sumActive,[1 1 3]);

            for kk = 1:4
                if showType(kk)==0; continue; end

                % add active cell boundaries to RBG plot
                boundaryMask = repmat(dataPlot{kk} == 1, [1 1 3]);
                colorArray = reshape(colors{kk}, [1 1 3]);
                plotData(boundaryMask) = repmat(colorArray, sum(boundaryMask(:))/3, 1);
            end

            % plot roi overlay
            figure
            imshow(plotData)
            title(['Active/Manual Overlay: d=' num2str(threshDist(ii)) ', o=' num2str(threshOvr(gg))])

            % save figure
            % savefig(gcf,[figPath 'overlay_FOV' num2str(curFOV) '_d' num2str(threshDist(ii)) '_o' replace(num2str(threshOvr(gg)),'.','-') '.fig'])


            %% Plot overlay on original image

            % % load image and make greyscale RBG
            % imBase = loadtiff([imDir 'refImage' num2str(curFOV) '.tif']);
            % imScale = imBase/max(imBase,[],'all')*1.2;
            % imOverlay = repmat(imScale,[1 1 3]);
            % 
            % figure; imshow(imOverlay)
            % 
            % imOverlay(boundaryMask) = repmat(colorArray, sum(boundaryMask(:))/3, 1);
            % imOverlay(boundaryMaskB) = repmat(colorArrayB, sum(boundaryMaskB(:))/3, 1);
            % imOverlay(boundaryMaskC) = repmat(colorArrayC, sum(boundaryMaskC(:))/3, 1);
            % 
            % figure; imshow(imOverlay)
        end

    end
end
