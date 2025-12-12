%% identifyPyrSteOverlay
% Generates example visualization of pyramidal/stellate cell selection and
% manual/active cell overlaps

clear; close all; clc

% set directory
foldMan = 'D:\AD_Project\imagingData\analysis_ManualSelection';
cd(foldMan);

% load image directory
imDir = 'D:\AD_Project\imagingData\data_FOV\refImages\';

% load mutaul index information
load('infoMatch.mat','infoMatch')

% load data folders
load('D:\AD_Project\imagingData\foldersLearning.mat','foldersLearning')

% set save directories
figPath = 'overlay_figures\';

% load morphology information
load('data_PyrSte.mat','dataPS');

% define threshold parameters
threshCenter = 158;        % threshold for morphology center
threshRange = 20;            % threshold for morphology range
threshDist = 4;             % threshold for cell matching. Selected based on visualization
threshOvr = 0.4;            % threshold for matching area overlap

% define use FOV
useFOV = 32;
useDay = 1;

% define colors
colorSumActive = [0 0 1];         % blue
colorSumManual = [1 0 0];         % red
colorSumManMatched = [1 0 1];     % purple
colorSTE = [217 89 26]/255;                     % cyan
colorPYR = [60 175 200]/255;                 % green

rangeX = 100:370;
rangeY = 90:430;

rng(42)


%% Calculate overlap percentages

% initialize data structures for current parameter set
dataIdxs = struct();


%% Load data and calculate matches

% load active rois
load([foldersLearning{useFOV}{useDay} '\allROIs.mat'])
roiActive = logical(roi);

% load manual rois
load(['FOV_data\allROIsManual' num2str(useFOV,'%02d') '.mat'])
roiManual = logical(allROIsManual);

% current match information
curMatch = infoMatch{useFOV}{useDay};

% apply distance threshold
useIdx = curMatch(:,4)<=threshDist & curMatch(:,5)>=threshOvr;
matchIdxManual = curMatch(useIdx,1);
matchIdxActive = curMatch(useIdx,2);


%% Generate active manual overlay

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

% get unmatched cell indices
diffIdxManual = setdiff(1:size(roiManual,3),(matchIdxManual));
diffIdxActive = setdiff(1:size(roiActive,3),(matchIdxActive));

% choose random subsets of cells
fractionMan = 1;
fractionAct = 1;
randIdxManual = randsample(diffIdxManual,round(fractionMan*length(diffIdxManual)));
randIdxActive = randsample(diffIdxActive,round(fractionAct*length(diffIdxActive)));

% define plot matrices
borderSumManual = max(borderManual(:,:,randIdxManual),[],3);
borderSumActive = max(borderActive(:,:,randIdxActive),[],3);
borderSumManMatched = max(borderManual(:,:,matchIdxManual),[],3);

% add active cell boundaries to RBG plot
boundaryMaskA = repmat(borderSumActive == 1, [1 1 3]);
colorArrayA = reshape(colorSumActive, [1 1 3]);

% add all manual cell boundaries to RBG plot
boundaryMaskB = repmat(borderSumManual == 1, [1 1 3]);
colorArrayB = reshape(colorSumManual, [1 1 3]);

% add matched manual cell boundaries to RBG plot
boundaryMaskC = repmat(borderSumManMatched == 1, [1 1 3]);
colorArrayC = reshape(colorSumManMatched, [1 1 3]);


% Plot overlay on original image

% load image and make greyscale RBG
peakValue = 1;
imBase = loadtiff([imDir 'refImage' num2str(useFOV) '.tif']);
imScale = imBase/max(imBase,[],'all')*1.4;
imBright = imScale/peakValue;
imBright(imBright>1) = 1;
imOverlay = repmat(imBright,[1 1 3]);

% figure; imshow(imOverlay)

% imOverlay(boundaryMaskA) = repmat(colorArrayA, sum(boundaryMaskA(:))/3, 1);
% imOverlay(boundaryMaskB) = repmat(colorArrayB, sum(boundaryMaskB(:))/3, 1);
% imOverlay(boundaryMaskC) = repmat(colorArrayC, sum(boundaryMaskC(:))/3, 1);

figure; imshow(imOverlay(rangeX,rangeY,:))
% figure; imshow(imOverlay)


%% Generate morphology overlay

usePS = dataPS.allArea(dataPS.fovIdx==useFOV);

combinedIdx = sort([matchIdxManual; randIdxManual']);
combinedArea = usePS(combinedIdx);

isSTE = combinedIdx(combinedArea>=threshCenter-threshRange);
isPYR = combinedIdx(combinedArea<threshCenter+threshRange);

% define plot matrices
borderSumSTE = max(borderManual(:,:,isSTE),[],3);
borderSumPYR = max(borderManual(:,:,isPYR),[],3);

% add stellate cell boundaries to RBG plot
boundaryMaskB = repmat(borderSumSTE == 1, [1 1 3]);
colorArrayB = reshape(colorSTE, [1 1 3]);

% add pyramidal cell boundaries to RBG plot
boundaryMaskC = repmat(borderSumPYR == 1, [1 1 3]);
colorArrayC = reshape(colorPYR, [1 1 3]);

% load image and make greyscale RBG
imOverlay2 = repmat(imBright,[1 1 3]);

% % overlay ROIs
imOverlay2(boundaryMaskB) = repmat(colorArrayB, sum(boundaryMaskB(:))/3, 1);
imOverlay2(boundaryMaskC) = repmat(colorArrayC, sum(boundaryMaskC(:))/3, 1);

figure; imshow(imOverlay2(rangeX,rangeY,:))
% figure; imshow(imOverlay2)




