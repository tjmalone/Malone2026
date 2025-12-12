%% identifyPyrSteOverlay
% Generates example visualization of pyramidal/stellate cell selection and
% manual/active cell overlaps

clear; close all; clc

% set directory
foldMan = 'D:\AD_Project\imagingData\analysis_ManualSelection\HistologyValidation\';
cd(foldMan);

% get images
imDir = 'TM240316-J\';

% set save directories
figPath = 'overlay_figures\';

% define rois and images
filesIm = {'calbindin_2.tif','reelin_2.tif','gcamp_2.tif'};
filesRoi = {'RoiSet_calbindin_cut.zip','RoiSet_reelin_cut.zip'};
nIms = length(filesIm);
nRoi = length(filesRoi);

% define brightness correction
brCorrect = [0.6 0.5 1.1];

% define colors (pyramidal, stellate)
colors2 = {[60 175 200]/255;[217 89 26]/255};


%% Load roi overlays

% define image size
baseIm = loadtiff([foldMan imDir filesIm{1}]);
rows = size(baseIm,1);
columns = size(baseIm,2);

% define dilation strel
se = strel('disk',2); % adjust radius for thickness

boundaryMask = cell(1,nRoi);
colorMask = cell(1,nRoi);
for rr = 1:nRoi
    % extract roi info
    strROIFilename = [foldMan imDir filesRoi{rr}];
    allROIsManual = extractImageJROI(rows,columns,strROIFilename);

    % dilate rois and calculate border
    borderDilated = imdilate(allROIsManual,se)-allROIsManual;

    % combine roi borders
    borderSum = max(borderDilated,[],3);

    % creat roi border mask
    boundaryMask{rr} = repmat(borderSum == 1, [1 1 3]);
    colorMask{rr} = reshape(colors2{rr}, [1 1 3]);
end


%% Plot histology overlays

% define plot range
rangeX = 200:775;
rangeY = 110:550;

figure
tiledlayout(1,nIms)

for ii = 1:nIms
    % load image and make grayscale RBG
    imBase = double(loadtiff([foldMan imDir filesIm{ii}]));
    imScale = imBase/max(imBase,[],'all')*1.4;
    imBright = imScale/brCorrect(ii);
    imBright(imBright>1) = 1;

    % generate RBG image and duplicate for overlay
    imRGB_base = repmat(imBright,[1 1 3]);
    imRGB_overlay = imRGB_base;

    if ii<nIms
        imRGB_overlay(boundaryMask{ii}) = repmat(colorMask{ii},sum(boundaryMask{ii}(:))/3, 1);
    else
        for rr = 1:nRoi
            imRGB_overlay(boundaryMask{rr}) = repmat(colorMask{rr},sum(boundaryMask{rr}(:))/3, 1);
        end
    end

    % figure
    % plot overlay image
    nexttile(ii); hold on
    imshow(imRGB_overlay(rangeX,rangeY,:))
end

