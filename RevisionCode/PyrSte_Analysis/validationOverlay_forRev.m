%% identifyPyrSteOverlay
% Generates example visualization of pyramidal/stellate cell selection and
% manual/active cell overlaps

clear; close all; clc

% set directory
foldMan = 'Z:\labMembers\TM\AD_Project\Revision\Histology\revisionHistology';
cd(foldMan);

% find all validation selections
foldsHV = findSubF('HV',2,[],0);

% use FOV
useFOV = [3:6 11:12 16];
imIdx = 12;
imIdxSub = find(useFOV==imIdx);

% load images
dHV = dir([foldsHV{imIdx} '\HV*.tif']);
imBase =double(loadtiff([foldsHV{imIdx} dHV(2).name]));
rows = size(imBase,1);
columns = size(imBase,2);
nIms = size(imBase,3);

% load true rois
load([foldsHV{imIdx} '\roiData.mat'],'roiData')
try
    ROIs = roiData.ROIs;
catch
    load([foldsHV{imIdx} '\ROI.mat'],'ROIs')
end
identsTrue = roiData.identity;
roiId = [1 -1];
roiCh = {1,[1 2],2};
nCh = length(roiId);

% load putative cell identities
load('Z:\labMembers\TM\AD_Project\Revision\Histology\histologyValidation\histValAccuracy.mat','dataAcc')
identsPut = dataAcc.ind.ident{imIdxSub};

% check roi number
if length(identsPut)~=length(identsTrue)
    error('Cell number mismatch')
end
idents = [identsTrue,identsPut];
nType = size(idents,2);

% set save directories
figPath = 'D:\AD_Project\Revision\PyrSteAnalysis\Results\Validation\';

% define colors (stellate, pyramidal)
colors2 = {[217 89 26]/255;[60 175 200]/255};

% define putative IDs and colors
roiIdPut = [1 0 -1];
colors3 = {[217 89 26]/255;[1 1 1];[60 175 200]/255};
nChPut = length(roiIdPut);
putIm = [1 3];


%% Load true roi overlays

% define dilation strel
se = strel('disk',2); % adjust radius for thickness

% dilate all rois
borderDilated = imdilate(ROIs,se)-ROIs;

boundaryMask = cell(nType,nCh);
colorMask = cell(nType,nCh);
for tt = 1:nType
    for cc = 1:nCh
        % get current indices
        idxCur = idents(:,tt)==roiId(cc);

        % combine roi borders
        borderSum = max(borderDilated(:,:,idxCur),[],3);

        % creat roi border mask
        boundaryMask{tt,cc} = repmat(borderSum == 1, [1 1 3]);
        colorMask{tt,cc} = reshape(colors2{cc}, [1 1 3]);
    end
end


%% Plot true histology overlays

% define plot range
rangeX = 1:2500;
rangeY = 1:1000;

% FOV 3
% rangeX = 350:1000;
% rangeY = 75:475;

% FOV 4
% rangeX = 540:1290;
% rangeY = 200:700;

% FOV 12
rangeX = 450:1100;
rangeY = 100:500;

% get scale bar conversion
pixelSize = 850.19/1024;    % um/pixel
szPixel = length(rangeX); % height in pixels
szUm = szPixel*pixelSize;   % height in um
szFrac = 50/szUm;           % fraction of total length of 50um section

% define brightness correction
brCorrect = [0.4 1 0.7];

for tt = 1%:nType
    figure
    tiledlayout(1,nIms)

    for ii = 1:nIms
        % make grayscale RBG
        imScale = imBase(:,:,ii)/max(imBase(:,:,ii),[],'all')*1.4;
        imBright = imScale/brCorrect(ii);
        imBright(imBright>1) = 1;

        % generate and duplicate RBG image
        imRGB_base = repmat(imBright,[1 1 3]);
        imRGB_overlay = imRGB_base;

        % add roi overlays
        for cc = 1:nCh
            if ismember(cc,roiCh{ii})
                imRGB_overlay(boundaryMask{tt,cc}) = repmat(colorMask{tt,cc},sum(boundaryMask{tt,cc}(:))/3, 1);
            end
        end

        % plot overlay image
        nexttile(ii); hold on
        imshow(imRGB_overlay(rangeX,rangeY,:))
    end
end


%% Load putative roi overlays

% define dilation strel
se = strel('disk',4); % adjust radius for thickness

% dilate all rois
borderDilated = imdilate(ROIs,se)-ROIs;

boundaryMaskPut = cell(nCh,nChPut);
colorMaskPut = cell(nCh,nChPut);
for cc = 1:nCh
    % get true IDs
    idxTrue = identsTrue==roiId(cc);

    for ccID = 1:nChPut
        % get current indices
        idxPut = identsPut==roiIdPut(ccID);

        % combine roi borders
        borderSum = max(borderDilated(:,:,idxPut & idxTrue),[],3);

        % creat roi border mask
        boundaryMaskPut{cc,ccID} = repmat(borderSum == 1, [1 1 3]);
        colorMaskPut{cc,ccID} = reshape(colors3{ccID}, [1 1 3]);
    end
end


%% Plot true histology overlays

brCorrect2 = [0.4 0.7];

figure
tiledlayout(1,nCh)

for cc = 1:nCh
    % make grayscale RBG
    imScale = imBase(:,:,putIm(cc))/max(imBase(:,:,putIm(cc)),[],'all')*1.4;
    imBright = imScale/brCorrect2(cc);
    imBright(imBright>1) = 1;

    % generate and duplicate RBG image
    imRGB_base = repmat(imBright,[1 1 3]);
    imRGB_overlay = imRGB_base;

    % add roi overlays
    for ccID = 1:nChPut
        imRGB_overlay(boundaryMaskPut{cc,ccID})...
            = repmat(colorMaskPut{cc,ccID},sum(boundaryMaskPut{cc,ccID}(:))/3, 1);
    end

    % plot overlay image
    nexttile(cc); hold on
    imshow(imRGB_overlay(rangeX,rangeY,:))
end