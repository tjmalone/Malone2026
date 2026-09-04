
%% Set inputs

clear; close all; clc

% p1 = 'Z:\labMembers\TM\AD_Project\Revision\Histology\cellpose_test';
% cd(p1)

allPref = 'TM240316-K_HV_';
baseType = '_Ch2';
probeTypes = {'_Ch3','_Ch1'};
probeVal = [-1 1];
probeZ = 3:5;
baseZ = probeZ(2);
nProbe = length(probeTypes);
nZ = length(probeZ);

% define thresholds
ovrThresh = 0.5;            % minimum overlap with MEC region
pixelSize = 850.19/1024;    % conversion from pixels to distance
threshDist = 4;             % threshold for cell matching.
threshOvr = 0.4;            % threshold for matching area overlap


%% Process reference selections

% define image size
A = loadtiff([allPref num2str(baseZ) baseType '.tif']);
rows = size(A,1);
cols = size(A,2);

% load MEC layer 2 region
MEC = extractImageJROI(rows,cols,'MEC.roi');

% load base rois
% baseFile = [allPref num2str(baseZ) baseType '_rois.zip'];
baseFile = 'base_rois_manual_2.zip';
baseSum = extractImageJROI_RoiSum(rows,cols,baseFile);

% define rois in layer 2 region and get roi properties
baseProps = findOverlap_v3(MEC,baseSum,ovrThresh);
baseList = {baseProps(:).PixelIdxList}';
baseCentroid = cat(1,baseProps(:).Centroid);


%% Process probe selections

% initialize probe roi storage
probeList = cell(1,nProbe);
probeCentroid = cell(1,nProbe);
for pp = 1:nProbe
    curProbe = probeTypes{pp};

    for zz = 1:nZ
        curZ = num2str(probeZ(zz));

        % current roi file
        probeFileCur = [allPref curZ curProbe '_rois.zip'];

        % load current rois
        curSum = extractImageJROI_RoiSum(rows,cols,probeFileCur);

        % get current overlap rois
        curProps = findOverlap_v3(MEC,curSum,ovrThresh);

        % store rois
        probeList{pp} = [probeList{pp};{curProps(:).PixelIdxList}'];
        probeCentroid{pp} = cat(1,probeCentroid{pp},cat(1,curProps(:).Centroid));
    end
end


%% Find matches


matchBase = cell(1,2);
for pp = 1:nProbe
    % calculate all pairwise distances
    PDist = pdist2(baseCentroid,probeCentroid{pp});

    % find best match to base cells
    [~,minCol] = min(PDist,[],2);    % mimimum column per row
    [minVals,minRow] = min(PDist,[],1);     % mimimum row per column

    % find mutual minimums
    mutMin = minCol(minRow)';
    mutIdxProbe = mutMin(mutMin==1:length(mutMin))';
    mutBase = minRow(mutIdxProbe)';
    mutMinDists = minVals(mutIdxProbe)'*pixelSize;

    % calculate matched rois overlap
    nMatch = length(mutMinDists);
    mutOvr = zeros(nMatch,1);
    for kk = 1:nMatch
        idxIntersect = length(intersect(probeList{pp}{mutIdxProbe(kk)},baseList{mutBase(kk)}));
        idxUnion = length(union(probeList{pp}{mutIdxProbe(kk)},baseList{mutBase(kk)}));
        mutOvr(kk) = idxIntersect/idxUnion;
    end

    % apply thresholds
    keepRois = mutMinDists<=threshDist & mutOvr>=threshOvr;

    % store results
    matchBase{pp} = mutBase(keepRois);
end

% find dual overlap cells
idxDual = cell(1,2);
[~,idxDual{1},idxDual{2}] = intersect(matchBase{1},matchBase{2});
for pp = 1:nProbe
    matchBase{pp}(idxDual{pp}) = [];
end


%% Generate final rois with labels

% generate final cell list
finalBase = cat(1,matchBase{:});

% generate cell identity reference list 
probeN = cellfun(@length,matchBase);
finalIDs = [probeVal(1)*ones(probeN(1),1);probeVal(2)*ones(probeN(2),1)];

% reconstruct ROI sums
finalROIs = false(rows,cols,length(finalBase));
for rr = 1:length(finalBase)
    curRoi = zeros(rows,cols);
    curRoi(baseList{finalBase(rr)}) = 1;
    finalROIs(:,:,rr) = curRoi;
end

% generate output struct
roiData = struct();
roiData.idx = finalBase;
roiData.identity = finalIDs;
roiData.ROIs = finalROIs;
roiData.centroid = baseCentroid(finalBase,:);
roiData.list = baseList(finalBase);

% calculate areas
finalAreas = cellfun(@length,roiData.list)*pixelSize^2;
roiData.area = finalAreas;

save('roiData.mat','roiData')


%% Process and plot size histogram

figure; hold on

xVals = 0:300;
wdth = 12;

% plot combined curve
ksdensity(finalAreas,xVals,'Bandwidth',wdth)

% plot pyramidal curve
ksdensity(finalAreas(roiData.identity==-1),xVals,'Bandwidth',wdth)

% plot stellate curve
ksdensity(finalAreas(roiData.identity==1),xVals,'Bandwidth',wdth)


%%

figure; tiledlayout(1,2)

nexttile()
imagesc(sum(finalROIs(:,:,roiData.identity==-1),3))

nexttile()
imagesc(sum(finalROIs(:,:,roiData.identity==1),3))


return
%%

outDir = {'match\calbindinMatch.zip','match\reelinMatch.zip'};
matchFold = 'match\';
if ~isfolder(matchFold)
    mkdir(matchFold)
end

for ii = 1:2
    rois = finalROIs(:,:,roiData.identity==probeVal(ii));
    for rr = 1:size(rois,3)
        B = bwboundaries(rois(:,:,rr), 'noholes');
        writeImageJROI(B{1}, 3, 0, 0, 0, matchFold);
    end

    zip(outDir{ii}, fullfile(matchFold,'*.roi'));
    delete(fullfile(matchFold,'*.roi'));
end

