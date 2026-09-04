%% identifyPyrSteMatch
% Extracts and combines the morphology parameters from a defined set of
% manually identified cells. Also, identifies identifies the closest manual
% and active cells based on centroid distance. Outputs the manual indices,
% active indices, common cell identity, and distance. Uses the
% "params_##.mat" and "allROIsManual##.mat" files output from roiSizeAll.
% Does not identify pyramidal/stellate cells or apply a distance threshold
% to matched cells.


%% Load folders/file infomation

clear; close all; clc

% set directory
foldMan = 'D:\AD_Project\Revision\PyrSteAnalysis\';
cd(foldMan);

% load data folders
load('D:\AD_Project\imagingData\foldersLearning.mat')
useDays = 1:11;
nDays = length(useDays);

% find alignment folders
foldsAlign = findSubF('AlignmentFinal',3,'D:\AD_Project\imagingData\',0);

% data path
dataPath = 'Z:\labMembers\TM\AD_Project\imagingData\data_FOV\refImages\highConfidenceFOVs\Final\';

% define use files
numUse = 1:42;
nFOV = length(numUse);

% define threshold parameters
threshRange = 10;            % threshold for morphology range
threshDist = 4;             % threshold for cell matching. Selected based on visualization
threshOvr = 0.4;            % threshold for matching area overlap


% whether to rerun morphology and activity collation
newMan = 1;


%% Extract morphology parameters

if newMan || ~isfile('data_PyrSte.mat')
    allArea = cell(nFOV,1);
    allRoiMan = cell(nFOV,1);
    allThresh = nan(nFOV,1);
    fovIdx = [];

    for ii = 1:nFOV
        numCur = num2str(numUse(ii),'%02d');
        fileCur = [dataPath 'params' numCur '.mat'];

        % store individual size information
        if isfile(fileCur)
            load([dataPath 'params' numCur '.mat'],...
                'params','allROIsManual','splitThr');

            % store global parameters information
            allArea{ii} = params.areaum;
            allRoiMan{ii} = allROIsManual;
            allThresh(ii) = splitThr;
        else
            continue
        end
    end

    % save global data
    dataPS = struct();
    dataPS.allArea = allArea;
    dataPS.allRoiMan = allRoiMan;
    dataPS.allThresh = allThresh;
    save('data_PyrSte_rev.mat','dataPS','-v7.3');
end


%% Calculate index overlaps

infoMatch = cell(nFOV,1);
nCells = struct();
nCells.manual = zeros(nFOV,1);
nCells.active = zeros(nFOV,nDays);
nCells.common = zeros(nFOV,1);

cellOvr = cell(nFOV,nDays);

% loop through FOV
for ii = 1:nFOV
    %% Get activity roi projections

    % skip empty FOV
    if isnan(dataPS.allThresh(ii))
        continue
    end
    disp(numUse(ii))

    % load aligned active cell centroids
    load([foldsAlign{numUse(ii)} '\alignmentProjections.mat'],...
        'centroidAct','pixelAct')

    % get current manual footprints
    curROIsManual = dataPS.allRoiMan{ii};

    % calculate properties of manual rois
    centroidMan = zeros(size(curROIsManual,3),2);
    pixelMan = cell(size(curROIsManual,3),1);
    for kk = 1:size(curROIsManual,3)
        rProps = regionprops(curROIsManual(:,:,kk),'centroid','PixelIdxList');
        centroidMan(kk,:) = rProps(1).Centroid;
        pixelMan{kk} = int64(rProps(1).PixelIdxList);
    end


    %% Identify match cells

    infoMatch{ii} = cell(1,nDays);

    for jj = 1:nDays
        % correct for missing imaging days
        trueJ = trueDays(numUse(ii),jj);
        if trueJ>nDays
            continue
        end

        % get active cell centroids
        curAct = centroidAct{jj};

        % calculate all pairwise distances
        curPDist = pdist2(curAct,centroidMan);

        % get minimum distances
        [~,minCol] = min(curPDist,[],2);    % mimimum column per row
        [minVals,minRow] = min(curPDist,[],1);     % mimimum row per column

        % find mutual minimums
        mutMin = minCol(minRow)';
        mutIdxMan = mutMin(mutMin==1:length(mutMin))';
        mutIdxAct = minRow(mutIdxMan)';
        mutMinDists = minVals(mutIdxMan)';

        % find common cell indices
        mutMinCommon = ismember(mutIdxAct,alignsLearning{numUse(ii)}(:,jj));

        nMin = length(mutIdxMan);
        mutOvr = zeros(nMin,1);
        for kk = 1:nMin
            idxIntersect = length(intersect(pixelMan{mutIdxMan(kk)},pixelAct{jj}{mutIdxAct(kk)}));
            idxUnion = length(union(pixelMan{mutIdxMan(kk)},pixelAct{jj}{mutIdxAct(kk)}));
            mutOvr(kk) = idxIntersect/idxUnion;
        end

        % store final matched indices
        infoMatch{ii}{trueJ} = [mutIdxMan,mutIdxAct,mutMinCommon,mutMinDists,mutOvr];
        nCells.active(ii,trueJ) = length(curAct);
    end

    % store total manual cell numbers
    nCells.manual(ii) = size(centroidMan,1);
    nCells.common(ii) = size(alignsLearning{numUse(ii)},1);

end

key = {'Manual ROI Indices','Active ROI Indices','Common Cell State','Matched Distances','Matched Overlaps'};
save('infoMatch_rev.mat','nCells','infoMatch','key')


%% Identify overlap indices

% initialize data structure
morphIdxs = struct();
morphIdxs.allCells.ste = cell(nFOV,nDays);
morphIdxs.allCells.pyr = cell(nFOV,nDays);
morphIdxs.commonCells.ste = cell(nFOV,nDays);
morphIdxs.commonCells.pyr = cell(nFOV,nDays);

% loop through FOV
for ff = 1:nFOV
    % threshold for morphology center
    threshCenter = dataPS.allThresh(ff);
    if isnan(threshCenter)
        continue
    end

    % current morphology edges
    thLow = threshCenter-threshRange;
    thHigh = threshCenter+threshRange;

    % calculate morphology classes
    diamCur = zeros(size(dataPS.allArea{ff},1),1);
    diamCur(dataPS.allArea{ff}<=thLow) = 2;
    diamCur(dataPS.allArea{ff}>=thHigh) = 1;

    % loop through days
    for gg = 1:nDays
        % find cell indices for base categories
        idxsSte = find(diamCur==1);
        idxsPyr = find(diamCur==2);

        for ii = 1:2
            % current match information
            if ii==1
                curMatch = infoMatch{ff}{gg};
            elseif ii==2
                curMatch = infoMatch{ff}{1};
            end

            if ~isempty(curMatch)
                % apply distance threshold
                useIdx = curMatch(:,4)<=threshDist & curMatch(:,5)>=threshOvr;

                % match stellate and pyramidal indices
                isStellate = ismember(curMatch(:,1),idxsSte);
                isPyramidal = ismember(curMatch(:,1),idxsPyr);

                if ii==1
                    % store all cell indices
                    morphIdxs.allCells.ste{ff,gg} = curMatch(useIdx & isStellate,2);
                    morphIdxs.allCells.pyr{ff,gg} = curMatch(useIdx & isPyramidal,2);
                else
                    % get matched common cell indices for day 1
                    idxSte = curMatch(useIdx & isStellate & curMatch(:,3)==1,2);
                    idxPyr = curMatch(useIdx & isPyramidal & curMatch(:,3)==1,2);

                    % identify positions of matched common cells
                    FOVCellsSte = ismember(alignsLearning{ff}(:,1),idxSte);
                    FOVCellsPyr = ismember(alignsLearning{ff}(:,1),idxPyr);

                    % store common cell indices unadjusted for "true" cells
                    morphIdxs.commonCells.ste{ff,gg} = alignsLearning{ff}(FOVCellsSte,gg);
                    morphIdxs.commonCells.pyr{ff,gg} = alignsLearning{ff}(FOVCellsPyr,gg);
                end
            end
        end
    end
end

% save data
save('morphIdxs.mat','morphIdxs')

