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
foldMan = 'D:\AD_Project\imagingData\analysis_ManualSelection';
cd(foldMan);

% load data folders
load('D:\AD_Project\imagingData\foldersLearning.mat')
useDays = 1:11;
nDays = length(useDays);
nFOV = length(foldersLearning);

% find alignment folders
foldsAlign = findSubF('AlignmentFinal',2,'D:\AD_Project\imagingData\',0);

% data path
dataPath = 'FOV_data\';

% find manual roi files
fileMan = dir([dataPath 'allROIsManual*.mat']);

% define use files
numUse = 1:42;
nFile = length(numUse);

% whether to rerun morphology and activity collation
newMan = 0;
newActive = 0;


%% Extract morphology parameters


if newMan || ~isfile('data_PyrSte.mat')
    allLong = [];
    allShort = [];
    allArea = [];
    fovIdx = [];

    for ii = 1:nFile
        numCur = num2str(numUse(ii),'%02d');

        % save individual size information
        load([dataPath 'params' numCur '.mat'],'params');

        % store global parameters information
        allLong = [allLong; params.longAxisum];
        allShort = [allShort; params.shortAxisum];
        allArea = [allArea; params.areaum];
        fovIdx = [fovIdx; ii*ones(size(params.longAxisum))];
    end

    % save global data
    dataPS = struct();
    dataPS.allLong = allLong;
    dataPS.allShort = allShort;
    dataPS.allArea = allArea;
    dataPS.fovIdx = fovIdx;
    save('data_PyrSte.mat','dataPS');
end


%% Calculate index overlaps

infoMatch = cell(nFile,1);
nCells = struct();
nCells.manual = zeros(nFile,1);
nCells.active = zeros(nFile,nDays);
nCells.common = zeros(nFile,1);

cellOvr = cell(nFile,nDays);

% loop through FOV
for ii = 1:nFile
    %% Get activity roi projections

    disp(numUse(ii))

    % get alignment projections
    cd(foldsAlign{numUse(ii)})

    % calculate or load aligned active cell centroids
    if newActive || ~isfile('alignmentProjections.mat')
        [projection_sum,SFC] = alignmentCombinedProjection(useDays,1);

        % load reference rois
        cd(foldersLearning{numUse(ii)}{1})
        load('allROIs.mat')

        % add reference projection
        SFC{1} = roi;

        % convert projections to reference size
        nX = 1:size(SFC{1},1);
        nY = 1:size(SFC{1},2);
        for jj = 2:nDays
            SFC{useDays(jj)} = SFC{useDays(jj)}(nX,nY,:);
        end

        % calculate projection sums for plotting
        sumSFC = cellfun(@(x) sum(x,3),SFC,'UniformOutput',false);

        % calculate properties of activity rois
        centroidAct = cell(1,nDays);
        pixelAct = cell(1,nDays);

        % loop through all days
        for jj = 1:nDays
            centroidAct{jj} = zeros(size(SFC{jj},3),2);
            pixelAct{jj} = cell(size(SFC{jj},3),1);

            % get centroids for all roi
            for kk = 1:size(SFC{jj},3)
                rProps = regionprops(SFC{jj}(:,:,kk),'centroid','PixelIdxList');
                if ~isempty(rProps)
                    centroidAct{jj}(kk,:) = rProps(1).Centroid;
                    pixelAct{jj}{kk} = int64(rProps(1).PixelIdxList);
                else
                    centroidAct{jj}(kk,:) = nan(1,2);
                    pixelAct{jj}{kk} = nan;
                end
            end
        end

        cd(foldsAlign{numUse(ii)})

        % save projection information
        save('alignmentProjections.mat','centroidAct','pixelAct','projection_sum','SFC','-v7.3')

    else
        % load projection information
        load('alignmentProjections.mat','centroidAct','pixelAct')
    end

    cd(foldMan)


    %% Get manual roi footprints

    % load manual rois
    load([fileMan(ii).folder '\' fileMan(ii).name])

    % calculate manual sums for plotting
    sumMan = sum(allROIsManual,3);

    % calculate properties of manual rois
    centroidMan = zeros(size(allROIsManual,3),2);
    pixelMan = cell(size(allROIsManual,3),1);
    for kk = 1:size(allROIsManual,3)
        rProps = regionprops(allROIsManual(:,:,kk),'centroid','PixelIdxList');
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
save('infoMatch.mat','nCells','infoMatch','key')

