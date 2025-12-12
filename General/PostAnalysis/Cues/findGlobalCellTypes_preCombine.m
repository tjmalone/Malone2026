%% findGlobalCellTypes
% Find cue, grid, and speed cells among a a set of common cells that meet
% the requirement in a sufficient number of sessions. Set up for AD
% learning analysis.
%


%% Initialize

clear; clc; close all;

p1 = 'D:\AD_Project\imagingData\data';
cd(p1)

% load alignment data
load('D:\AD_Project\imagingData\foldersLearning.mat')

reRun = 1;


%% Identify global cue cells

sideName = {'Left','Right','All'};
suffName = '_sig';
suffName2 = '_dfof_sig';
commonTypes = {'cueL','cueR','cueAll','grid','speedPos','speedNeg'};
nSide = length(sideName);
nTypes = length(commonTypes);

globalCellTypeMat = cell(length(foldersLearning),1);

% cycle through FOV
for ii=1:length(foldersLearning)

    % initialize sessions folders
    dfolders = foldersLearning{ii};

    globalCellTypeMat{ii} = zeros(size(alignsLearning{ii},1),nTypes,length(dfolders));

    % cycle through sessions
    for jj=1:length(dfolders)
        % move to data directory
        disp([num2str(ii) '-' num2str(jj)])
        cd(dfolders{jj});

        % aligned cells
        curAligns = alignsLearning{ii}(:,jj);

        % initialize save variables
        commonTypeIdx = cell(1,nTypes);
        commonTypeMat = zeros(length(curAligns),nTypes);


        %% Find cue cells

        % load cue data
        load(['cueAnalysis' suffName '\newScoreShuffleTemplate\cueCellsAll.mat'],'cueCellsAll')

        % find cue cells among common cells
        for mm = 1:nSide
            [cueCur,cueCommonIdx] = intersect(curAligns,cueCellsAll.(sideName{mm}).cueCellRealIdx);

            commonTypeIdx{mm} = cueCur;
            commonTypeMat(cueCommonIdx,mm) = 1;
        end


        %% Find grid cells

        % load grid data
        load(['gridAnalysis' suffName '\allCellsCorrected.mat'],'allCellsCorrected')

        % find grid cells among common cells
        [~,gridCommonIdx] = intersect(curAligns,allCellsCorrected.gridIdx);

        % remove cue cells from grid cell list
        commonTypeMat(gridCommonIdx,nSide+1) = 1;

        commonTypeMat(:,nSide+1) = commonTypeMat(:,nSide+1) & ~any(commonTypeMat(:,1:nSide),2);
        commonTypeIdx{nSide+1} = curAligns(commonTypeMat(:,nSide+1)==1);


        %% Find speed cells

        % load speed cells
        load(['speed' suffName2 '\speedCellsUniThresh.mat'])

        % define positive speed cells
        [~,speedPCommonIdx] = intersect(curAligns,speedCellsUniThresh.speedCellPost);
        commonTypeIdx{nSide+2} = speedPCommonIdx;
        commonTypeMat(speedPCommonIdx,nSide+2) = 1;

        % define negative speed cells
        [~,speedNCommonIdx] = intersect(curAligns,speedCellsUniThresh.speedCellNegt);
        commonTypeIdx{nSide+3} = speedNCommonIdx;
        commonTypeMat(speedNCommonIdx,nSide+3) = 1;


        %% Save common cell info

        globalCellTypeMat{ii}(:,:,jj) = commonTypeMat;
    end
end

cd(p1)
save('globalCellTypes_PreComb.mat','suffName','commonTypes','globalCellTypeMat')


%% Calculate global cell types

clear; close all; clc

p1 = 'D:\AD_Project\imagingData\data';
cd(p1)

% load data
load('D:\AD_Project\imagingData\foldersLearning.mat','alignsLearning')
load('globalCellTypes_PreComb.mat')

% day threshold
dayThresh = 6;

% define combination categories
catIdx = {[1 2 3],[5 6]};
catsN = length(catIdx);

% define others categories
OIdx = {1:6,[4 7 8]};
ONames = {'OthersInd','OthersComb'};
ON = length(OIdx);

% define category names
commonTypes = [commonTypes {'cueComb','speedComb','OthersInd','OthersComb'}];
nTypes = length(commonTypes);

% sum categorization across days
globalCellTypeSum = cellfun(@(x) sum(x,3),globalCellTypeMat,'UniformOutput',0);

% add combination categories
for cc = 1:catsN
    globalCellTypeSum = cellfun(@(x) cat(2,x,sum(x(:,catIdx{cc}),2)),globalCellTypeSum,'UniformOutput',0);
end

% apply day threshold
globalCellTypeThresh = cellfun(@(x) x>=dayThresh,globalCellTypeSum,'UniformOutput',0);

% get cell type indices
globalCellTypeIdx = struct();
globalCellTypeLogical = struct();

% store cell type logical and indices
for jj = 1:nTypes
    globalCellTypeLogical.(commonTypes{jj}) = ...
        cellfun(@(x) x(:,jj)==1,globalCellTypeThresh,'UniformOutput',0);
    globalCellTypeIdx.(commonTypes{jj}) =...
        cellfun(@(x,y) y(x(:,jj)==1,:),globalCellTypeThresh,alignsLearning,'UniformOutput',0);
end

% store others logical and indices
for nn = 1:ON
    globalCellTypeLogical.(ONames{nn}) = ...
        cellfun(@(x) ~any(x(:,OIdx{nn})==1,2),globalCellTypeThresh,'UniformOutput',0);
    globalCellTypeIdx.(ONames{nn}) = ...
        cellfun(@(x,y) y(~any(x(:,OIdx{nn})==1,2),:),globalCellTypeThresh,alignsLearning,'UniformOutput',0);
end

save('globalCellTypes_PreComb.mat','suffName','commonTypes',...
    'globalCellTypeMat','globalCellTypeIdx','globalCellTypeLogical')

