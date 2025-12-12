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
suffName = {'_dfof','_sig'};
suffName2 = {'_dfof','_dfof_sig'};
commonTypes = {'cueL','cueR','cueAll','grid','speed'};
nSide = length(sideName);
nSuff = length(suffName);

globalCellTypeMat = cell(length(foldersLearning),1);

% cycle through FOV
for ii=1:length(foldersLearning)
    
    % initialize sessions folders
    dfolders = foldersLearning{ii};
    
    for kk = 1:nSuff
        globalCellTypeMat{ii,kk} = zeros(size(alignsLearning{ii},1),nSide+2,length(dfolders));
    end
    
    % cycle through sessions
    for jj=1:length(dfolders)
        % move to data directory
        disp([num2str(ii) '-' num2str(jj)])
        cd(dfolders{jj});
        
        % aligned cells
        curAligns = alignsLearning{ii}(:,jj);
        
        for kk = 1:nSuff
            %% Find cue cells

            % skip previously run analysis
            if reRun==0 & isfile(['commonCellTypes' suffName{kk} '.mat'])
                load(['commonCellTypes' suffName{kk} '.mat'])
                globalCellTypeMat{ii,kk}(:,:,jj) = commonTypeMat;
                continue
            end

            % move to suffix directory
            p2 = pwd;
            cd(['cueAnalysis' suffName{kk} '\newScoreShuffleTemplate'])
            
            % initialize save variables
            commonTypeIdx = cell(1,nSide+2);
            commonTypeMat = zeros(length(curAligns),nSide+2);
            
            % load cue data
            load('cueCellsAll.mat','cueCellsAll')
            
            % find cue cells among common cells
            for mm = 1:nSide
                [cueCur,cueCommonIdx] = intersect(curAligns,cueCellsAll.(sideName{mm}).cueCellRealIdx);
                
                commonTypeIdx{mm} = cueCur;
                commonTypeMat(cueCommonIdx,mm) = 1;
            end
            
            cd(p2)
            
            
            %% Find grid cells
            
            % move to suffix directory
            cd(['gridAnalysis' suffName{kk}])
            
            % load grid data
            load('allCellsCorrected.mat','allCellsCorrected')

            % find grid cells among common cells
            [~,gridCommonIdx] = intersect(curAligns,allCellsCorrected.gridIdx);
            
            % remove cue cells from grid cell list
            commonTypeMat(gridCommonIdx,nSide+1) = 1;
            
            commonTypeMat(:,nSide+1) = commonTypeMat(:,nSide+1) & ~any(commonTypeMat(:,1:nSide),2);
            commonTypeIdx{nSide+1} = curAligns(commonTypeMat(:,nSide+1)==1);

            cd(p2)
            

            %% Find speed cells

            % load speed cells
            load(['speed' suffName2{kk} '\speedCellsUniThresh.mat'])
            speedIdxCur = [speedCellsUniThresh.speedCellNegt,speedCellsUniThresh.speedCellPost]';
            [speedTrueCommon,speedCommonIdx] = intersect(curAligns,speedIdxCur);

            commonTypeIdx{nSide+2} = speedCommonIdx;
            commonTypeMat(speedCommonIdx,nSide+2) = 1;

            
            %% Save common cell info
            
            globalCellTypeMat{ii,kk}(:,:,jj) = commonTypeMat;
            save(['commonCellTypes' suffName{kk}],'commonTypes','commonTypeIdx','commonTypeMat')
            
        end
    end
end

cd(p1)
save('globalCellTypes.mat','suffName','commonTypes','globalCellTypeMat')


%% Calculate global cell types

clear; close all; clc

p1 = 'D:\AD_Project\imagingData\data';
cd(p1)

% load data
load('D:\AD_Project\imagingData\foldersLearning.mat')
load('globalCellTypes.mat')

sideName = {'Left','Right','All'};
suffName = {'_dfof','_sig'};
nSide = length(sideName);
nSuff = length(suffName);

% day threshold
dayThresh = 6;

% identify cells by type 
globalCellTypeSum = cellfun(@(x) sum(x,3),globalCellTypeMat,'UniformOutput',0);
globalCellTypeThresh = cellfun(@(x) x>=dayThresh,globalCellTypeSum,'UniformOutput',0);

% get cell type indices
globalCellTypeIdx = struct();
globalCellTypeLogical = struct();
for ii = 1:nSuff
    for jj = 1:nSide+2
        globalCellTypeLogical.([commonTypes{jj} suffName{ii}]) = ...
            cellfun(@(x) x(:,jj)==1,globalCellTypeThresh(:,ii),'UniformOutput',0);
        globalCellTypeIdx.([commonTypes{jj} suffName{ii}]) =...
            cellfun(@(x,y) y(x(:,jj)==1,:),globalCellTypeThresh(:,ii),alignsLearning,'UniformOutput',0);
    end
    globalCellTypeLogical.(['others' suffName{ii}]) = ...
        cellfun(@(x) ~any(x==1,2),globalCellTypeThresh(:,ii),'UniformOutput',0);
    globalCellTypeIdx.(['others' suffName{ii}]) = ...
        cellfun(@(x,y) y(~any(x==1,2),:),globalCellTypeThresh(:,ii),alignsLearning,'UniformOutput',0);
end

save('globalCellTypes.mat','suffName','commonTypes','globalCellTypeMat','globalCellTypeIdx','globalCellTypeLogical')

