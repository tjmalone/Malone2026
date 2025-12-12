%% findLocalCellTypes
% Find cue, grid, and speed cells that meet the requirement for given
% sessions and saves a local cell type file. Set up for AD learning
% analysis.
%



%% Initialize

clear; clc; close all;

p1 = 'D:\AD_Project\imagingData\data';
cd(p1)

% load alignment data
load('D:\AD_Project\imagingData\foldersLearning.mat')

reRun = 1;


%% Identify local cue cells

sideName = {'Left','Right','All'};
suffName = {'_dfof','_sig'};
cellTypes = {'cueL','cueR','cueAll','grid'};
nSide = length(sideName);
nSuff = length(suffName);

% cycle through FOV
for ii=1:length(foldersLearning)
    
    % initialize sessions folders
    dfolders = foldersLearning{ii};
    
    % cycle through sessions
    for jj=1:length(dfolders)
        % move to data directory
        disp([num2str(ii) '-' num2str(jj)])
        cd(dfolders{jj});

        curRoi = matfile('allRois.mat');
        nCells = size(curRoi,'roi',3);
        
        for kk = 1:nSuff
            %% Find cue cells

            % move to suffix directory
            p2 = pwd;
            cd(['cueAnalysis' suffName{kk} '\newScoreShuffleTemplate'])
            
            % load cue data
            load('cueCellsAll.mat','cueCellsAll')
            
            % initialize save variables
            localTypeIdx = cell(1,nSide+1);
            localTypeMat = zeros(length(nCells),nSide+1);
            
            % find cue cells among common cells
            for mm = 1:nSide
                cueCur = cueCellsAll.(sideName{mm}).cueCellRealIdx;
                
                localTypeIdx{mm} = cueCur;
                localTypeMat(cueCur,mm) = 1;
            end
            
            cd(p2)
            
            
            %% Find grid cells
            
            % move to suffix directory
            cd(['gridAnalysis' suffName{kk}])
            
            % load grid data
            load('allCellsCorrected.mat','allCellsCorrected')

            % find grid cells among common cells
            gridIdx = allCellsCorrected.gridIdx;
            
            % remove cue cells from grid cell list
            localTypeMat(gridIdx,nSide+1) = 1;
            
            localTypeMat(:,nSide+1) = localTypeMat(:,nSide+1) & ~any(localTypeMat(:,1:nSide),2);
            localTypeIdx{nSide+1} = localTypeMat(:,nSide+1)==1;

            cd(p2)
            
            
            %% Save cell type info
            
            save(['localCellTypes' suffName{kk}],'cellTypes','localTypeIdx','localTypeMat')

            copyfile('D:\AD_Project\imagingData\code\findLocalCellTypes.m','postAnalysis\findLocalCellTypes.m')
        end
    end
end

