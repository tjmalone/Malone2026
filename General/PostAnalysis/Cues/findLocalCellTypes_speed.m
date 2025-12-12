%% findLocalCellTypes
% Find cue, grid, and speed cells that meet the requirement for given
% sessions and saves a local cell type file. Set up for AD learning
% analysis.
%



%% Initialize

clear; clc; close all;

p1 = 'D:\AD_Project\imagingData';
cd(p1)

% load alignment data
load('foldersLearning.mat')


%% Identify local cell types

suffName = {'_dfof','_sig'};
suffName2 = {'_dfof','_dfof_sig'};

cellTypes = {'speedPos','speedNeg'};
nSuff = length(suffName);

% cycle through FOV
for ii=1:length(foldersLearning)
    
    % initialize sessions folders
    dfolders = foldersLearning{ii};
    
    % cycle through sessions
    for jj=1:length(dfolders)
        % move to data directory
        disp([num2str(ii) '-' num2str(jj)])

        cd(dfolders{jj})

        curRoi = matfile('allROIs.mat');
        nCells = size(curRoi,'roi',3);
        
        for kk = 1:nSuff
            %% Find speed cells

            localTypeMat = zeros(nCells,2);
            localTypeIdx = cell(1,2);

            % load speed cells
            load(['speed' suffName2{kk} '/speedCellsUniThresh.mat'])
            speedIdxPos = speedCellsUniThresh.speedCellPost;
            speedIdxNeg = speedCellsUniThresh.speedCellNegt;

            % make cell type matrix
            localTypeMat(speedIdxPos,1) = 1;
            localTypeMat(speedIdxNeg,2) = 1;

            % save cell type indices
            localTypeIdx{1} = speedIdxPos;
            localTypeIdx{2} = speedIdxNeg;


            %% Save cell type info
            
            save(['localCellTypes_speed' suffName{kk}],'cellTypes','localTypeIdx','localTypeMat')

            copyfile('D:/AnalysisCode/PostAnalysis/Cues/findLocalCellTypes.m','postAnalysis\findLocalCellTypes.m')
        end
    end
end

