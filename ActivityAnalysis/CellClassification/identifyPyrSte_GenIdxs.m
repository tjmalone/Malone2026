%% identifyPyrSte_GenIdxs
% Calculates final active cell morphologies for learning data and exports
% to data folder

clear; close all;clc

% set directory
foldMan = 'D:\AD_Project\imagingData\analysis_ManualSelection';
cd(foldMan);

% load morphology information
load('data_PyrSte.mat','dataPS');

% load mutaul index information
load('infoMatch.mat','nCells','infoMatch')

% load alignment data
load('D:\AD_Project\imagingData\foldersLearning.mat')
nFOV = size(trueDays,1);
nDays = size(trueDays,2);

% set save directories
outPath = 'D:\AD_Project\imagingData\data\morph\';

% define threshold parameters
threshCenter = 158;        % threshold for morphology center
threshRange = 20;            % threshold for morphology range
threshDist = 4;             % threshold for cell matching. Selected based on visualization
threshOvr = 0.4;            % threshold for matching area overlap


%% Define morphology groups

% current morphology edges
thLow = threshCenter-threshRange;
thHigh = threshCenter+threshRange;

% calculate morphology classes
morph = zeros(size(dataPS.allArea));
morph(dataPS.allArea<=thLow) = 2;
morph(dataPS.allArea>=thHigh) = 1;


%% Identify overlap indices

% initialize data structure
morphIdxs = struct();
morphIdxs.allCells.ste = cell(nFOV,nDays);
morphIdxs.allCells.pyr = cell(nFOV,nDays);
morphIdxs.commonCells.ste = cell(nFOV,nDays);
morphIdxs.commonCells.pyr = cell(nFOV,nDays);

dataIdxs = struct();
dataN = struct();
dataRBR = struct();

% loop through FOV
for ff = 1:nFOV
    % current FOV morphology classes
    diamCur = morph(dataPS.fovIdx==ff);

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
save([outPath 'morphIdxs.mat'],'morphIdxs')

