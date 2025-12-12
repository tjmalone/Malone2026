%% calculateFinalCellCount_Manual
% for 1.5 mm
% for reelin/calbindin/tau analysis

clear; close all; clc

dataFolder = 'Z:\labMembers\KC\AD_Project\histology_KC1.5mm\reelCalTau';
foldSave = 'Z:\labMembers\KC\AD_Project\Histology\tau_RC';
tauFiles ='tau_*';
numMice = 14;
numSlice = numMice*2;
miceInfoLoc = 'Z:\labMembers\KC\AD_Project\Histology\tau_RC';

% Individual cell type specific settings
key = {'reelin_','calbindin_','tau_', 'overlay_reelin-tau-','overlay_calbindin-tau-', 'full_overlay_'};
ending = {'1st','2nd'};

images = 'MEC_*';

% %identify data folders
cd(dataFolder)

% Find folders with Data
searchFile = 'MEC_1st.tif';
foldsNewDataFold = findSubF(searchFile,2,[],0)';

% remove searchFile
stringName = ['\' searchFile '\'];
folds = erase(foldsNewDataFold,stringName);
fNum = length(folds);

rois = {};

cd(miceInfoLoc)
load('ADMouseOrder.mat')

cellCountInitial = [];

%% Loop through all FOVS
% need to adjust
roiXIndex = cell(numSlice,1); % indices of cfos rois that overlap with X
sizeXIndex = zeros(numSlice,1);
roiCmb = {};

count = 1;
nearDistance = {};

% loop through all mice
for f = 1:fNum
% for f = 1 %%%testing
%%%
    cd(folds{f})

    tFiles = dir(tauFiles);
    if size(tFiles,1) >= 1
        
        MECImages = dir(images);
    
        % loop through all planes for each mouse
        for g = 1:size(tFiles,1)    
            cd(folds{f})
            filePath{count,1} = pwd;
            tiffName{count,1} = MECImages(g).name;
                             
            % loop through each of the roi 6 files
            for ii = 1:length(key)
                roiFilename = [key{ii} ending{g}];

                if isfile([roiFilename '.zip'])
                    [cvsROIs] = ReadImageJROI([roiFilename '.zip']);
                    nCur = length(cvsROIs);
                elseif isfile([roiFilename '.roi'])
                    nCur = 1;
                else
                    nCur = 0;
                end
            
                cellCountInitial(count,ii) = nCur;
                
            end
            count = count + 1;
        end
    end
end % needed for the for loop


%% Adjust for full overlay cells
cellCountFinal(:,1) = cellCountInitial(:,1)-cellCountInitial(:,6); %reel
cellCountFinal(:,2) = cellCountInitial(:,2)-cellCountInitial(:,6); %cal
cellCountFinal(:,3) = cellCountInitial(:,3); %tau
cellCountFinal(:,4) = cellCountInitial(:,4)-cellCountInitial(:,6); %reel tau
cellCountFinal(:,5) = cellCountInitial(:,5)-cellCountInitial(:,6); %cal tau
cellCountFinal(:,6) = cellCountInitial(:,6); %reel cal tau

% Calculate percentages

for ii = 1:2
    perOfTau(:,ii) = cellCountFinal(:,ii+3)./cellCountFinal(:,3);
    perOfExcit(:,ii) = cellCountFinal(:,ii+3)./cellCountFinal(:,ii);
end


%% Save data
cd(foldSave);
mkdir('FinalVersion')
cd('FinalVersion')
save('calculateFinalCellCount.mat','filePath','cellCountInitial','cellCountFinal','perOfTau','perOfExcit','key');


%% Check accuracy by comparing to auto
% cellCountManual = cellCountFinal;
% 
% foldAuto = 'Z:\SNM\labMembers\KC\AD_Project\Histology\tau_RC\B3456\1.5mm\Threshold_5um';
% cd(foldAuto)
% load('calculateFinalCellCount.mat','cellCountFinal')
% cellCountAuto = cellCountFinal; 
% 
% diffManVsAuto = cellCountAuto - cellCountManual;
% diffTotal = sum(diffManVsAuto);
% cd(foldSave);
% save('accuracyCheck.mat','diffManVsAuto','diffTotal')