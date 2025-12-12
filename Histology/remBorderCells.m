%remBorderCells
% Output: indeces of rois that are a certain percentage within the MEC layer 2 roi

clear; close all; clc

% set values
versionName = 'FinalVersion';
foldNeuroden = 'Z:\SNM\labMembers\KC\AD_Project\Histology\reelCalNeurodegen';
borderCutoff = 90;
plotTestA = 0;

dataFolder = 'Z:\SNM\labMembers\KC\AD_Project\histology_KC1.5mm\reelCalTau';
types = {'reelin','calbindin'};
images = 'MEC_brightCal.tif';
cellFiles = {[types{1} '_neurodegen_brightCal.zip'], [types{2} '_neurodegen_brightCal.zip']};
mecL2 = 'degen_brightCal_v4.roi';

useT = 1;
useIdx = 3;
width = 0.001;

cd(dataFolder)
% Find folders
foldsNewDataFold = findSubF(images,2,[],0)';

% remove searchFile
stringName = ['\' images '\'];
folds = erase(foldsNewDataFold,stringName);
fNum = length(folds);

cellsInsideBorder = {};

% make fold for plots
cd(foldNeuroden);
mkdir(versionName);
cd(versionName);
mkdir('remBorderCells')
cd('remBorderCells')
foldSave = pwd;

    %% Loop through all mice

% for f = 1
for f = 1:fNum
     % Load Image
    
    cd(folds{f})

    curIm = tiffreadVolume(images);
    
    % normalize image
    curIm = double(curIm);
   
    curIm = curIm(:,:,useIdx(useT));
    
    % Extract nissl ROIs
    [rows,columns] = size(curIm,1:2);

    for cc = 1:length(types)
        testFile = [types{cc} '_neurodegen_brightCal_v3.zip'];
        if exist(testFile,'file') == 2
            cells{f,cc} = testFile;
        else
            cells{f,cc} = cellFiles{cc};
        end

        [allROIs] = extractImageJROI(rows,columns,cells{f,cc});
        nRoi = size(allROIs,3);

        % Extract MEC roi and find indices
        MECreg = extractImageJROI(rows,columns,mecL2);  
        MECreg_indices = find(MECreg);

        %check if there are pixels outside of MEC roi
        count = 1;
        count2 = 1;
        for ii = 1:nRoi
            % Find the indices of the nissl roi
            roi_indices = find(allROIs(:,:,ii));

            % check how many pixels overlap
            isPresent = ismember(roi_indices, MECreg_indices);

            % Calculate percentage of matching elements
            matchPercentage = sum(isPresent) / numel(roi_indices) * 100;

            %check if all indices are in MEC roi
            if matchPercentage >= borderCutoff
                cellsInsideBorder{f,cc}(count,1) = ii;
                count = count + 1;
            else
                cellsOutsideBorder{f,cc}(count2,1) = ii;
                count2 = count2 + 1;
            end
        end

        if plotTestA
            %Plot with imagesc
            figure;
            roiCmb1 = sum(allROIs(:,:,1:nRoi),3) + MECreg;
            subplot(1,2,1);
            imagesc(roiCmb1)
            axis('equal','off')

            roiCmb2 = sum(allROIs(:,:,cellsInsideBorder{f,cc}),3) + MECreg;
            subplot(1,2,2);
            imagesc(roiCmb2)
            axis('equal','off')

            ax1 = subplot(1,2,1);
            ax2 = subplot(1,2,2);
            linkaxes([ax1, ax2], 'xy');  % Link both x and y axes

            nameFig = [mice{f}(:,:) '-plane' num2str(pp) '-removeBorderNissl-' num2str(borderCutoff) '%'];
            sgtitle(nameFig);
            cd(foldSave);
            savefig([nameFig '.fig']);
        end

    end
end
%%
cd(foldSave);

% folder where code is saved
baseFolder = 'D:\MATLAB\Code\AnalysisCode\Analysis\HistologyCode\ReelCalNeurodegen\withRemoveBorderCells';

% copy self to current directory
copyfile([baseFolder '\remBorderCells.m'],'remBorderCells.m');

save('removeBorderCells.mat','cellsInsideBorder','cellsOutsideBorder','folds','cells');
