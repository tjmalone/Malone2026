%% generateRefAuto
% generates a reference image for given FOV using the activity data of
% active ROIs

clear; close all; clc

% move to base folder
p1 = 'D:\AD_Project\imagingData\';
cd(p1)

% define save folder
svFolder = [p1 'data_FOVauto\'];

% define stack source folder
stackFolders = {'E:\tempStacks\','G:\Mice\MC\','H:\ADBatch6_MotionCorrect\'};
nSource = length(stackFolders);

% load activity folder
load('foldersLearning.mat','foldersLearning')

useNums = 31%1:42;

%%

for useNum = useNums
    %% Load all data

    % if isfile(['data_FOVauto\refAuto_' num2str(useNum) '.tif'])
    %     continue
    % end

    curActFolder = foldersLearning{useNum}{1};

    % load rois
    load([curActFolder 'allROIs.mat'],'roi')
    nX = size(roi,2);
    nY = size(roi,1);
    nRoi = size(roi,3);

    % load activty
    d = dir([curActFolder 'dfof_sig*.mat']);
    load([curActFolder d(1).name],'dfof_sig')

return
    %% Process stacks

    % get reference folder
    curFolderSplit = strsplit(curActFolder,'\');
    curRefFolder = fullfile(curFolderSplit{5:8},'motionCorrected5times');

    % get reference files
    for ff = 1:nSource
        curRefFiles = dir([stackFolders{ff} curRefFolder '\*.tif']);
        if ~isempty(curRefFiles)
            break
        end
    end
    nMC = length(curRefFiles);

    % check for source file success
    if nMC==0
        warning(['Motion corrected source not found: Roi ' num2str(useNum)])
        continue
    end

    roiBuff = 15;

    % initialize blank image
    refImAll = nan(size(roi));
    actNAll = nan(nRoi,1);

    curIdx = 1;

    for mc = 1:nMC

        % load reference stack
        curRef = loadtiff([curRefFiles(mc).folder '\' curRefFiles(mc).name]);

        curSpan = curIdx:curIdx+size(curRef,3)-1;
        curIdx = curIdx+size(curRef,3);


        %% Loop across rois

        for rr = 1:nRoi
            % get roi bounding box
            roiProp = regionprops(roi(:,:,rr),'BoundingBox');
            bbox = roiProp.BoundingBox;

            % define x range (2nd dimension)
            x1 = max(1,floor(bbox(1)-roiBuff));
            x2 = min(nX,ceil(bbox(1)+bbox(3)+roiBuff));

            % define y range (1st dimension)
            y1 = max(1,floor(bbox(2)-roiBuff));
            y2 = min(nY,ceil(bbox(2)+bbox(4)+roiBuff));

            % define activity range
            idxAct = find(dfof_sig(curSpan,rr)>0);
            nAct = length(idxAct);
            actNAll(rr) = sum([actNAll(rr),nAct],'omitnan');

            % skip cells with no activity
            if nAct>0
                % get active roi range
                roiRef = curRef(y1:y2,x1:x2,idxAct);
                roiSum = sum(roiRef,3,'omitnan');

                % store all subset
                roiCatCur = cat(3,refImAll(y1:y2,x1:x2,rr),roiSum);
                refImAll(y1:y2,x1:x2,rr) = sum(roiCatCur,3,'omitnan');
            end
        end
    end

    refImInd = refImAll./reshape(actNAll,1,1,[]);
    refImFinal = mean(refImInd,3,'omitnan');
    refImAllFinal = sum(refImAll,3,'omitnan')/sum(actNAll,'omitnan');

    save(['data_FOVauto\refAuto_' num2str(useNum)],'refImInd','refImFinal','refImAllFinal')


    %% Sort roi by position and save as tif

    % get first non-zero pixel per roi
    idxBound = zeros(nRoi,1);
    for rr = 1:nRoi
        firstIdx = find(refImInd(:,:,rr)>0,1);
        if ~isempty(firstIdx)
            idxBound(rr) = firstIdx;
        else
            idxBound(rr) = 1;
        end
    end

    % sort pixels and rois
    [~,idxSort] = sort(idxBound);
    refImIndSort = refImInd(:,:,idxSort);

    saveastiff(refImIndSort,['data_FOVauto\refAuto_' num2str(useNum) '.tif']);

end

return


%%

figure
imshow(refImFinal/max(refImFinal,[],'all','omitnan'))


%%

% imagsc(refImInd(:,:,rr)/max(refImInd(:,:,rr),[],'all','omitnan'))

figure

for rr=1:nRoi
    imshow(refImInd(:,:,rr)/max(refImInd(:,:,rr),[],'all','omitnan'))
    axis('off','square')

    pause
end




%%

close all

figure; imagesc(sum(roi,3)); axis('square','off')
figure; imagesc(mean(curRef,3)); axis('square','off')



