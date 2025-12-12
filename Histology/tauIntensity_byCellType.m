%% tauIntensity_byCellType

clear; close all; clc

% move to directory
p1 = 'Z:\labMembers\KC\AD_Project\Histology\tau_RC\FinalVersion';
cd(p1)

% set results directory
p2 = 'Z:\labMembers\KC\AD_Project\histology_KC1.5mm\reelCalTau';

% find folders with data
searchFile = 'tau_1st.zip';
foldsNewDataFold = findSubF(searchFile,2,p2,0);

% remove searchFile
stringName = ['\' searchFile '\'];
folds = erase(foldsNewDataFold,stringName);
fNum = length(folds);


%%

% define cell type groups
cellTypes = {'rellin','calbindin'};
cellCodes = {'R','C'};
nTypes = length(cellTypes);
repCodes = {'1st','2nd'};
nReps = length(repCodes);
nBKG = 10;

% initilize data array
cellData = cell(fNum,1);
bkgAll = zeros(fNum,nReps);

for ff = 1:fNum
    curFolder = folds{ff};
    cellData{ff} = cell(nReps,nTypes);

    for rr = 1:nReps
        % load image
        curIm = tiffreadVolume([curFolder '\MEC_' repCodes{rr} '.tif']);
        useIm = double(curIm(:,:,4));
        [rows,cols] = size(curIm,1:2);

        % load background rois
        bkgFile = [curFolder '\trueBKG_' num2str(rr)];
        [bkgROIs,~] = extractImageJROI(rows,cols,[bkgFile '.zip']);
        if size(bkgROIs,3)~=nBKG
            error('Incorrect background ROI number')
        end

        % calculate background intensity
        bkgData = zeros(nBKG,1);
        for roi = 1:nBKG
            curVals = useIm(bkgROIs(:,:,roi)==1);
            bkgData(roi) = mean(curVals,'omitnan');
        end
        curBkg = mean(bkgData);
        bkgAll(ff,rr) = curBkg;
        if isnan(curBkg)
            error('Invalid background')
        end

        for tt = 1:nTypes
            % load current ROIs
            curFile = [curFolder '\true' cellCodes{tt} 'T_' num2str(rr)];
            if isfile([curFile '.zip'])
                [allROIs,~] = extractImageJROI(rows,cols,[curFile '.zip']);
            elseif isfile([curFile '.roi'])
                [allROIs,~] = extractImageJROI(rows,cols,[curFile '.roi']);
            else
                cellData{ff}{rr,tt} = NaN;
                continue
            end

            % extract intensity from ROIs
            nROIs = size(allROIs,3);
            cellData{ff}{rr,tt} = zeros(nROIs,1);
            for roi = 1:nROIs
                curVals = useIm(allROIs(:,:,roi)==1);
                cellData{ff}{rr,tt}(roi) = mean(curVals,'omitnan')-curBkg;
            end
        end
    end
end

save('fullCellData.mat','cellData','bkgAll')

%% Validate cell numbers

clear; close all; clc

load('fullCellData.mat','cellData')

% get calculated cell numbers
cellNums = cellfun(@(x) cellfun(@(y) sum(~isnan(y)),x), cellData,'UniformOutput',false);
cellNumsAll = cat(1,cellNums{:});

% get true cell numbers
load('Z:\labMembers\KC\AD_Project\Histology\tau_RC\FinalVersion\calculateFinalCellCount.mat','cellCountFinal')
cmpNumsAll = cellCountFinal(:,4:5);

if ~isequal(cellNumsAll,cmpNumsAll)
    error('Mismatched cell numbers')
end


%% Calculate intensites per FOV

% load sex info
load('Z:\labMembers\KC\AD_Project\Histology\tau_RC\FinalVersion\ByFOV\plotRCTBySex\sexIDs.mat',...
    'sexID','sexNum')
useSex = repelem(sexNum,2,1);
nSex = length(sexID);

% calculate final intensity
cellInts = cellfun(@(x) cellfun(@(y) mean(y,'omitnan'),x), cellData,'UniformOutput',false);
cellIntsAll = cat(1,cellInts{:});
cellTypes = {'stellate','pyramidal'};
nTypes = size(cellIntsAll,2);

% split data by sex and cell type
finalInt = cell(nTypes,nSex);
for ss = 1:2
    for cc = 1:2
        finalInt{cc,ss} = cellIntsAll(useSex==ss,cc);
    end
end

% plot data
colors2 = {[0.6 0 0],[1 0.3 0.3]};
pShow = generatePShow(2,[1 2]);
figure; hold on
[h,p] = barGroup(cellTypes,finalInt,'violin',colors2,pShow);

% set labels
legend(h,sexID)
ylabel('Arbitrary Units')
title('Tau Intensity in STE/PYR cells (per FOV)')

meanInt = cellfun(@(x) mean(x,'omitnan'),finalInt);


%% Calculate intensites per cell

% calculate cell intensity
cellIntsCell = cat(1,cellData{:});

% split data by sex and cell type
finalIntCell = cell(nTypes,nSex);
for ss = 1:2
    for cc = 1:2
        finalIntCell{cc,ss} = cat(1,cellIntsCell{useSex==ss,cc});
    end
end

% plot data
pShow = generatePShow(2,[1 2]);
figure; hold on
[h,p] = barGroup(cellTypes,finalIntCell,'violin',colors2,pShow);

% set labels
legend(h,sexID)
ylabel('pTau intensity (A.U.)')
title('Tau Intensity in STE/PYR cells (per FOV)')

% perform full statistics
nUnits = 'cells';
testName = 'two-tailed paired Students t-test';
testPair = 0;
testMC = 0;
testLimitP = 0;
statData = finalIntCell';
testCat = 'Tau intensity by cell';
outStats = [testCat ttestEffectSize(statData(1,:),statData(2,:),testName,nUnits,testPair,testMC,testLimitP)];



%% Plot heat map

ttls = 'pTau intensity';

% calculate group means
meanIntCell = cellfun(@(x) mean(x,'omitnan'),finalIntCell);

curData = num2cell(fliplr(meanIntCell))';

morphData = cell(2,1);
for kk = 1:2
    morphData{kk} = cat(1,curData{:,kk});
end

% define data
curMap = cell(3,3,2);
for kk = 1:numel(curMap)
    curMap{kk} = nan;
end
curMap(1,2:3,:) = cat(2,{0;0},morphData);
curMap(2:3,2:3,:) = cat(3,{0,0;0,0},curData);

% make heat map
h = heatmapGridSplit(curMap,{'Sex','All','Female','Male'},...
    {'Morphology','All','ste','pyr'},-1,{2:3},{2:3},{'eastoutside'});
sgtitle(ttls,'FontSize',18)

