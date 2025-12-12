%% findLocalCellTypes
% Find cue, grid, and speed cells that meet the requirement for given
% sessions and saves a local cell type file. Set up for AD learning
% analysis.
%


%% Initialize

clear; clc; close all;

p1 = '/MATLAB Drive/FY2025/imagingData';
cd(p1)

% load alignment data
load('foldersLearning.mat','foldersLearning','trueDays')
nFOV = length(foldersLearning);
nDays = length(foldersLearning{1});


%% Identify local cell types

cellNums = struct();
cellNums.speedPos = nan(nFOV,nDays);
cellNums.speedNeg = nan(nFOV,nDays);
cellNums.total = nan(nFOV,nDays);

fNames = {'speedPos','speedNeg'};

% cycle through FOV
for ii=1:nFOV
    disp(num2str(ii))

    % initialize sessions folders
    dfolders = foldersLearning{ii};
    
    % cycle through sessions
    for jj=1:nDays
        % identify true day index
        trueDay = find(trueDays(ii,:)==jj);

        % skip days past day limits
        if isempty(trueDay)
            continue
        end

        % move to data directory
        newFolderName = replace(dfolders{trueDay},'D:\AD_Project\imagingData\data','/MATLAB Drive/FY2025/SubData2');
        newFolderName = replace(newFolderName,'\','/');
        cd(newFolderName)

        % load speed cell numbers
        load('localCellTypes_speed_sig','localTypeMat','localTypeIdx')
        for ff = 1:2
            cellNums.(fNames{ff})(ii,jj) = sum(localTypeMat(:,ff));
        end
        cellNums.total(ii,jj) = length(localTypeMat);

    end
end

cd(p1)


%%

% calculate percentages
cellNums.perPos = cellNums.speedPos./cellNums.total*100;
cellNums.perNeg = cellNums.speedNeg./cellNums.total*100;

cellNums.perPosMean = mean(cellNums.perPos,'all','omitnan');
cellNums.perNegMean = mean(cellNums.perNeg,'all','omitnan');

cellNums.perPosGlob = sum(cellNums.speedPos,'all','omitnan')/sum(cellNums.total,'all','omitnan')*100;
cellNums.perNegGlob = sum(cellNums.speedNeg,'all','omitnan')/sum(cellNums.total,'all','omitnan')*100;


cellNums.perPosTime = mean(cellNums.perPos,1,'omitnan');
cellNums.perNegTime = mean(cellNums.perNeg,1,'omitnan');

cellNums.relPosNeg = cellNums.speedPos./cellNums.speedNeg;


%% Split by genotype

% set plot fields
plotFields = {'perPos','perNeg'};
nFlds = length(plotFields);

% load groups
load('groupIDs.mat')
nSexes = length(sexes);
nGeno = length(groups);

% set plot parameters
colors6 = {[0 0 1],[1 0 0];[102 105 255]/255,[255 102 102]/255;[0 0 175]/255,[175 0 0]/255};
X = ["FE",string(1:10)];
dayCats = {1,2:11};
xLims = [0.5 11.5];
yLims = [0 20];

% set save folder name
svFile = [p1 '/Figures/Figures' char(datetime('today','format','yyyyMMdd')) '/cellNumbers'];
mkdir(svFile);

% initilaize figure
figure
tiledlayout(nFlds,nSexes)
sgtitle('Cell numbers')

for ff = 1:nFlds
    % set current data
    curData = cellNums.(plotFields{ff});

    for ss = 1:nSexes
        % set panel
        nexttile(ss + nSexes*(ff-1)); hold on
        title([plotFields{ff} ': ' sexIDs{ss}])

        % split data by type
        plotData = cell(nGeno,nDays);
        for gg = 1:nGeno
            plotData(gg,:) = num2cell(curData(intersect(sexes{ss},groups{gg}),:),1);
        end

        % plot data
        errorSig(X,plotData,colors6(ss,:),groupIDs,dayCats);

        % set labels
        xlim(xLims)
        ylim(yLims)
    end
end

% save figure
savefig([svFile '/speedCellNumbers'])

