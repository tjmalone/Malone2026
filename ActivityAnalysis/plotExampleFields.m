%%

clear; close all; clc

p1 = 'D:\AD_Project\imagingData';
cd(p1)

load('groupIDs.mat')
% disp(intersect(groups{1},sexes{2}));

load('data\dataByCell.mat','dataAll')
spatialSelectivity = dataAll.spatialSelectivity;

load('foldersLearning.mat','foldersLearning','alignsLearning')

% AD female examples (96 is primary, 180 is alternative)
% useMouse = 19;
% useDay = 11;

% WT female examples (19 is primary)
useMouse = 30;
useDay = 11;

% move to correct folder
cd(foldersLearning{useMouse}{useDay})

% get spatial selectivity for common cells only
useSS = spatialSelectivity{useMouse,useDay};
useAlign = alignsLearning{useMouse}(:,useDay);
useSkip = setdiff(1:length(useSS),useAlign);
useSS(useSkip) = nan;

open('D:\AD_Project\imagingData\data\TM230304-2\230512\loc1\TSeries_920_1380um-926_1\suite2p\gridAnalysis_sig\allCorrected\cell_96.fig')
subplot(4,1,1)
colorbar
colormap('parula')


open('D:\AD_Project\imagingData\data\TM240316-I\240514\loc1\TSeries-2373_2\suite2p\gridAnalysis_sig\allCorrected\cell_19.fig')


useMouse = [19 30];
useDay = [11 11];
useCell = [96 19];

figure; hold on
tiledlayout(1,2)
X = 2.5:5:597.5;

for ii = 1:2
    % load dfof and fields
    cd(foldersLearning{useMouse(ii)}{useDay(ii)})
    load('gridAnalysis_sig\allCellsCorrected.mat')

    % define dfof and fields
    curDfof = allCellsCorrected.dfofaveragesmooth(:,useCell(ii))';
    curFields = allCellsCorrected.inFieldBins{useCell(ii)};

    % get spatial selectivity
    curSS = spatialSelectivity{useMouse(ii),useDay(ii)}(useCell(ii));

    nexttile(ii)
    plot(X,curDfof)
    xlim([0 600])

    title(num2str(curSS,2))

end






