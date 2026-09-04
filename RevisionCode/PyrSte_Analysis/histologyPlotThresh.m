%% Set inputs

clear; close all; clc
p1 = 'Z:\labMembers\TM\AD_Project\Revision\Histology\revisionHistology';


%% Process histology FOV

% set base folder
pHist = 'Z:\labMembers\TM\AD_Project\Revision\Histology\revisionHistology';
cd(pHist)

% find all validation selections
foldsHV = findSubF('HV',2,[],0);

% define use FOV
% useFOV = [3:6 10:12 16:17];
% useFOV = [3:6 10:12 16];
useFOV = [3:6 11:12 16];

nFOV = length(useFOV);

% get all selection thresholds
threshHist = zeros(nFOV,1);
for ff = 1:nFOV
    curF = useFOV(ff);
    load([foldsHV{curF} '\indThresh.mat'],'splitThr')
    threshHist(ff) = splitThr;
end


%% Process imaging FOV

% set base folder
pImg = 'Z:\labMembers\TM\AD_Project\imagingData\data_FOV\refImages\highConfidenceFOVs\Final';
cd(pImg)

% define use FOV
useNums = [9,13,14,15,17,19,20,21,22,23,24,27,29,30,31,32,35,37,38,39,40];
nFOV = length(useNums);

% get all selection thresholds
threshImg = zeros(nFOV,1);
for ff = 1:nFOV
    numCur2 = num2str(useNums(ff),'%02d');
    load(['params' numCur2 '.mat'],'splitThr')
    threshImg(ff) = splitThr;
end
cd(p1)


%% Plot thresholds

% define inputs
x = {'Histology','Imaging'};
threshAll = {threshHist,threshImg};
pShow = [1 2];
plotInd = 2;
offsets = [-10 0 10];
nOffsets = length(offsets);

% set t-test settings
nUnits = ' FOVs';
testName = 'two-tailed unpaired Students t-test';
testPair = 0;
testMC = 0;
testLimitP = 0;

figure; tiledlayout(1,nOffsets)

outStats = {};
for ii = 1:nOffsets
    nexttile(ii); hold on

    % apply threshold offset
    curThresh = cellfun(@(x) x+offsets(ii),threshAll,'UniformOutput',false);

    barGroup(x,curThresh,'violin',[],pShow,[],[],plotInd)

    % set labels
    ylim([120 220])
    ylabel('Areas (\mum^2)')


    outStats(end+1,:) = ['' ttestEffectSize(curThresh(:,1),curThresh(:,2),testName,nUnits,testPair,testMC,testLimitP)];
end
