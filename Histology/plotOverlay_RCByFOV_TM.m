%% plotOverlay_RCByFOV
% graph of Reel/cal/tau overlap
% num = numerator
% den = denominator
%tauNum = tau/excit
%tauDen = excit/tau

clear; close all; clc

baseFold1 = 'Z:\labMembers\KC\AD_Project\Histology\tau_RC\FinalVersion\';
load([baseFold1 'calculateFinalCellCount.mat']);

cd([baseFold1 'ByFOV'])
saveFold2 = pwd;

% define labels
titleString = {'Percent of tau+ cells that are stellate or pyramidal',...
    'Percent of stellate or pyramidal cells that are tau+'};
yString = {'Percent of tau cells (%)','Percent of stel or pyr cells (%)'};
yLims = [0 90; 0 60];
xCellType = {'Ste', 'Pyr'};
pShow = [1 2];


%% Categorize FOVs based on Mouse

load('Z:\labMembers\KC\AD_Project\Histology\tau_RC\ADMouseOrder.mat');
numFOV = length(cellCountFinal);

dataUse = {perOfTau,perOfExcit};
colorsCells = {[0 1 0],[1 0 1]};

% initiaize figue and output data
mousePercentFinal = cell(1,2);
pVal = cell(1,2);
figure; hold on

% initialize stats parameters and array
outStats = {};
nUnits = 'FOV';
testName = 'two-tailed paired Students t-test';
testPair = 1;
testMC = 0;
testLimitP = 0;

% loop through calculation types
for ii = 1:2
    % calculate percetage
    mousePercentFinal{ii} = dataUse{ii}*100;

    % plot bar graph
    subplot(1,2,ii);
    [~,pVal{ii}] = barGroup(xCellType,mousePercentFinal{ii},'violin',colorsCells,pShow,'pair');
    ylim(yLims(ii,:))
    ylabel (yString{ii})
    title(titleString{ii})

    % perform full statistics
    statData = mousePercentFinal{ii};
    testCat = titleString{ii};
    outStats(end+1,:) = [testCat ttestEffectSize(statData(:,1),statData(:,2),testName,nUnits,testPair,testMC,testLimitP)];
end

% save figure and data
savefig('plotOverlay_RCByFOV.fig')
save('tau_RC_overlapByMouse.mat','mousePercentFinal','pVal');

