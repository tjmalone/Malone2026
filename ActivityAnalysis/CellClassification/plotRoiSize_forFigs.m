%% plotRoiSize_forFigs

clear; close all; clc

p1 = 'D:\AD_Project\imagingData\analysis_ManualSelection';
cd(p1)

% load data
load('roiSizeFinal.mat','allArea')
useData = allArea;
len = length(useData);

% define ks plot parameters
width = 4.5;
pts = 0:1:400;
legs = 'combined';

% define separation parameters
center = 158;
cRange = 20;


%% Plot size curves

figure; hold on

% plot ks density
[k,xi] = ksdensity(useData,pts,'width',width);
plot(xi,k*len)

% plot x lines
xline(center)
xline(center+cRange)
xline(center-cRange)

% set labels
legend(legs)
ylim([0 100])
ylabel('Cell count density')
xlabel('Cell area (\mum^2)')

sgtitle('two-photon FOVs')
savefig('2pPS_forFig.fig')
