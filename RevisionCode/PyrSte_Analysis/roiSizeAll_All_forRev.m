%%

clear; close all; clc

p1 = 'D:\AD_Project\Revision\PyrSteAnalysis\';
cd(p1)

% image raw image and data path
imPath = 'Z:\labMembers\TM\AD_Project\imagingData\data_FOV\refImages\highConfidenceFOVs\Final\';
dataPath = imPath;

useNums = [9,13,14,15,17,19,20,21,22,23,24,27,29,30,31,32,35,37,38,39,40];
nFile = length(useNums);

pixelSize = 750/512; %for 750um FOV under  512x512 pixels
areaWidth = [6,8,10,12,14,16];
buff = 10;

allLong = [];
allShort = [];
allArea = [];


%% Collect global area

globArea = [];

for ii = 1:nFile
    numCur = num2str(useNums(ii));
    numCur2 = num2str(useNums(ii),'%02d');


    % load parameters
    load([dataPath 'params' numCur2 '.mat'],'params')
    globArea = [globArea; params.areaum];
end


%% Plot global area

ksSize = 7;

figure; hold on

% plot current ksdensity
[p,x] = ksdensity(globArea,0:1:400,'width',ksSize);

plot(x,p);
xlim([min(x) max(x)]);
title(['Area: width = ' num2str(ksSize)]);

xline([152 162 172])
ylim([0 0.01])
