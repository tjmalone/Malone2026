%% plotCuePairs
% plot sample selections of cue and anti-cue zones

clear; close all; clc

p1 = 'D:\AD_Project\imagingData\analysis_Decoding';
cd(p1)

% load cue info
load('cuePairInfo_NE.mat','cuePairs')
load('antiCuePairInfo_NE.mat','antiCuePairs')
nPlot = size(cuePairs,1);
colors = {[0 1 0],[1 0 1]};

% load cue templates
cueFolder = 'D:\AnalysisCode\PostAnalysis\Cues\6mEnv2\';

% load cue templates
load([cueFolder 'tempL.mat'])
load([cueFolder 'tempR.mat'])
colorsCue = [0.5 0.5 0.5];
cueTemp = [tempL,tempR];
cueX = 1:length(cueTemp);


%% Plot cue and anti-cue zones

figure; hold on
tiledlayout(1,2)

for ii = 1:2
    nexttile(); hold on

    if ii==1
        curPairs = cuePairs;
        nCur = nPlot;
    else
        curPairs  = antiCuePairs;
        nCur = 2*nPlot;
    end

    for jj = 1:nCur
        for kk = 1:2
            plot(curPairs{jj,kk},jj*ones(size(curPairs{jj,kk})),'Color',colors{ii});
        end
    end

    % plot cues
    for cc = 1:size(cueTemp,2)
        plotCues(cueX,cueTemp(:,cc),nCur,colorsCue);
    end

    xlim([0 120])
    ylim([0.5 nCur+0.5])
    set(gca,'XColor','none','YColor','none')
end

savefig('cuePairsExample.fig')

