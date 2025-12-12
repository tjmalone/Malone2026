%% plotOverlays
% Plots an overlay of PCAICA identified ROIs on an average projection
% image.

clear all; close all; clc;


%% Generate average projection overlays

p = pwd;

load('allROIs.mat')
roiN = size(roi,3);

fMC = 'Z:\labMembers\TM\novelAll\draw\ID_210209\old\ID20210209_0305Loc1_oldEnv.jpg';
mcTif = mat2gray(imread(fMC));

mkdir('overlays')
cd('overlays')

parfor k = 1:roiN
    bord = bwperim(roi(:,:,k));
    
    curTif = mcTif;
    curTif(bord) = 1;
    
    fig = figure;
    imshow(curTif)
    title(['ROI ' num2str(k)])
    savefig(fig,['ROI_' num2str(k)])
    
    close(fig)
end

cd(p)


%% Display figures

close all

p = pwd;

load('allROIs.mat')
cd('overlays')

roiN = size(roi,3);

for k = 1:roiN
    
    fig = openfig(['ROI_' num2str(k)]);
    
%     set(gcf, 'Position',  [0, 50, 550, 550])
    movegui(fig,[25 200]);

%     imcontrast(fig)
    cd(p)
    pause
    cd('overlays')
    
    close(fig)
end

cd(p)

