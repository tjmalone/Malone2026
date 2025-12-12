%% Prep alignment

clear; close all; clc

st = [2,12];
load('folders.mat')
cats = [];
fovCells = cell(1,2);

saveFolder = fullfile(pwd,num2str(st(1)),[num2str(st(1)) '_' num2str(st(2))]);

for ii = 1:2
    curF = folders{st(ii)};
    
%     if isfile([curF '\allROIsNew.mat'])
%         load([curF '\allROIsNew.mat']);
%     else
        load([curF '\allROIs.mat']);
%     end
    
    if ~isempty(cats)
        load([curF 'byCat\catIDs.mat']);
        
        IDs = cat(2,catIDs{cats});
        roi = roi(:,:,IDs);
    end

    roiSum = 0.5*(sum(roi,3)>0);
    imwrite(roiSum,[saveFolder '\all_rois' num2str(ii) '.tif'])
    
    fovCells{ii} = roi;
end

cd(saveFolder)


%% Set Inputs

goAhead = 1;    % save manual alignment
close all
degRotation = 0;
xdiff = 257;
ydiff = 247;
resizeFactor_W = 506;
resizeFactor_H = 502;

% Process Inputs
% DON'T TOUCH
szScale = size(fovCells{2},1:2);
resizeFactor_W = resizeFactor_W / szScale(2);
resizeFactor_H = resizeFactor_H / szScale(1);
xdiff = xdiff - szScale(2)/2*resizeFactor_W;
ydiff = ydiff - szScale(1)/2*resizeFactor_H;


% Align FOV
if ~goAhead
    FOV_manual_alignment_new(fovCells,xdiff,ydiff,degRotation,goAhead,resizeFactor_W,resizeFactor_H);
else
    [spatial_footprints_corrected,centroid_locations_corrected,footprints_projections_corrected,centroid_projections_corrected,maximal_cross_correlation,alignment_translations,overlapping_FOV] = ...
        FOV_manual_alignment_new(fovCells,xdiff,ydiff,degRotation,goAhead,resizeFactor_W,resizeFactor_H);
    
    save([saveFolder '\manualAlignment.mat'],'spatial_footprints_corrected','centroid_locations_corrected','footprints_projections_corrected','centroid_projections_corrected','maximal_cross_correlation','alignment_translations','overlapping_FOV','-v7.3')
    savefig([saveFolder '\manualAligned.fig']);
%     close
end


%% Run alignment

cd ..\..
alignmentSingle_suite2p(folders,st(1),st(2))

