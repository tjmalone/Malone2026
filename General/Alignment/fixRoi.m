function bROIs = fixRoi(fold,roiFile,actFile)
%% fixROIs.m
% Removes ROIs with no significant transients. For use within pcaica
% folders.
%
% Inputs: (all optional)
%       fold = folder in which to process ROIs. Default is current
%       directory.
%
%       roiFile = path/file name to matrix containing roi data. Variable
%       within file must be called "roi" with rois in 3rd dim. Default is
%       "allROIs.mat"
%
%       actFile = path/file name to matrix containing activity measurement
%       for all cells. Variable within file must be called dfof_sig with
%       rois in 3rd dim. Default is "dfof_sig_#_cell.mat" with variable
%       cell number
% 

%% Process inputs

% determine path
if nargin==0
    fold = pwd;
end

% change path
p = pwd;
cd(fold)

% determine roi file
if nargin<2 || isempty(roiFile)
    roiFile = 'allROIs.mat';
end

% load rois
load(roiFile,'roi')

% determine activity file
if nargin<=3 || isempty(actFile)
    nCells = size(roi,3);
    actFile = ['dfof_sig_' num2str(nCells) '_cells.mat'];
end

% load activity
load(actFile,'dfof_sig')


%% Identify and process bad rois

bROIs = find(all(dfof_sig==0,1));

i = setdiff(1:size(roi,3),bROIs);
roi=roi(:,:,i);


%% Save procesed rois

save('allROIsNew.mat','roi');
save('bROIs.mat','bROIs')

cd(p)

end