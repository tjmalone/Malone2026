%% overlay.m

clear all; close all; clc


%% Load pcaica

f='TSeries*';

d=dir(f);
fPath=d.name;

load('pcaica/allROIs.mat','roi')

sz = size(roi,1);
dataPr = zeros(sz);

load([fPath '/suite2p/plane0/Fall.mat'],'ops');

mxProj = ops.max_proj;
xRng = ops.xrange;
yRng = ops.yrange;

dataPr(yRng(1)+1:yRng(2),xRng(1)+1:xRng(2)) = mxProj;

load('motionCorrected5times/matrix/matrix_001.mat','data');

mvPcaica = data;


%% Load suite2p

sz = 512;

load([fPath '/suite2p/plane0/Fall.mat'],'stat','iscell');

roiS2P = zeros(sz);

nRoi = length(stat);

for i = 1:nRoi
    if ~iscell(i,1)
        %         continue
    end
    
    curRoi = [stat{i}.ypix;stat{i}.xpix]';
    nPix = size(curRoi,1);
    
    for j = 1:nPix
        roiS2P(curRoi(j,1),curRoi(j,2)) = 1;
    end
end


%% Create S2P movie

% load tif stack
fileID = fopen([fPath '/suite2p/plane0/data.bin'],'r');


mvS2P = fread(fileID,512*512*1000,'*int16');

mvS2P = reshape(mvS2P, 512, 512, []);

fclose(fileID);

mvNormS2P = permute(mat2gray(mvS2P),[2 1 3]);

clear mvS2P


%% Generate overlay

% combine rois
roiPcaica = sum(roi,3);

% create outline
bordPcaica = bwperim(roiPcaica);
bordS2P = bwperim(roiS2P);

% scale original data
normPrPcaica = mat2gray(dataPr);
normPrS2P = mat2gray(dataPr);

% add projection outlines
normPrPcaica(bordPcaica==1) = 1;
normPrS2P(bordS2P==1) = 1;


%% Generate S2P Movie

% scale original data
mvNormS2P(1:yRng(1),:,:) = 0;
mvNormS2P(yRng(2)+1:end,:,:) = 0;
mvNormS2P(:,1:xRng(1),:) = 0;
mvNormS2P(:,xRng(2)+1:end,:) = 0;

% add movie outlines
[row,col] = find(bordS2P==1);
ind = sub2ind(size(bordS2P),row,col);

for i = 1:size(mvNormS2P,3)
    tmp = mvNormS2P(:,:,i);
    tmp(ind) = 1;
    mvNormS2P(:,:,i) = tmp;
end


%% Generate PCAICA Movie

% scale original data
mvNormPcaica = mat2gray(mvPcaica);

mvNormPcaica(1:yRng(1),:,:) = 0;
mvNormPcaica(yRng(2)+1:end,:,:) = 0;
mvNormPcaica(:,1:xRng(1),:) = 0;
mvNormPcaica(:,xRng(2)+1:end,:) = 0;

% add movie outlines
[row,col] = find(bordPcaica==1);
ind = sub2ind(size(bordPcaica),row,col);

for i = 1:size(mvNormPcaica,3)
    tmp = mvNormPcaica(:,:,i);
    tmp(ind) = 1;
    mvNormPcaica(:,:,i) = tmp;
end


%% Display max projection

figure
imshow(normPrPcaica);

figure
imshow(normPrS2P);


%% Play movie

implay(mvNormPcaica,100);

implay(mvNormS2P,100);

%% Clear excess variables

% clear d data f mxProj ops xRng yRng
%



