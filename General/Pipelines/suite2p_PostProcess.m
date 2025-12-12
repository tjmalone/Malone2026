function suite2p_PostProcess(filePath,trackEnd,sizeThresh)
%% Use motion corrected data for Suit2p
% run suite2p directly using motion corrected tif stacks. Suite2p will
% output Fall.mat
%

if nargin==0 || isempty(filePath)
    filePath = pwd;
end

if nargin<2 || isempty(trackEnd)
    trackEnd = 600;
end

if nargin<3 || isempty(sizeThresh)
    sizeThresh = 25;
end

cd(filePath)


%% Load data

binWidth = 5;
trackStart = 0;

d = dir('*\*\Fall.mat');
fold = [d.folder '\' d.name];
if ~isfile('Fall.mat')
    movefile(fold,filePath)
end

load('Fall.mat','stat','ops','F','iscell');
load('abf.mat','abf')
load('abfFake.mat','abfFake')

[sz1,sz2] = size(ops.meanImg);


%% Extract ROIs

% identify true cells
npix = zeros(length(stat),1);
for ii = 1:length(stat)
    npix(ii) = stat{ii}.npix;
end

cellIdxA = find(iscell(:,1)==1);
cellIdxB = find(npix>sizeThresh);
cellIdx = intersect(cellIdxA,cellIdxB);

% plot cells
cellStat=stat(cellIdx);
roi=zeros(sz1,sz2,length(cellStat));
for n=1:length(cellStat)
    roiCur = zeros(sz1,sz2);
    
    for m=1:length(cellStat{n}.xpix)
        % extract roi pixels (convert from python 0 indexing)
        roiCur(cellStat{n}.ypix(m)+1,cellStat{n}.xpix(m)+1) = 1;
    end
    
%     % process rois
%     se = strel('diamond',2);
%     roiCur = imdilate(roiCur,se);
%     roiCur = imfill(roiCur,'holes');
%     roiCur = imerode(roiCur,se);
    
    roi(:,:,n) = roiCur;
end

save('allROIs.mat','roi');
plotROIShape(roi,1:size(roi,3));
axis('off','equal')
title('all rois');

saveas(gcf,'all_rois.tif');
saveas(gcf,'all_rois.fig');
close

F=F(cellIdx,:)';
filenameF=sprintf('%s_%d_%s.mat', 'F', size(roi,3), 'cells');
disp('saving F');
save(filenameF,'F','-v7.3');


%% Calculate dfof

disp('calculate dfof');

%in case the imaging finished after volatge recording
if size(F,1)>length(abf.imageIndex)
    F = F(1:length(abf.imageIndex),:);
    save(filenameF,'F','-v7.3');
end

% calculate dfof
dfof=zeros(size(F));

for cellindex=1:size(F,2)
    f=F(:,cellindex);
    f(f<0)=0;
    dfof(:,cellindex)=dfof_transients_interactive_batch_automatic(f,abfFake.t);
    disp(cellindex);
end

dfof(isnan(dfof))=0;
dfof(isinf(dfof))=0;

filenamedfof = sprintf('%s_%d_%s.mat', 'dfof', size(F,2), 'cells');
save(filenamedfof,'dfof');

%plot activities
%speedThreshold
disp('calculate speed threshold')

[speedThreshold]= speedThreshold1D( 1,length(abfFake.t),abfFake);
save('speedThreshold','speedThreshold');
saveas(gcf,'speed.fig');
saveas(gcf,'speed.tif');
close

imageStartNumber=1;
imageEndNumber=length(abfFake.t);
abfStartNumber=imageStartNumber;
abfEndNumber=imageEndNumber;
cellStart=1;
cellEnd=size(dfof,2);

disp('plot run by run')
plotPCAICA_dfof(imageStartNumber,imageEndNumber,cellStart,cellEnd,abfStartNumber,abfEndNumber,binWidth,trackStart,trackEnd,speedThreshold,dfof,abfFake);
saveas(gcf,'allroiactivity.fig');
saveas(gcf,'allroiactivity.tif');
close

% get dfofaveragesmooth
disp('calculate dfofaveragesmooth')
clear dfofaveragesmooth
for n=1:size(dfof,2)
    [dfofaveragesmooth(:,n)]=get_dfofaveragesmooth(n,imageStartNumber,imageEndNumber,abfStartNumber,abfEndNumber,trackStart,trackEnd,binWidth,speedThreshold,dfof,abfFake);
end
filenamedfofA = sprintf('%s_%d_%s.mat', 'dfofaveragesmooth', size(F,2), 'cells');
save(filenamedfofA ,'dfofaveragesmooth');


%% Calculate significant transients

disp('calculate significant transients')
dfof_sig=zeros(size(dfof));
fpr=[];
t_on_off={};
thresh_off_sigma=0.5;
thresh_on_sigma=[];
thresh_len_sec=[];
fprThresh=0.01;

for n=1:size(dfof,2)
        [ thresh_on_sigma(n),thresh_len_sec(n) ] = getSigTranThreshold( F(:,n),dfof(:,n),abfFake.t,fprThresh,0 );
    [ dfof_sig(:,n),fpr(n,1),t_on_off{n},sigAmp{n} ] = getSigTrans( F(:,n),dfof(:,n),abfFake.t,thresh_on_sigma(n),thresh_off_sigma,thresh_len_sec(n),0 );
end

filenamedfofS = sprintf('%s_%d_%s.mat', 'dfof_sig', size(F,2), 'cells');
save(filenamedfofS ,'dfof_sig');
save('sigTransInfo.mat','fpr','t_on_off','thresh_on_sigma','thresh_len_sec','sigAmp');

%get frequency of significant transients
trackLength=trackEnd-trackStart;
for n=1:size(dfof,2)
    t_on_offThisCell=t_on_off{n};
    [ Freq(n,1) ] = sigTransFrequency( t_on_offThisCell,abfFake,trackLength,abfStartNumber,abfEndNumber );
end

%plot siginificant only run by run
disp('plot run by run_significant transients only');

plotPCAICA_dfof(imageStartNumber,imageEndNumber,cellStart,cellEnd,abfStartNumber,abfEndNumber,binWidth,trackStart,trackEnd,speedThreshold,dfof_sig,abfFake);
saveas(gcf,'allroiactivity_sigTrans.fig');
saveas(gcf,'allroiactivity_sigTrans.tif');
close

%compute dfofaveragesmooth_sig
disp('calculate dfofaveragesmooth_sigTrans')
clear dfofaveragesmooth_sig
for n=1:size(dfof,2)
    [dfofaveragesmooth_sig(:,n)]=get_dfofaveragesmooth(n,imageStartNumber,imageEndNumber,abfStartNumber,abfEndNumber,trackStart,trackEnd,binWidth,speedThreshold,dfof_sig,abfFake);
end
filenamedfofAS = sprintf('%s_%d_%s.mat', 'dfofaveragesmooth_sig', size(F,2), 'cells');
save(filenamedfofAS ,'dfofaveragesmooth_sig');


end
