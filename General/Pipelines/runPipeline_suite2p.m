%% runPipeline_dual.m
%
% Runs PCAICA pipeline on Z-series imaging data in current directory. Includes
% converting images to a tifstack, motion correction, and pcaica analysis.
% Compatible with automation by runPipeAuto
%

clear all; close all;

% folder where pipeline code is saved
baseFolder = 'D:\AnalysisCode\Pipelines\suite2p\';

% copy self to current directory
if ~exist([pwd '\runPipeline_suite2p.m'],'file')
    copyfile([baseFolder 'runPipeline_suite2p.m'],'runPipeline_suite2p.m');
end

% identify folder with original tiffs
f='TSeries*';
d = dir(f);
folderPath = d(1).name;

% generate Z folder names
outFolders = [folderPath '_1'];

% make internal folders/files
mkdir([outFolders '\tifStack'])
mkdir([outFolders '\suite2p'])
%     copyfile([baseFolder 'PCAICA_resonantGalvoNewNew.m'],[outFolders '\pcaica\PCAICA_resonantGalvoNewNew.m'])


%% Set parameters

limLaps = 0;

filter = '*Ch2*';   % tiff name filter
NPerStack = 1002;   % number of tiffs per stack
dwnSample = 2;      % downsampling factor
cutPix = 1;


%% Load/save glavo trace

thresh=4.8;

f2 = '*VoltageRecording*.csv';
d2=dir([folderPath '\' f2]);

m = csvread([folderPath '\' d2(1).name],2,0);
save([outFolders '\m.mat'],'m','-v7.3');

% identify synchronized points
a=m(:,2)>thresh;
b=contiguous(a,[1]);
c=b{1,2};
middleFrameIdx=ceil(mean(c,2));

% split galvo trace
md = m(middleFrameIdx,:);
fRate = mean(diff(md(:,1)))/1000;   % frame rate (sec)
save([outFolders '\md.mat'],'md','fRate','dwnSample','-v7.3')


%% Identify last frame before max lap number.

%threshold for identifying lap number by voltage
lapThresh = 0.6;

if limLaps>=1
    % max lap value
    maxLap = max(m(middleFrameIdx,4))*100;
    
    % voltage threshold for max lap number
    lastLapThresh = (floor(maxLap+lapBuffer)-lapBuffer)/100;
    
    % index of first frame after last lap
    flIdx = find(m(middleFrameIdx,4)>lastLapThresh,1,'first');
    
    % index of second frame after last lap (so each Z plane has one frame)
    lastIdx = min(length(middleFrameIdx),flIdx+1);
else
    % index of last image frame
    lastIdx = length(middleFrameIdx);
end


%% Convert Images to Tiff Stacks

tic

% NPerStack must be divisible by dwnSample (unnecessary with rolling
% average)
if mod(NPerStack,dwnSample)~=0
    disp('Error: Stack Mismatch')
    return
end

% generate tiff stacks
indivTiffsToStack_single(folderPath,outFolders,filter,NPerStack,dwnSample,lastIdx,cutPix);

toc


%% Run motion correction for each Z plane

%% Start Miji

% identify folder with original tiffs
f='TSeries*';
d = dir(f);
folderPath = d(1).name;

% generate Z folder names
outFolders = [folderPath '_1'];

cd(outFolders)

%motion correction using Andrea's code. modified based on the
%agExampleScriptXMovie.m
%this code add a bit number to the movie so that all the numbers are
%positive.

%run Miji first
% Miji
mkdir('motionCorrected5times');
cd('motionCorrected5times');
Miji

clear all classes
try
    MIJ
    disp('you are good to go')
catch
    error('type Miji before starting')
end


p = pwd;
filePath = [p '\tifStack\'];

cd('motionCorrected5times');


%% Initialize stacks

mkdir('ref');
mkdir('matrix')
mkdir('correctionInfo')

bigNumber=0;
offset=0;
save('offset.mat','offset');
frameRate=.067; % mvie frame rate. unecessary but have to put a number here
%find all files in the same imaging session in the folder
% filePath='D:\data\ID20201118\20210104\2ndloc2\TSeries-01012021-0935-069\tifStack\';
imageName='TS_*.tif';

fn= dir([filePath imageName]);
stackNumber=length(fn);
imPaths = cellfun(@(x)([filePath x]),{fn(1:stackNumber).name},'uniformoutput',0);

p=pwd;

filePathCorrectedRef=[p '\ref\'];
filePathCorrected=[p '\'];
%stay in this folder when do the refrence correction

%% first correct one frame in the middle of the stacks
%stay in the ref folder
%this step just to get a good reference frame doesn't care what is the real
%movie. so do not need add big number to make image number to be consistant
%across stacks
tic
cd('ref');
refStackNum=round(stackNumber/2);
%start to correct that refrence stack
%create a security copy
movOrig=XMovie(imPaths{refStackNum},frameRate); % create movie object (load movie)j
filenameref=sprintf('%s_%03d.tif','refStack',refStackNum);
movOrig.saveToTiff(filenameref) % create security copy to not affect original
[min_,max_,mean_,median_,std_]=movOrig.getStatistics(); % get some basic info about the movie
movOrig=movOrig-(offset);
%after subtracting this offset, there will still be residual negative
%numbers, which are noises. so we should zero all these numbers. the step
%below will do this as well as converting everything to 16 bit (look at
%this website for why it does this:
%http://imagej.1557.x6.nabble.com/How-ImageJ-handle-with-negative-value-during-image-subtraction-td5000706.html)
movOrig=movOrig.setBitDepth(16);

%if you want to increase the SNR you can downsample in the z direction
%     if 1
%         newNumberOfFrames=round(movOrig.numFrames/3);
%         movOrig=movOrig.resize(0,0,newNumberOfFrames);
%     end
newNumberOfFrames = movOrig.numFrames;
movOrig.saveToTiff(filenameref) % save movie to compare

%motion correct the reference movie to generate a reference frame for all
%stacks to use.
xtot=[];
ytot=[];
correlationsRef=[];
maxShift=10;
% numIterations=13;%normally need 10 rounds of motion correction here. this number can be modified based on how the final correction looks like
for kk=1:2
    %use the median of the frame as a reference frame. mean doesn't
    %work!!!!!!!!!!!!!!!!
    imgref=movOrig.zProject('median');
    [movOrig,xshifts,yshifts]=movOrig.IJmotionCorrect(maxShift,imgref);
    correlationsRef(kk,1)=corr2(imgref,movOrig.zProject('median'));
    xtot=[xtot; xshifts];
    ytot=[ytot; yshifts];
end

while correlationsRef(end) - correlationsRef(end-1)> 0.00005;
    imgref=movOrig.zProject('median');
    [movOrig,xshifts,yshifts]=movOrig.IJmotionCorrect(maxShift,imgref);
    correlationsRef(end+1,1)=corr2(imgref,movOrig.zProject('median'));
    xtot=[xtot; xshifts];
    ytot=[ytot; yshifts];
    plot([xtot ytot])
    drawnow
    axis tight
end
%close the figure
close

numIterations=numel(xtot)/newNumberOfFrames;


%after this step, there are three possibilities
%(1) correlations(end)=correlations(end-1);
%(2) correlations(end)<correlations(end-1);
%(3) correlations(end)>correlations(end-1),but the difference is smaller
%than 0.00005;

%in either case, we just wnat to keep the results from the second last correction
correlationsn=[];
xtotn=[];
ytotn=[];

if correlationsRef(end)>correlationsRef(end-1);
    correlationsn=correlationsRef;
    xtotn=xtot;
    ytotn=ytot;
    numIterationsUse=numIterations;
else
    correlationsn=correlationsRef(1:end-1);
    xtotn=xtot(1:end-newNumberOfFrames);
    ytotn=ytot(1:end-newNumberOfFrames);
    numIterationsUse=numIterations-1;
end
% %if do this all corrections will take longer than 5 rounds
% ref2=movOrig.zProject('median');

%close the plot of xtot and ytot
close
%save this corrected movie and
filenamerefmc=sprintf('%s_%03d.tif','refStackMC',refStackNum);
movOrig.saveToTiff(filenamerefmc)
movOrig.movieId.close();

%get information about the correction. and also apply the shifts to movie and save shifts
movMC=XMovie([filePathCorrectedRef filenameref]);
stdShift=[];
xx=[];
yy=[];
eachXX=[];
eachYY=[];
%x or y shift in each round of correction, rows are number of frames,
%each column is one round of motion correction
eachXX=reshape(xtotn,[movMC.numFrames numIterationsUse]);
eachYY=reshape(ytotn,[movMC.numFrames numIterationsUse]);
xx=sum(eachXX,2);
yy=sum(eachYY,2);
for n=1:numIterationsUse;
    stdShift(n,1)=std(eachXX(:,n));
    stdShift(n,2)=std(eachYY(:,n));
end
%plot some info about this correction
figure,
%correlations
subplot(221)
plot([1:1:numIterationsUse],correlationsn);
xlabel('number of corrections');
ylabel('correlations');
axis tight
title('correlations to reference')

%std of xtshift in each correction
subplot(222);
plot([1:1:numIterationsUse],stdShift(:,1),'r');
hold on
plot([1:1:numIterationsUse],stdShift(:,2),'g');
axis tight
title('std of shifts')
xlabel('number of corrections');
ylabel('std');
legend('x shifts','y shifts','Location','Northeast');

%sum of x y shifts
subplot(223);
plot([1:1:movMC.numFrames],xx,'r');
hold on
plot([1:1:movMC.numFrames],yy,'g');
axis tight
title('sum of shifts')
xlabel('frames after downsampling');
ylabel('shifts in pixel');
legend('x shifts','y shifts','Location','Northeast');

%shifts of xy in each round of correction
subplot(224);
plot(xtotn,'r');
hold on
plot(ytotn,'g');
axis tight
title('shifts in each correction')
ylabel('shifts in pixel');
legend('x shifts','y shifts','Location','Northeast');

savefig('infoOfCorrections.fig');
save('shiftsAndCorrsOfRef.mat','correlationsRef','xtot','ytot','xx','yy','eachXX','eachYY')
close

%apply these shift on the original movie
subpixelreg=1; % =1 if you want subpixel registration.DISCLAIMER:it
% will change the fluorescence of pixels because of the
% interpolation
movMC.translateFrame(xx,yy,subpixelreg);
%get the 3D matrix of the movie
data=movMC.getMovie();
%generating referece frame
ref=median(data,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%add later 20201204 for Bruker data
ref=double(ref);
%this is the previous one
% ref=int16(ref);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save the reference frame
save('ref.mat','ref');
save('data.mat','data');

%save
% filenamereffinal=sprintf('%s_%d.tif','refStackMCFinal',refStackNum
filenamereffinal=sprintf('%s_%03d.tif','refStackMCFinal',refStackNum);
movMC.saveToTiff(filenamereffinal)
movMC.movieId.close();

% now use the refrence frame to correct all movie. stay in the correctedMiji folder

% MIJ.exit
% %restart miji
% Miji;
%

cd ../
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%the line below was removed on 20201205
% load('ref/ref.mat');
%somehow reloading ref change the ref file. so do not reload, just use the
%one saved in workplace.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nOriginalFrames=0;
%
%
for n=1:stackNumber;
    %%
    
    disp(['correcting',num2str(n),' / ', num2str(stackNumber)])
    clear movOrig
    clear movMC
    clear newNumberOfFrames
    clear data
    
    movOrig=XMovie(imPaths{n},frameRate); % create movie object (load movie)j
    % filename=sprintf('%s_%d.tif','corrected',n);
    filename=sprintf('%s_%03d.tif','corrected',n);
    movOrig.saveToTiff(filename) % create security copy to not affect original
    [min_,max_,mean_,median_,std_]=movOrig.getStatistics(); % get some basic info about the movie
    movOrig=movOrig-(offset);
    %after subtracting this offset, there will still be residual negative
    %numbers, which are noises. so we should zero all these numbers. the step
    %below will do this as well as converting everything to 16 bit (look at
    %this website for why it does this:
    %http://imagej.1557.x6.nabble.com/How-ImageJ-handle-with-negative-value-during-image-subtraction-td5000706.html)
    movOrig=movOrig.setBitDepth(16);
    
    
    
    %if you want to increase the SNR you can downsample in the z direction
    %         if 1;
    %             newNumberOfFrames=round(movOrig.numFrames/3);
    %             movOrig=movOrig.resize(0,0,newNumberOfFrames);
    %         end
    newNumberOfFrames = movOrig.numFrames;
    movOrig.saveToTiff(filename) % save movie to compare
    %motion correct the reference movie to generate a reference frame for all
    %stacks to use.
    maxShift=10;
    %do 2 rounds of correction first
    correlations=[];
    xtot=[];
    ytot=[];
    numIterations=5;
    for kk=1:numIterations
        [movOrig,xshifts,yshifts]=movOrig.IJmotionCorrect(maxShift,ref);
        correlations(end+1,1)=corr2(ref,movOrig.zProject('median'));
        xtot=[xtot; xshifts];
        ytot=[ytot; yshifts];
        plot([xtot ytot])
        drawnow
        axis tight
    end
    
    close
    %no saving to this corrected movie
    movOrig.close(1);%force to close no changes
    %we only get shifts from this correction and then apply it to the original
    %movie
    
    disp('got xyshifts')
    disp('save correction info')
    %get information about the correction. and also apply the shifts to movie and save shifts
    
    numIterations=numel(xtot)/newNumberOfFrames;
    %after this step, there are three possibilities
    %(1) correlations(end)=correlations(end-1);
    %(2) correlations(end)<correlations(end-1);
    %(3) correlations(end)>correlations(end-1),but the difference is smaller
    %than 0.0001;
    
    %in either case, we just wnat to keep the results from the second last correction
    correlationsn=[];
    xtotn=[];
    ytotn=[];
    
    
    correlationsn=correlations;
    xtotn=xtot;
    ytotn=ytot;
    numIterationsUse=numIterations;
    
    
    stdShift=[];
    xx=[];
    yy=[];
    eachXX=[];
    eachYY=[];
    %x or y shift in each round of correction, rows are number of frames,
    %each column is one round of motion correction
    
    eachXX=reshape(xtotn,[newNumberOfFrames numIterationsUse]);
    eachYY=reshape(ytotn,[newNumberOfFrames numIterationsUse]);
    xx=sum(eachXX,2);
    yy=sum(eachYY,2);
    for m=1:numIterationsUse;
        stdShift(m,1)=std(eachXX(:,m));
        stdShift(m,2)=std(eachYY(:,m));
    end
    %plot some info about this correction
    figure,
    %correlations
    subplot(221)
    plot([1:1:numIterationsUse],correlationsn);
    xlabel('number of corrections');
    ylabel('correlations');
    axis tight
    title('correlations to reference')
    
    %std of xtshift in each correction
    subplot(222);
    plot([1:1:numIterationsUse],stdShift(:,1),'r');
    hold on
    plot([1:1:numIterationsUse],stdShift(:,2),'g');
    axis tight
    title('std of shifts')
    xlabel('number of corrections');
    ylabel('std');
    legend('x shifts','y shifts','Location','Northeast');
    
    %sum of x y shifts
    subplot(223);
    plot([1:1:newNumberOfFrames],xx,'r');
    hold on
    plot([1:1:newNumberOfFrames],yy,'g');
    axis tight
    title('sum of shifts')
    xlabel('frames after downsampling');
    ylabel('shifts in pixel');
    legend('x shifts','y shifts','Location','Northeast');
    
    %shifts of xy in each round of correction
    subplot(224);
    plot(xtotn,'r');
    hold on
    plot(ytotn,'g');
    axis tight
    title('shifts in each correction')
    ylabel('shifts in pixel');
    legend('x shifts','y shifts','Location','Northeast');
    
    % filenameFigure=sprintf('%s_%d.fig','correctionInfo/CorrectionFigure',n);
    filenameFigure=sprintf('%s_%03d.fig','correctionInfo/CorrectionFigure',n);
    savefig(filenameFigure);
    close
    % filenameInfo=sprintf('%s_%d.mat','correctionInfo/CorrectionInfo',n);
    filenameInfo=sprintf('%s_%03d.mat','correctionInfo/CorrectionInfo',n);
    save(filenameInfo,'correlationsn','xtotn','ytotn','xx','yy','eachXX','eachYY')
    
    disp('apply shifts to original downsampled movie')
    %load the original movie again
    movMC=XMovie(imPaths{n},frameRate);
    movMC.saveToTiff(filename)%replace the movie before
    %subtract the big number
    movMC=movMC-bigNumber;
    %downsample
    if 1;
        movMC=movMC.resize(0,0,newNumberOfFrames);
    end
    
    %apply these shift on the original movie
    subpixelreg=1; % =1 if you want subpixel registration.DISCLAIMER:it
    % will change the fluorescence of pixels because of the
    % interpolation
    movMC.translateFrame(xx,yy,subpixelreg);
    movMC.saveToTiff(filename)
    
    %read and save the 3dmatrix of the movie
    data=movMC.getMovie();
    %this is single of the data just to save space. The data is a matrix with all 0 around the edge and big positive numbers
    %at the center (e.g. start from 850..). We should be able to subtract this
    %1000 (bigNumber) in the following analysis.
    data=single(data);
    % filenameData=sprintf('%s_%d.mat','matrix/matrix',n);
    %filenameData=sprintf('%s_%03d.mat','matrix/matrix',n);
    %save(filenameData,'data');
    
    nOriginalFrames=nOriginalFrames+size(data,3);
    movMC.movieId.close();
    disp('done')
    
    save('nOriginalFrames.mat','nOriginalFrames');
    
end

%%

toc;
clear all
close all
MIJ.exit

correctMCArea()

cd ..\..


%% Reestablish folders

% identify folder with original tiffs
f = 'TSeries*';
d = dir(f);
folderPath = d(1).name;

% generate Z folder names
outFolders = [folderPath '_1'];


%% Generate abf

cd(outFolders)

% nFrameDetected is the number of volatge trace time points
load('md.mat')
nFrameDetected = size(md,1);

% nOriginalFrames is the number of imaging frames after motion correction
load('motionCorrected5times/nOriginalFrames.mat');

if nOriginalFrames>nFrameDetected;
    nOriginalFrames=nFrameDetected;
end

abf.imageIndex = (1:1:nOriginalFrames)';
abf.t = fRate*(abf.imageIndex-1);

abf.y = zeros(nOriginalFrames,1);
for jj = 1:nOriginalFrames
    abf.y(jj) = mean(md(jj:jj+dwnSample-1,3)*100);
end

cd('suite2p')
save('abf.mat','abf');

cd ..\..


%% Plot behavior
%5th column is lick
%6th column is reward

load([outFolders{1} '\m.mat'])

% identify synchronized points
thresh=4.8;
a=m(:,2)>thresh;
b=contiguous(a,1);
c=b{1,2};
middleFrameIdx=ceil(mean(c,2));

% calculate all imaging indices
cImaging=b{1,2};

% calculate all lick indices
lickThresh=0.8;
lick=m(:,5)>lickThresh;
lickAllIdx=find(lick);

% calculate all reward indices
rewardThresh=0.8;
reward=m(:,6)>rewardThresh;
b=contiguous(reward,1);
c=b{1,2};
if ~isempty(c)
    rewardAllIdx = c(:,1);
else
    rewardAllIdx = c;
end

% rewardAllIdx=find(reward);

% cycle through Z planes
%% Initialize

cd(outFolders)

% load abf
load('suite2p\abf.mat')

nOriginalFrames = length(abf.imageIndex);

curIdx = middleFrameIdx(ii:2:end);

% calculate current imaging indices
cImagingFrame = [];
cImagingIdx = [];
for n = ii:2:nOriginalFrames*2
    cImagingFrame = [cImagingFrame; (cImaging(n,1):cImaging(n,2))'];
    cImagingIdx = [cImagingIdx; n*ones(diff(cImaging(n,:))+1,1)];
end


%% Calculate licks
% Calculate the overlap betweeen the lick trace and imaging trace. Uses
% all lick indices (not just those that coincide with peak lick signal).
% lickIdx - all lick indices matched to the nearest imaging frame.
% lickFrame - lick indices that overlap an imaging frame

% calculate nearest image frame to licks
closestIndex = zeros(length(lickAllIdx),1);
for n = 1:length(lickAllIdx)
    [~,closestIndex(n)] = min(abs(lickAllIdx(n)-curIdx));
end

% calculate unique lick indices up to max frame number
closestIndexUnq = unique(closestIndex);
closestIndexUnq(closestIndexUnq>=nOriginalFrames) = [];
closestIndexUnq(closestIndexUnq<1) = [];

% save lick indices (all lick indices matched to the nearest frame)
lickIdx = false(nOriginalFrames,1);
lickIdx(closestIndexUnq) = 1;
abf.lickIdx = lickIdx;

% calculate lick overlaps
[~,~,lickInt] = intersect(lickAllIdx,cImagingFrame);
lickFrameIdx = unique(cImagingIdx(lickInt));

% save lick frames (lick indices that overlap an imaging frame)
lickFrame = zeros(nOriginalFrames,1);
lickFrame(lickFrameIdx) = 1;
abf.lickFrame = lickFrame;


%% Calculate rewards
% Calculate the overlap betweeen the reward trace and imaging trace.
% Uses all reward indices (not just those that coincide with start of
% reward signal).
% rewardIdx - all reward indices matched to the nearest imaging frame.
% rewardFrame - reward indices that overlap an imaging frame

% calculate nearest image frame to rewards
closestIndex = zeros(length(rewardAllIdx),1);
for n = 1:length(rewardAllIdx)
    [~,closestIndex(n)] = min(abs(rewardAllIdx(n)-curIdx));
end

% calculate unique reward indices up to max frame number
closestIndexUnq = unique(closestIndex);
closestIndexUnq(closestIndexUnq>=nOriginalFrames) = [];
closestIndexUnq(closestIndexUnq<1) = [];

% save reward indices (all reward indices matched to the nearest frame)
rewardIdx = false(nOriginalFrames,1);
rewardIdx(closestIndexUnq) = 1;
abf.rewardIdx = rewardIdx;

% calculate reward overlaps
[~,~,rewardInt] = intersect(rewardAllIdx,cImagingFrame);
rewardFrameIdx = unique(cImagingIdx(rewardInt));

% save reward frames (reward indices that overlap an imaging frame)
rewardFrame = zeros(nOriginalFrames,1);
rewardFrame(rewardFrameIdx) = 1;
abf.rewardFrame = rewardFrame;


%% Plot behavior

% calculate y positions
ylick = abf.y(lickIdx);
yreward = abf.y(rewardIdx);

% calculate times
t = abf.t;
tlick = abf.t(lickIdx);
treward = abf.t(rewardIdx);

scaleX = 1;
figure; hold on
plot(t(1:scaleX:end,1),abf.y(1:scaleX:end),'k');
plot(tlick,ylick,'r.','MarkerSize',15);
plot(treward,yreward,'g.','MarkerSize',15);

saveas(gcf,'suite2p\behavior.fig');
close


%% Save

save('suite2p\abf.mat','abf')

abfFake = abf;
save('suite2p\abfFake.mat','abfFake')

cd ..



toc