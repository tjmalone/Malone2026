%% plotExampleRBR

clear; close all; clc

p1 = 'D:\AD_Project\imagingData';
cd(p1)

binWidth=5;
trackStart = 0;
trackEnd = 600;

load('foldersLearning.mat')

% set save folder name
svFile = [p1 '/Figures/Figures' char(datetime('today','format','yyyyMMdd')) '/activityExamples'];
if ~isfolder(svFile)
    mkdir(svFile);
    mkdir([svFile '/data']);
end


%% Define example data

exFOV = [23, 21, 31, 8];
exTypes = {'WT-Male','PS19-Male','WT-Female','PS19-Female'};
exCells = [11, 8, 7, 13];
nEx = length(exFOV);

% WT-Male: 23-11
% PS19-Male: 21-16
% WT-Female: 31-7
% PS19-Female: 7-18

% define use days
useDays = 1:11;
nDays = length(useDays);


%% Plot high RBR consistency

figure
tiledlayout(nEx,nDays+1)

trueRBR = zeros(nEx,nDays);

for ff = 1:nEx
    curFOV = exFOV(ff);
    for dd = 1:nDays
        %%
        curDay = useDays(dd);

        curFolderName = foldersLearning{curFOV}{curDay};
        cd(curFolderName)

        % load dfof
        d = dir('dfof_sig_*_cells.mat');
        load(d(1).name)

        % load abf and calculate speed threshold
        load('abf.mat','abf')
        [speedThreshold] = speedThreshold1D(1,length(abf.t),abf,0);

        % select current cell
        curIdx = alignsLearning{curFOV}(exCells(ff),curDay);
        dfof_sig = dfof_sig(:,curIdx);

        % plot binned dfof run by run
        nexttile(dd+(ff-1)*(nDays+1))
        dfofM = plotCurrentDFOF(dfof_sig,abf,trackStart,trackEnd,binWidth,speedThreshold);
        % axis('square')

        % calculate within-day conistency for cell
        cellMean = mean(dfofM,1,'omitnan')';
        corrInfo = dfofMCorrelationSubset({dfofM},cellMean);
        trueRBR(ff,dd) = corrInfo.meantoOthers;

        title(num2str(size(dfofM,1)))
    end

    % plot daily correlation
    nexttile(ff*(nDays+1))
    plot(trueRBR(ff,:))

    % set limits
    xlim([0.5 11.5])
    ylim([0 1])
    set(gca,'XTick',[])
    set(gca,'YTick',[])
end

cd(p1)

% save figure
savefig([svFile '/intraday_examples.fig'])


%% Plot individual cell dfof

function dfofMInterp = plotCurrentDFOF(dfof_sig,abfFake,trackStart,trackEnd,binWidth,speedThreshold)
% adapted from plotdfofrunbyrun. Detailed commenting is found in that
% function

% bin y positions
trackLength=trackEnd-trackStart;
track = trackStart:binWidth:trackEnd;
[~,nbin]=histc(abfFake.y,track);

% apply speed threshold
speed=diff(abfFake.y);
fastEnough=[false ;speed>speedThreshold];

% identify teleports
speed=diff(abfFake.y);
[row,col]=find(speed<-100);

% identify start and end indices
endIndicesP=[row' abfFake.imageIndex(end)]';
startIndicesP=[1 (row+1)']';
startIndices=startIndicesP(diff(startIndicesP)>1);
startIndices=[startIndices' startIndicesP(end)]';
endIndices=endIndicesP(diff([1 endIndicesP'])>1);

% bin dfof and perform guassian smooth
dfofM=zeros(size(endIndices,1),trackLength/binWidth);
for r=1:size(endIndices,1);
    FcuseThisRun=dfof_sig(startIndices(r,1):endIndices(r,1));
    fastEnoughThisRun=fastEnough(startIndices(r,1):endIndices(r,1));
    nThisRun=nbin(startIndices(r,1):endIndices(r,1));
    for binNumber=1:trackLength/binWidth;
        dfofM(r,binNumber)=mean(FcuseThisRun(nThisRun==binNumber&fastEnoughThisRun));
    end
    gaussianWindow=gausswin(3,1);
    [dfofM(r,:)]=gaussianSmoothWithNan(dfofM(r,:),gaussianWindow);
end

dfofMInterp = getDfofMInterp(dfofM);

% plot activity
dfofMPlot = dfofMInterp;
imAlpha = ones(size(dfofMPlot));
imAlpha(isnan(dfofMPlot))=0;
imagesc(dfofMPlot,'AlphaData',imAlpha);
set(gca,'color',0*[1 1 1]);
clim([0 max(dfofMPlot,[],'all')])

% add labels
set(gca,'XTick',[])
set(gca,'YTick',[])

end

function dfofMInterp = getDfofMInterp(dfofM)
% based on dfofMCorrelationFull

% find if end is nan
nanEnd = all(dfofM(end,end-1:end));

dfofMInterp = [];

for m=2:size(dfofM,1)
    % include last lap if one of last two bins is non-nan (TM 2/22/24)
    if m==size(dfofM,1) && nanEnd
        continue
    end
    if length(find(~isnan(dfofM(m,:))))>1 %at least two numbers are not nan
        dfofMInterp(m-1,:)=naninterp(dfofM(m,:));
        %if there are less than 2 non-nan numbers, keep nan
    else
        dfofMInterp(m-1,:)=dfofM(m,:);
    end
end
end

function X = naninterp(X)
% Interpolate over NaNs
% See INTERP1 for more info
X(isnan(X)) = interp1(find(~isnan(X)), X(~isnan(X)), find(isnan(X)),'PCHIP');
end