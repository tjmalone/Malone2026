%% plotExampleRBR

clear; close all; clc

p1 = '/MATLAB Drive/FY2025/imagingData';
cd(p1)

binWidth=5;
trackStart = 0;
trackEnd = 600;

% set save folder name
svFile = [p1 '/Figures/Figures' char(datetime('today','format','yyyyMMdd')) '/activityExamples'];
if ~isfolder(svFile)
    mkdir(svFile);
    mkdir([svFile '/data']);
end


%% plot high RBR consistency

cd('D:/AD_Project/imagingData/data/TM240316-J/240519/loc1/TSeries-2388_1/suite2p')

useIdx = 38;

load('dfof_sig_318_cells.mat')
load('abfFake.mat','abfFake')
load('speedThreshold.mat')

dfof_sig = dfof_sig(:,useIdx);

% plot binned dfof run by run
plotCurrentDFOF(dfof_sig,abfFake,trackStart,trackEnd,binWidth,speedThreshold);
c = clim;

% save figure
savefig([svFile '/intraday_highAll.fig'])

% load mean activity
load('dfofaveragesmooth_sig_interp.mat')
dfofMean = dfofaveragesmooth_sig_interp(:,useIdx);

% normalize mean activity
mn = min(dfofMean);
mx = max(dfofMean);
dfofMeanNorm  = (dfofMean-mn)/(mx-mn);

% plot activity
figure
plot(dfofMeanNorm)
ylim([0 1])

% add labels
xLabs = string(get(gca,'XTick')*binWidth);
set(gca,'XTickLabels',xLabs)
xlabel('Track Position')
ylabel('Normalized dfof')

% save figure
savefig([svFile '/intraday_highMean.fig'])


%% plot low RBR consistency

cd('D:/AD_Project/imagingData/data/TM240316-J/240510/loc1/TSeries-2348_2/suite2p/')

useIdx = 105;
useIdx = 140;

load('dfof_sig_379_cells.mat')
load('abfFake.mat','abfFake')
load('speedThreshold.mat')

dfof_sig = dfof_sig(:,useIdx);

% plot bined dfof run by run
plotCurrentDFOF(dfof_sig,abfFake,trackStart,trackEnd,binWidth,speedThreshold);
c = clim;

% save figure
savefig([svFile '/intraday_lowAll.fig'])

% load mean activity
load('dfofaveragesmooth_sig_interp.mat')
dfofMean = dfofaveragesmooth_sig_interp(:,useIdx);

% normalize mean activity
mn = min(dfofMean);
mx = max(dfofMean);
dfofMeanNorm  = (dfofMean-mn)/(mx-mn);

% plot activity
figure
plot(dfofMeanNorm)
ylim([0 1])

% add labels
xLabs = string(get(gca,'XTick')*binWidth);
set(gca,'XTickLabels',xLabs)
xlabel('Track Position')
ylabel('Normalized dfof')

% save figure
savefig([svFile '/intraday_lowMean.fig'])

cd(p1)


%% Plot individual cell dfof

function plotCurrentDFOF(dfof_sig,abfFake,trackStart,trackEnd,binWidth,speedThreshold)
% adapted from plotdfofrunbyrun. Detailed commenting is found in that
% function

% bin y positions
trackLength=trackEnd-trackStart;
track=[trackStart:binWidth:trackEnd];
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

% plot activity
figure
dfofMPlot = dfofM(2:end,:);
imAlpha = ones(size(dfofMPlot));
imAlpha(isnan(dfofMPlot))=0;
imagesc(dfofMPlot,'AlphaData',imAlpha);
set(gca,'color',0*[1 1 1]);
clim([0 max(dfofMPlot,[],'all')])

% add labels
xLabs = string(get(gca,'XTick')*binWidth);
set(gca,'XTickLabels',xLabs)
xlabel('Track Position')
ylabel('Lap')

end