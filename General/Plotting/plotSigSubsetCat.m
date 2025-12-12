function plotSigSubsetCat(cellStart,cellEnd,cat,trackEnd)
%%

binWidth=5;
trackStart=0;

if nargin<4 || isempty(trackEnd)
    trackEnd=1000;
end

load('abfFake.mat','abfFake')

fold = 'byCat\';
load([fold 'dfof_sigCat.mat'],'dfof_sigCat')
dfof_sig = dfof_sigCat{cat};

cellEnd = min(cellEnd,size(dfof_sig,2));

imageStartNumber=1;
imageEndNumber=length(abfFake.t);

abfStartNumber=imageStartNumber;
abfEndNumber=imageEndNumber;

[speedThreshold]= speedThreshold1D(1,length(abfFake.t),abfFake,0);

plotPCAICA_dfof(imageStartNumber,imageEndNumber,cellStart,cellEnd,...
    abfStartNumber,abfEndNumber,binWidth,trackStart,trackEnd,...
    speedThreshold,dfof_sig,abfFake);

