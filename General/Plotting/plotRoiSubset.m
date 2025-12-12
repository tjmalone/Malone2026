function plotRoiSubset(cellStart,cellEnd)
%%

binWidth=5;
trackStart=0;
trackEnd=1000;

d = dir('dfof*.mat');
load(d(1).name)
load('abfFake.mat','abfFake')

imageStartNumber=1;
imageEndNumber=length(abfFake.t);

abfStartNumber=imageStartNumber;
abfEndNumber=imageEndNumber;

[speedThreshold]= speedThreshold1D(1,length(abfFake.t),abfFake);

plotPCAICA_dfof(imageStartNumber,imageEndNumber,cellStart,cellEnd,...
    abfStartNumber,abfEndNumber,binWidth,trackStart,trackEnd,...
    speedThreshold,dfof,abfFake);

