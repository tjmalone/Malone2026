function plotSigSubset(cellStart,cellEnd,trackEnd)
%%

binWidth=5;
trackStart=0;

if nargin<3 || isempty(trackEnd)
    trackEnd=400;
end

d = dir('dfof_sig*.mat');
load(d(1).name)
load('abfFake.mat','abfFake')

imageStartNumber=1;
imageEndNumber=length(abfFake.t);

abfStartNumber=imageStartNumber;
abfEndNumber=imageEndNumber;

[speedThreshold]= speedThreshold1D(1,length(abfFake.t),abfFake,0);

nCellMax = size(dfof_sig,2);
cellStart = min(cellStart,nCellMax);
cellEnd = min(cellEnd,nCellMax);

plotPCAICA_dfof(imageStartNumber,imageEndNumber,cellStart,cellEnd,...
    abfStartNumber,abfEndNumber,binWidth,trackStart,trackEnd,...
    speedThreshold,dfof_sig,abfFake);

end
