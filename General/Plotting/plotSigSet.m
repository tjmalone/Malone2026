function plotSigSet(cellSet,trackEnd,dim,rat)
%%

binWidth=5;
trackStart=0;

if nargin<2 || isempty(trackEnd)
    trackEnd=1000;
end

d = dir('dfof_sig*.mat');
load(d(1).name,'dfof_sig')
load('abfFake.mat','abfFake')

dfof_sig = dfof_sig(:,cellSet);

imageStartNumber=1;
imageEndNumber=length(abfFake.t);

cellStart = 1;
cellEnd = length(cellSet);

abfStartNumber=imageStartNumber;
abfEndNumber=imageEndNumber;

[speedThreshold]= speedThreshold1D(1,length(abfFake.t),abfFake,0);

plotPCAICA_dfof(imageStartNumber,imageEndNumber,cellStart,cellEnd,...
    abfStartNumber,abfEndNumber,binWidth,trackStart,trackEnd,...
    speedThreshold,dfof_sig,abfFake);


%%

if nargin>=3
    hFig    = gcf;
    hAxes   = findobj(allchild(hFig), 'flat', 'Type', 'axes');
    axesPos = get(hAxes, 'Position');
    
    axesPos = cat(1,axesPos{:});
    xStart = min(axesPos(:,1));
    xStop = max(axesPos(:,1)) + axesPos(1,3);
    xRng = xStop-xStart;
    
    if dim(2)==1
        xRat = 1;
    else
        if nargin>=4
            xRat = rat(2);
        else
            xRat = 0.8;
        end
    end
    newXStart = xStart;
    xSz = xRng*xRat/dim(2);
    if dim(2)==1
        xStep = 1;
    else
        xStep = xSz + xRng*(1-xRat)/(dim(2)-1);
    end
    
    yStart = min(axesPos(:,2)) - axesPos(1,4);
    yStop = max(axesPos(:,2));
    yRng = yStop-yStart;
    
    if dim(1)==1
        yRat = 1;
    else
        if nargin>=4
            yRat = rat(1);
        else
            yRat = 0.8;
        end
    end
    newYStart = yStart + axesPos(1,4);
    ySz = yRng*yRat/dim(1);
    if dim(1)==1
        yStep = 1;
    else
        yStep = ySz + yRng*(1-yRat)/(dim(1)-1);
    end
    
    for ii = 1:dim(2)
        for jj = 1:dim(1)
            if dim(1)*(ii-1)+jj <= size(axesPos,1)
            newPos = [newXStart+(ii-1)*xStep, newYStart+(jj-1)*yStep, xSz, ySz];
            
            set(hAxes(dim(1)*(ii-1)+jj), 'Position', newPos);
            end
        end
    end
end

end
