function dfofaveragesmooth_interp = get_dfofaveragesmoothSubset(dfof,abf,endIdx,trackLength,binWidth)

% load data
load('speedThreshold.mat','speedThreshold');

% set parameters
imStartN = 1;
imEndN = endIdx;
cellN = size(dfof,2);
trackStart = 0;

% get dfofaveragesmooth_interp
dfofaveragesmooth_interp = zeros(trackLength/binWidth,cellN);

for n=1:cellN
    A = get_dfofaveragesmooth(n,imStartN,imEndN,imStartN,imEndN,...
        trackStart,trackLength,binWidth,speedThreshold,dfof,abf);

    if sum(~isnan(A))>1%at least two numbers are not nan
        dfofaveragesmooth_interp(:,n)=naninterp(A);
    else
        A(isnan(A))=0;
        dfofaveragesmooth_interp(:,n)=A;
    end
end