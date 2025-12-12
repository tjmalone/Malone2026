function corrInfo = dfofMCorrelationSubset(dfofMInterpM,dfofaveragesmooth_interp,binStart,binEnd)
%% dfofMCorrelationSubset
% calculate RBR consistency for a given track subset

%% Process inputs

if nargin<3
    binStart = 1;
end

if nargin<4
    binEnd = size(dfofMInterpM{1},2)-binStart+1;
end


%% Calculate correlation info

% initialize correlation information struct
corrInfo = struct();

nCells = length(dfofMInterpM);
nRuns = size(dfofMInterpM{1},1);

% truncate inputs
for n=1:nCells
    dfofMInterpM{n} = dfofMInterpM{n}(:,binStart:binEnd);
end
dfofaveragesmooth_interp = dfofaveragesmooth_interp(binStart:binEnd,:);

% identify runs without nans
corrInfo.noNaN = zeros(nRuns,nCells);
for n=1:nCells
    corrInfo.noNaN(:,n) = ~any(isnan(dfofMInterpM{n}),2);
end

% compute the correlation of each run to the mean dfofaveragesmooth_interp
corrInfo.toMean = zeros(nRuns,nCells);
corrInfo.meantoMean = zeros(1,nCells);
for n=1:nCells
    meanCur = dfofaveragesmooth_interp(:,n);
    corrCur = corr(dfofMInterpM{n}',meanCur);
    corrInfo.toMean(:,n) = corrCur;
    corrInfo.meantoMean(:,n) = mean(corrCur,'omitnan');
end

% compute the correlation for adjescent and all other runs
corrInfo.toNext = zeros(nRuns-1,nCells);
corrInfo.meantoNext = zeros(1,nCells);
corrInfo.toOthers = zeros(nRuns*(nRuns-1)/2,nCells);
corrInfo.meantoOthers = zeros(1,nCells);
for n=1:nCells
    corrCur = corr(dfofMInterpM{n}',dfofMInterpM{n}');

    % adjacent run correlation
    corrNext = diag(corrCur,-1);
    corrInfo.toNext(:,n) = corrNext;
    corrInfo.meantoNext(1,n) = mean(corrNext,'omitnan');

    % other run correlation
    corrOther = corrCur(tril(true(size(corrCur)),-1));
    corrInfo.toOthers(:,n) = corrOther;
    corrInfo.meantoOthers(1,n) = mean(corrOther,'omitnan');
end

end

