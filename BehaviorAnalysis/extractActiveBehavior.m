function behavior = extractActiveBehavior(params,envParams,useLaps)
%% extractActiveBehavior
% Calculates basic behavioral parameters for a virmen training log. Called
% in behaviorAD_all.
%
% Inputs:
%   params - virmen log data (output of logParams.m)
%   envParams - environment parameters (reward zone, false alarm zone, and
%       track length)
%   useLaps - range of laps to be analzyed
%

% set bin width
binWidth = 5;

% set stop cluster size
szCluster = 5;

% read environment parameters
rewStart = envParams.rewStart;
rewStop = envParams.rewStop;
rewWin = rewStop-rewStart;
falseStart = envParams.falseStart;
falseStop = falseStart+rewWin;
trackLength = envParams.trackLength;
yBins = 0:binWidth:trackLength;
nBins = length(yBins);

% calculate rolling d-prime paramters
yStarts = 1:trackLength-rewWin;
yMids = yStarts+rewWin/2;
nDPrime = length(yStarts);
goalIdx = find(yStarts>=rewStart,1);

% read use laps
if isempty(useLaps) || isequal(useLaps,0)
    startLap = 1;
    stopLap = 0;
else
    startLap = useLaps(1);
    stopLap = useLaps(end);
end

% set analysis parameters
speedThresh = [1 5];    % speed threshold [stop try]
timeThresh = [60 30];   % stop threshold for smoothing (sec*60) [stop try]
winBuff = [0 25];       % zone window buffer [stop try]
postRewardTh = 30;      % post reward threshold for removing stops


%% Calculate rolling speed

% extract behavior data
t = params.t;
y = params.y;
rIdx = params.rewardIdx;

% calculate number of runs
[runIdx,nRuns] = identifyLaps(y,startLap,stopLap);

% segment data by lap
yLaps = cell(nRuns,1);
tLaps = cell(nRuns,1);
rIdxLaps = zeros(nRuns,1);
for m = 1:nRuns
    yLaps{m} = y(runIdx(m,1):runIdx(m,2));
    tLaps{m} = t(runIdx(m,1):runIdx(m,2));

    for n = 1:length(rIdx)
        if inrange(rIdx(n),runIdx(m,:))
            rIdxLaps(m) = rIdx(n)-runIdx(m,1)+1;
            continue
        end
    end
end

% calculate velocity
vel = cellfun(@(a,b) diff(a)./diff(b),yLaps,tLaps,'UniformOutput',0);


%% Calculate stop/try statistics

% function to limit d-prime range (normally undefined if inputs are 0 or 1)
constrainRange = @(x) max(min(mean(x),1-1/(length(x)*2)),1/(length(x)*2));
constrainRange2D = @(x) max(min(mean(x,1),1-1/(size(x,1)*2)),1/(size(x,1)*2));

% stop versus try
for ii = 1:2
    %% Identify stops
    % An additional stop is added and then removed from each lap to prevent
    % errors caused by laps with no stops.
    
    % velocity with time threshold averaging
    velTh = cellfun(@(a) movmean(abs(a),[timeThresh(ii)-1 0]),vel,...
        'UniformOutput',0);
    velTh = cellfun(@(a) [(speedThresh(ii)+1)*ones(60,1);a(61:end)],...
        velTh,'UniformOutput',0);
    
    % time steps
    velDiff = cellfun(@diff,tLaps,'UniformOutput',0);

    % all stop indices
    velThStop = cellfun(@(a) a<=speedThresh(ii),velTh,'UniformOutput',0);

    % all stop index lengths
    velStopLen = cellfun(@(a,b) [0; a.*b],velDiff,velThStop,'UniformOutput',0);
    velStopLenCat = cat(1,velStopLen{:});

    % y bin along track
    yLapCat = cat(1,yLaps{:});
    yDisc = discretize(yLapCat,yBins);

    % identify the start index of all stops
    idxStopAll = cellfun(@(a) contiguous(a<=speedThresh(ii),1),velTh,'UniformOutput',0);
    idxStopEnds = cellfun(@(a) [[1 1]; a{1,2}],idxStopAll,'UniformOutput',0);

    % identfy start index of all stops
    idxStopStart = cellfun(@(a) a(:,1),idxStopEnds,'UniformOutput',0);

    % identify the start position of all stops
    yStop = cellfun(@(a,b) b(a),idxStopStart,yLaps,'UniformOutput',0);
    yStop = cellfun(@(a) a(2:end),yStop,'UniformOutput',0);
    yStopCat = vertcat(yStop{:});

    % stop times along track
    stopTimeDist = zeros(1,nBins);
    sumStop = sum(velStopLenCat);
    for jj = 1:nBins
        stopTimeDist(jj) = sum(velStopLenCat(yDisc==jj));
        % stopTimeDist(jj) = sum(velStopLenCat(yDisc==jj))/sumStop;
    end
    
    % post-reward y positions
    postRewY = cell(nRuns,1);
    for jj = 1:nRuns
        postRewY{jj} = zeros(size(yLaps{jj}));

        % skip no reward laps
        if rIdxLaps(jj)~=0
            % find current reward index
            rIdxCur = idxStopStart{jj}(find(inrange(yStop{jj},[rewStart rewStop]),1)+1);
            if isempty(rIdxCur)
                rIdxCur = rIdxLaps(jj);
            else
                rIdxCur = min(rIdxCur,rIdxLaps(jj));
            end

            % calculate post-reward indices
            yRew = yLaps{jj}(rIdxCur);
            L1 = 1:length(yLaps{jj})>rIdxCur;
            L2 = yLaps{jj}'<=yRew+postRewardTh;
            rPostIdx = find(L1 & L2,1,'last');
            
            % store post-reward indices
            postRewY{jj}(rIdxCur+1:rPostIdx) = 1;
        end
    end

    % stop times along track (removing post-reward zone)
    postRewYCat = cat(1,postRewY{:});
    velStopLen_NPR = velStopLenCat;
    velStopLen_NPR(postRewYCat==1) = 0;

    % stop times along track
    stopTimeDist_NPR = zeros(1,nBins);
    sumStop_NPR = sum(velStopLen_NPR);
    for jj = 1:nBins
        % stopTimeDist_NPR(jj) = sum(velStopLen_NPR(yDisc==jj))/sumStop_NPR;
        stopTimeDist_NPR(jj) = sum(velStopLen_NPR(yDisc==jj));
    end

    % calculate thresholded start index of all stops
    idxStopStart_NPR = idxStopStart;
    for jj = 1:nRuns
        for kk = 1:length(idxStopStart_NPR{jj})
            if postRewY{jj}(idxStopStart_NPR{jj}(kk))==1
                idxStopStart_NPR{jj}(kk) = 0;
            end
        end
        idxStopStart_NPR{jj}(idxStopStart_NPR{jj}==0) = [];
    end
    
    % identify the start position of all stops with no post-reward stops
    yStop_NPR = cellfun(@(a,b) b(a),idxStopStart_NPR,yLaps,'UniformOutput',0);
    yStop_NPR = cellfun(@(a) a(2:end),yStop_NPR,'UniformOutput',0);
    
    % calculate total number of stops
    numStops = length(yStopCat);
    
    % calculate number of stops per lap
    numStopsPerLap = cellfun('length',yStop);

    numStopClustersPerLap = zeros(nRuns,1);
    for jj= 1:nRuns
        % sort y positions
        [L] = sort(yStop{jj});

        % iteratively remove next lick if is within minLickDist
        kk = 1;
        while kk<length(L)
            if L(kk+1)-L(kk)>szCluster
                kk = kk+1;
            else
                L(kk+1) = [];
            end
        end

        numStopClustersPerLap(jj) = length(L);
    end

    % calculate stop distribution
    if ~isempty(yStopCat)
        stopDensity = ksdensity(yStopCat,yBins);
        stopDist = stopDensity*100;
    else
        stopDist = zeros(size(yBins));
    end

    % find unique discretization of the start position of all stops (max of
    % 1 stop per lap per bin
    yStartDiscIdx = zeros(nRuns,length(yBins));
    for jj = 1:nRuns
        curLap = yStop{jj};
        yStartDisc = unique(discretize(curLap,yBins));
        yStartDisc(isnan(yStartDisc)) = [];
        yStartDiscIdx(jj,yStartDisc) = 1;
    end

    % find mean start of stop distributions
    yStartDistribution = mean(yStartDiscIdx,1);


    % find unique discretization of the start position of all stops (max of
    % 1 stop per lap per bin
    yStartDiscIdx_NPR = zeros(nRuns,length(yBins));
    for jj = 1:nRuns
        curLap = yStop_NPR{jj};
        yStartDisc_NPR = unique(discretize(curLap,yBins));
        yStartDisc_NPR(isnan(yStartDisc_NPR)) = [];
        yStartDiscIdx_NPR(jj,yStartDisc_NPR) = 1;
    end

    % find mean start of stop distributions
    yStartDistribution_NPR = mean(yStartDiscIdx_NPR,1);
    
    % idenitfy all indices in true stops
    idxStopTrueAll = cellfun(@(a) arrayfun(@(b,c) b:c,a(:,1),a(:,2),...
        'UniformOutput',0),idxStopEnds,'UniformOutput',0);
    
    % idenitfy all y positions in true stops  (+1 tries to account for index shifting by diff)
    yStopAll = cellfun(@(a,b) cellfun(@(c) b(c+1),a,...
        'UniformOutput',0),idxStopTrueAll,yLaps,'UniformOutput',0);
    yStopAll = cellfun(@(a) a(2:end),yStopAll,'UniformOutput',0);
    
    % calculate number of stops in reward zone
    numStopsInReward = cellfun(@(a) sum(cellfun(@(b) any(inrange(b,...
        [rewStart-winBuff(ii) rewStop+winBuff(ii)])),a)),yStopAll);
    
    % calculate percent of stops in reward zone
    perStopsInReward = 100*sum(numStopsInReward)/sum(numStopsPerLap);
    
    
    %% Caclulate statistics on all laps
    
    % identify all success and false alarm laps
    success = numStopsInReward>0;
    false = cellfun(@(a) any(cellfun(@(b) any(inrange(b,...
        [falseStart-winBuff(ii) falseStop+winBuff(ii)])),a)),yStopAll);
    
    % calculate percent success
    perSuccess = mean(success)*100;
    
    % calculate dprime
    conSuccess = constrainRange(success);
    conFalse = constrainRange(false);
    dprime = norminv(conSuccess)-norminv(conFalse);
    
    
    %% Caclulate statistics on laps with stops
    
    % idenify laps with stops
    lapsWithStops = numStopsPerLap>0;
    
    % identify all success and false alarm laps
    successWithStops = success(lapsWithStops);
    falseWithStops = false(lapsWithStops);
    
    % calculate percent success
    perSuccessWithStops = mean(successWithStops);
    
    % calculate dprime
    conSuccess = constrainRange(successWithStops);
    conFalse = constrainRange(falseWithStops);
    dprimeWithStops = norminv(conSuccess)-norminv(conFalse);
    

    %% Calculate rolling d-prime
    
    % calculate inclusion indices
    excludeIdx1 = inrange(yStarts,[rewStart-winBuff(ii) rewStop+winBuff(ii)]);
    excludeIdx2 = inrange(yStarts+rewWin,[rewStart-winBuff(ii) rewStop+winBuff(ii)]);
    includeIdx = ~excludeIdx1 & ~excludeIdx2;

    % number of stops in a rolling window
    numStopsRoll = zeros(nRuns,nDPrime);

    for jj = 1:length(yStarts)
        curStart = yStarts(jj);
        curEnd = yStarts(jj)+rewWin;

        % calculate number of stops in current
        for kk = 1:length(yStopAll)
            curStops = zeros(1,length(yStopAll{kk}));
            for mm = 1:length(yStopAll{kk})
                curStops(mm) = any(inrange(yStopAll{kk}{mm},[curStart-winBuff(ii) curEnd+winBuff(ii)]));
            end
            numStopsRoll(kk,jj) = sum(curStops);
        end
    end

    % identify all success laps
    successRoll = numStopsRoll>0;
    successRollMean = mean(successRoll,1);
    
    % calculate constrained percent success
    conSuccessRoll = constrainRange2D(successRoll);

    % calculate d-prime
    dprimeRoll = norminv(conSuccessRoll(goalIdx))-norminv(conSuccessRoll(includeIdx));
    
    % calculate mean d-prime
    dprimeGlobal = mean(dprimeRoll,'omitnan');

    % calculate success percentile
    successTile = prctileInv(successRollMean(includeIdx),successRollMean(goalIdx));
    
    
    %% Save results
    
    % results struct
    dataCur = struct();
    dataCur.yStop = yStop;
    dataCur.stopTimeDist = stopTimeDist;
    dataCur.stopTimeDist_NPR = stopTimeDist_NPR;
    dataCur.yStartDistribution = yStartDistribution;
    dataCur.yStartDistribution_NPR = yStartDistribution_NPR;
    dataCur.numStops = numStops;
    dataCur.numStopsPerLap = numStopsPerLap;
    dataCur.numStopClustersPerLap = numStopClustersPerLap;
    dataCur.perStopsInReward = perStopsInReward;
    dataCur.stopDist = stopDist;
    dataCur.successPerLap = success;
    dataCur.perSuccess = perSuccess;
    dataCur.successRoll = successRoll;
    dataCur.dprime = dprime;
    dataCur.perSuccessWithStops = perSuccessWithStops;
    dataCur.dprimeWithStops = dprimeWithStops;
    dataCur.includeIdx = yMids(includeIdx);
    dataCur.dprimeRoll = dprimeRoll;
    dataCur.dprimeGlobal = dprimeGlobal;
    dataCur.successTile = successTile;
    dataCur.idxInclude = includeIdx;
    dataCur.idxGoal = goalIdx;
    
    % store results
    if ii==1
        stopStats = dataCur;
        
        % find true rewards in use laps
        nRewards = sum(inrange(params.rewardIdx,[runIdx(1,1) runIdx(end,2)]));
        
        % true percent correct for all laps (for debugging)
        percentSuccess = nRewards/nRuns*100;
        percentSuccessWithStops = nRewards/sum(successWithStops)*100;
        
    elseif ii==2
        tryStats = dataCur;
    end
end


%% Calculate velocity statistics

velStats = struct();

% velocity save fields
svFields = {'raw','abs'};

% absolute velocity
velAbs = cellfun(@abs,vel,'UniformOutput',0);

% indices when mouse is moving by lap
isMoving = cellfun(@(a) a>speedThresh(1),velAbs,'UniformOutput',0);
isMovingCat = cat(1,isMoving{:});

% percent of time spent stopped
perTimeStopped = (1-sum(isMovingCat)/length(isMovingCat))*100;

% calculate raw velocity per lap
avgVelMovingRawPerLap = cellfun(@(a,b) mean(b(a),'omitnan'),isMoving,vel);
stopStats.avgVelMovingRawPerLap = avgVelMovingRawPerLap;

% raw velocity versus absolute value
for ii = 1:2
    if ii==1
        velUse = vel;
    else
        velUse = velAbs;
    end
    
    % calculate all velocity
    velAll = cat(1,velUse{:});	% concatenated velocity
    meanVelAll = nanmean(velAll);
    
    % calculate moving velocity
    velAllMoving = velAll(isMovingCat);
    meanVelAllMoving = nanmean(velAllMoving);
       
    
    %% Save results
    
    velStats.(svFields{ii}).velAll = velAll;
    velStats.(svFields{ii}).meanVelAll = meanVelAll;
    velStats.(svFields{ii}).velAllMoving = velAllMoving;
    velStats.(svFields{ii}).meanVelAllMoving = meanVelAllMoving;
    
end

% calculate velocity distribution (mean velocity in each bin)
% discretize y positions
yBinIdx = cellfun(@(a) discretize(a,yBins),yLaps,'UniformOutput',0);

%calcualte velocity distribution (all times)
    velDist = zeros(nRuns,length(yBins));
    for jj = 1:length(yBins)
        velDist(:,jj) = cellfun(@(a,b,c) mean(a(b(1:end-1)==jj & c)),vel,yBinIdx,isMoving);
%         velDist(:,jj) = cellfun(@(a,b,c) nanmean(a(b(1:end-1)==jj)),vel,yBinIdx);
    end

% calcilate velocity means
velDistMean = nanmean(velDist,1);
velDistMean(isnan(velDistMean)) = 0;

% calculate velocity distribution (percent of time spent in each bin)
dwellDensity = cellfun(@(a,b) ksdensity(a(b),yBins),yLaps,isMoving,'UniformOutput',0');
dwellDist = cat(1,dwellDensity{:})*100;

% calcilate dwell means
dwellDistMean = nanmean(dwellDist,1);
dwellDistMean(isnan(dwellDistMean)) = 0;

% save results
velStats.perTimeStopped = perTimeStopped;
velStats.velDistMean = velDistMean;
velStats.dwellDistMean = dwellDistMean;


%% Save outputs

behavior = struct();

behavior.stopStats = stopStats;
behavior.tryStats = tryStats;
behavior.velStats = velStats;
behavior.percentSuccess = percentSuccess;
behavior.percentSuccessWithStops = percentSuccessWithStops;

end

