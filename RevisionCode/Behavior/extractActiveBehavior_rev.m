function stopStats = extractActiveBehavior_rev(params,envParams,useLaps,rewCut)
%% extractActiveBehavior
% Calculates basic behavioral parameters for a virmen training log. Called
% in behaviorAD_all.
%
% Inputs:
%   params - virmen log data (output of logParams.m)
%   envParams - environment parameters (reward zone, false alarm zone, and
%       track length)
%   useLaps - range of laps to be analzyed
%   rewCut - reduction values for reward window size
%

% set bin width
binWidth = 5;

% read environment parameters
rewStart = envParams.rewStart;
rewStop = envParams.rewStop;
rewLen = rewStop-rewStart;
trackLength = envParams.trackLength;
yBins = 0:binWidth:trackLength;

% read use laps
if isempty(useLaps) || isequal(useLaps,0)
    startLap = 1;
    stopLap = 0;
else
    startLap = useLaps(1);
    stopLap = useLaps(end);
end

% set analysis parameters
speedThresh = 1;    % speed threshold
timeThresh = 60;   % stop threshold for smoothing (sec*60)
nCut = length(rewCut);       % number of reward window sizes


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


%% Identify stops
% An additional stop is added and then removed from each lap to prevent
% errors caused by laps with no stops.

% velocity with time threshold averaging
velTh = cellfun(@(a) movmean(abs(a),[timeThresh-1 0]),vel,...
    'UniformOutput',0);
velTh = cellfun(@(a) [(speedThresh+1)*ones(60,1);a(61:end)],...
    velTh,'UniformOutput',0);

% identify the start index of all stops
idxStopAll = cellfun(@(a) contiguous(a<=speedThresh,1),velTh,'UniformOutput',0);
idxStopEnds = cellfun(@(a) [[1 1]; a{1,2}],idxStopAll,'UniformOutput',0);

% identfy start index of all stops
idxStopStart = cellfun(@(a) a(:,1),idxStopEnds,'UniformOutput',0);

% identify the start position of all stops
yStop = cellfun(@(a,b) b(a),idxStopStart,yLaps,'UniformOutput',0);
yStop = cellfun(@(a) a(2:end),yStop,'UniformOutput',0);
yStopCat = vertcat(yStop{:});

% calculate stop distribution
if ~isempty(yStopCat)
    stopDensity = ksdensity(yStopCat,yBins);
    stopDist = stopDensity*100;
else
    stopDist = zeros(size(yBins));
end

% calculate y position of reward stops only
yStopRew = cellfun(@(a) a(find(inrange(a,[rewStart rewStop]),1)),yStop,'UniformOutput',0);

% idenitfy all indices in true stops
idxStopTrueAll = cellfun(@(a) arrayfun(@(b,c) b:c,a(:,1),a(:,2),...
    'UniformOutput',0),idxStopEnds,'UniformOutput',0);

% idenitfy all y positions in true stops  (+1 tries to account for index shifting by diff)
yStopAll = cellfun(@(a,b) cellfun(@(c) b(c+1),a,...
    'UniformOutput',0),idxStopTrueAll,yLaps,'UniformOutput',0);
yStopAll = cellfun(@(a) a(2:end),yStopAll,'UniformOutput',0);

perSuccess = zeros(nCut,1);
for cc = 1:nCut
    % get current cut length
    if rewLen-rewCut(cc)<0
        curCut = rewLen;
    else
        curCut = rewCut(cc);
    end

    % calculate number of stops in reward zone
    numStopsInReward = cellfun(@(a) sum(cellfun(@(b) any(inrange(b,...
        [rewStart rewStop-curCut])),a)),yStopAll);

    % identify all success and false alarm laps
    success = numStopsInReward>0;

    % calculate percent success
    perSuccess(cc) = mean(success)*100;
end

% results struct
stopStats = struct();
stopStats.stopDist = stopDist;
stopStats.yStopRew = yStopRew;
stopStats.perSuccess = perSuccess;

end

