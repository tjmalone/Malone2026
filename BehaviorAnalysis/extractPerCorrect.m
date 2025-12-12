function behavior = extractPerCorrect(params,envParams,useLaps)
%% extractActiveBehavior
% A light version of extractActiveBehavior that only calculates percent
% success for a given set of criteria Called in behaviorAD_compareThresh
%
% Inputs:
%   params - virmen log data (output of logParams.m)
%   envParams - environment parameters (reward zone, false alarm zone, and
%       track length, threshold criteria)
%   useLaps - range of laps to be analzyed
%

% read environment parameters
rewStart = envParams.rewStart;
rewStop = envParams.rewStop;

% read use laps
if isempty(useLaps) || isequal(useLaps,0)
    startLap = 1;
    stopLap = 0;
else
    startLap = useLaps(1);
    stopLap = useLaps(end);
end

% set analysis parameters
speedThresh = envParams.speedThresh;
timeThresh = envParams.timeThresh;
winBuff = envParams.winBuff;


%% Calculate rolling speed

% extract behavior data
t = params.t;
y = params.y;

% calculate number of runs
[runIdx,nRuns] = identifyLaps(y,startLap,stopLap);

% segment data by lap
yLaps = cell(nRuns,1);
tLaps = cell(nRuns,1);
for m = 1:nRuns
    yLaps{m} = y(runIdx(m,1):runIdx(m,2));
    tLaps{m} = t(runIdx(m,1):runIdx(m,2));
end

% calculate velocity
vel = cellfun(@(a,b) diff(a)./diff(b),yLaps,tLaps,'UniformOutput',0);


%% Calculate stop statistics with all thresholds

behavior = zeros(length(timeThresh),1);

% stop versus try
for ii = 1:length(timeThresh)
    %% Identify stops
    % An additional stop is added and then removed from each lap to prevent
    % errors caused by laps with no stops.

    % velocity with time threshold averaging
    velTh = cellfun(@(a) movmean(abs(a),[timeThresh(ii)-1 0]),vel,...
        'UniformOutput',0);
    velTh = cellfun(@(a) [(speedThresh(ii)+1)*ones(60,1);a(61:end)],...
        velTh,'UniformOutput',0);

    % idenitify the start index of all stops
    idxStopAll = cellfun(@(a) contiguous(a<=speedThresh(ii),1),velTh,'UniformOutput',0);
    idxStopEnds = cellfun(@(a) [[1 1]; a{1,2}],idxStopAll,'UniformOutput',0);

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

    %% Caclulate statistics on all laps

    % identify all success and false alarm laps
    success = numStopsInReward>0;

    % calculate percent success
    perSuccess = mean(success)*100;


    %% Save results

    behavior(ii) = perSuccess;

end

end

