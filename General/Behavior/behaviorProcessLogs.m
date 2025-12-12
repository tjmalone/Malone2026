function behaviorProcessLogs(startNew,nRuns,lickType)
%% behaviorProcessLogs
% generate a time course of training data for an individual mouse over
% multiple days


%% Set parameters

% whether to reuse earlier data
if nargin<1 || isempty(startNew)
    startNew = 1;
end

% number of runs to use
if nargin<2 || isempty(nRuns)
    nRuns = 0;
end

% set lick calculation type prediction among all licks ('all') or percent
% of predictive licking laps ('raw')
if nargin<3 || isempty(lickType)
    lickType = 'all';
end

% define files
d = dir('*T*.txt');
numDays = size(d,1);            % number of training days


%% Intialize/load data matrix

if startNew
    numRuns = zeros(numDays,1);     % number of runs per day
    perCorrect = zeros(numDays,1);
    timeRuns = zeros(numDays,1);    %
    runsPerHour = zeros(numDays,1);
    tLength = zeros(numDays,1);
    pc10 = zeros(numDays,1);
    pLicking = zeros(numDays,1);
    pSlowing = zeros(numDays,1);

    mDates = cell(length(d),1);

    stDay = 1;

elseif ~startNew
    load('trData.mat');
    numRuns = trData.numRuns;
    perCorrect = trData.perCorrect;
    timeRuns = trData.timeRuns;
    runsPerHour = trData.runsPerHour;
    tLength = trData.tLength;
    pc10 = trData.pc10;
    pLicking = trData.pLicking;
    pSlowing = trData.pSlowing;
    mDates = trData.mDates;

    stDay = length(numRuns)+1;
end


%% Process all logs

for f = stDay:numDays
    %% Calculate general behavior statistics
    logData = readLog(d(f).name, 9);
    params = logParams(logData);

    % record number of runs completed
    [runIdx,numRuns(f)] = identifyLaps(params.y,[],nRuns);

    lastLapIdx = runIdx(end,2);
    numRewards = sum(params.rewardIdx<=lastLapIdx);
    perCorrect(f) = numRewards/numRuns(f);

    % record length of training session
    time = (logData(:,1)-logData(1,1))*24*60;
    timeRuns(f) = time(lastLapIdx);
    runsPerHour(f) = numRuns(f)/timeRuns(f)*60;

    difY10 = runIdx(1:min(10,size(runIdx,1)),2);
    numRewards10 = find(params.rewardIdx<=difY10(end),1,'last');
    if ~isempty(numRewards10)
        pc10(f) = numRewards10/length(difY10);
    end

    % calculate date
    mFile = d(f).name;
    mDates{f} = [num2str(mFile([5 6])) '/' num2str(mFile([7 8]))...
        '/' num2str(mFile([3 4]))];

    % calculate track length
    tLength(f) = round(max(params.y));
    if tLength(f)<105
        pSlowing(f) = NaN;
        pLicking(f) = NaN;
        continue
    end


    %% Calculate predict slowing/licking

    % track length
    trEnd = tLength(f);

    % set variable calculation parameters
    if trEnd==600
        binSlow = 90;
    elseif trEnd>600
        binSlow = 150;
    elseif trEnd<600
        binSlow = 50;
    end

    % set constant parameters
    trStart = 0;
    postRewSlow = 30;
    preRewLick=20;
    postRewLick = 30;

    % calculate predictive slowing
    try
        curSlow = speedChange_Percentile_TM(d(f).name,binSlow,postRewSlow,...
            nRuns,trStart,trEnd,0);
        pSlowing(f) = curSlow.meanCentile/100;

    catch
        pSlowing(f) = NaN;
    end

    % calculate predictive licking
    try
        try
            curLick = predLicking_allRuns_allLicking(params.y,params.lickIdx,...
                params.rewardIdx,preRewLick,postRewLick,nRuns,0);
        catch
            curLick = predLickingTwo_allRuns_allLicking(params.y,...
                params.lickIdx,params.rewardIdx,preRewLick,postRewLick,nRuns,0);
        end
        if strcmp(lickType,'all')
            pLicking(f) = curLick.perPredInAll;   % norm predictive licking
        elseif strcmp(lickType,'raw')
            pLicking(f) = curLick.perPredLick;  % raw perdicitve licking
        end
    catch
        pLicking(f) = NaN;
    end

end

% convert nans to zeros
numRuns(isnan(numRuns)) = 0;
perCorrect(isnan(perCorrect)) = 0;
timeRuns(isnan(timeRuns)) = 0;
pLicking(isnan(pLicking)) = 0;
pSlowing(isnan(pSlowing)) = 0;


%% save full data structure

trData = struct();
trData.numRuns = numRuns;
trData.perCorrect = perCorrect;
trData.timeRuns = timeRuns;
trData.runsPerHour = runsPerHour;
trData.tLength = tLength;
trData.pc10 = pc10;
trData.pLicking = pLicking;
trData.pSlowing = pSlowing;
trData.mDates = mDates;

save('trDataByRun.mat','trData')


%% Plot behavior

close all
figure

subplot(1,2,1); hold on

xidx = 1:numDays;

plot(xidx,trData.pSlowing*100,'LineWidth',2,'Color',[0,0.9,0]);
plot(xidx,trData.pLicking*100,'LineWidth',2,'Color',[0,0.4,0])
plot(xidx,trData.perCorrect*100,'LineWidth',2,'Color',[0.7,0,0]);

set(gca,'fontsize', 18);

hold off

ylim([0 100]);
xlim([0.9 xidx(end)+0.1]);

leg1 = legend('%slow','%lick','%correct',...
    'Orientation', 'horizontal', 'Location', 'northoutside');

xlabel('day');
ylabel('percent of runs')

subplot(1,2,2); hold on

plot(xidx,trData.runsPerHour, 'r');
set(gca,'fontsize', 18);
hold off

xlabel('day');
ylabel('runs per hour')

xlim([0.9 xidx(end)+0.1]);
% ylim([0 100])

savefig(figure(1), 'trainingByRun.fig');
saveas(figure(1), 'trainingByRun.tif');


end

