%% training_timecourse
% generate a time course of training data for an individual mouse over
% multiple days


%% Set paramters

d = dir('*T*.txt');

% d = d([1:11 13 16:19]);
% d = d(12:20);
% d = d([12 14:15]);

% 1202
% d = d(14:end);

% 1207
% d = d(9:end);

% 1209
d = d(23:end);

numDays = size(d,1);            % number of training days
numRuns = zeros(numDays,1);     % number of runs per day
timeRuns = zeros(numDays,1);    %
predLick = zeros(numDays,1);
perSlow = zeros(numDays,1);

rwLoc = 360;

% rwLoc = [35,535];
% rwLoc = 890;

distBeforeReward = 10;
distBeforeSlowDown = 190;

% 1m track
% rwLoc = 70;
% distBeforeSlowDown = 40;


dst = [10 100];     % distances for calculation


%%

for f = 1:numDays
    logData = readLog(d(f).name, 9);
    params = logParams(logData);
    
    % record length of training session
    t = range(logData(:,1))*24*60;
    timeRuns(f) = t;
    
    % record number of runs completed
    numRuns(f) = size(params.yReward,1)/length(rwLoc);

    % calculate predictive licking
    licking = predLicking(params.y,params.lickIdx,params.rewardIdx,dst);
    
    % record predictive licking
    predLick(f) = licking.perPredLick;
    
    loc1 = rwLoc(end)-dst(1);
    loc2 = loc1-dst(2);
    
    speedCompare = slowingDown(params.y,params.rewardIdx,loc1,loc2);
    
    perSlow(f) = speedCompare.perSlow;
    
    % plot example data
    if f==numDays
        predLicking(params.y,params.lickIdx,params.rewardIdx,dst,1);
        slowingDown(params.y,params.rewardIdx,loc1,loc2,1);
    end
    
end

numRuns(isnan(numRuns)) = 0;
timeRuns(isnan(timeRuns)) = 0;
predLick(isnan(predLick)) = 0;
perSlow(isnan(perSlow)) = 0;


%%

figure

subplot(1,2,1); hold on

xidx = 1:numDays;

plot(xidx,predLick*100,'g')
plot(xidx,perSlow*100,'b');
set(gca,'fontsize', 18);

hold off

ylim([0 100]);
xlim([0.9 xidx(end)+0.1]);

leg1 = legend('Predictive licking', 'Percentage slow',...
    'Orientation', 'horizontal', 'Location', 'northoutside');

xlabel('day');
ylabel('percent of runs')

subplot(1,2,2); hold on

plot(xidx,numRuns./timeRuns*60, 'r');
set(gca,'fontsize', 18);
hold off

xlabel('day');
ylabel('runs per hour')

xlim([0.9 xidx(end)+0.1]);

% savefig(figure(1), 'trainingUpdate.fig');

