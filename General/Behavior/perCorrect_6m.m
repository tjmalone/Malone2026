%% training_timecourse
% generate a time course of training data for an individual mouse over
% multiple days

clear; close all; clc

%% Set parameters

d = dir('*T*.txt');
%
numDays = size(d,1);            % number of training days

% laps to include
startLap = 2;
stopLap = 20;

perCorrect = zeros(numDays,1);
mDates = cell(numDays,1);


%%

for f = numDays-10:numDays
    logData = readLog(d(f).name, 9);
    params = logParams(logData);
    
    %% Calculate rolling speed
    
    y = params.y;
    
    % calculate number of runs
    [runIdx,nRuns] = identifyLaps(y,startLap,stopLap);
    
    rewardIdxCur = inrange(params.rewardIdx,[runIdx(1,1) runIdx(end,2)]);
    numRewards = sum(rewardIdxCur);

    perCorrect(f) = numRewards/nRuns*100;
    
    % calculate date
    mFile = d(f).name;
    mDates{f} = [num2str(mFile([5 6])) '/' num2str(mFile([7 8]))...
        '/' num2str(mFile([3 4]))];
    
end

perCorrect(isnan(perCorrect)) = 0;


%%

data = struct();
data.perCorrect = perCorrect;
data.mDates = mDates;
