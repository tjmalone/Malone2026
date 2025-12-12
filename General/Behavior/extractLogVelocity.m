%% extractLogVelocity
% generate a time course of training data for an individual mouse over
% multiple days

clear; close all; clc


%% Set parameters

% filesRun = {'20240829T125749_V1-fam.txt'};

% define files
d = dir('*T*.txt');
filesRun = {d(:).name};
numDays = length(filesRun);            % number of training days

speedThresh = 1;


%% Process all logs

velMean = zeros(1,numDays);

for ii = 1:numDays
    %% Calculate general behavior statistics
    logData = readLog(filesRun{ii}, 9);
    params = logParams(logData);

    % record number of runs completed
    [runIdx,numRuns] = identifyLaps(params.y);

    % get all use indices
    velAll = [];
    for jj = 1:numRuns
        % get current indices
        curIdx = runIdx(jj,1):runIdx(jj,2);

        % calculate velocity
        curY = params.y(curIdx);
        curT = params.t(curIdx);
        curVel = diff(curY)./diff(curT);

        % apply speed thresholding
        curVelMove = curVel(curVel>speedThresh);
        % curVelMove = curVel;

        % store velocities
        velAll = [velAll; curVelMove];
    end

    velMean(ii) = mean(velAll,'omitnan');

end


%%

figure; hold on

plot(velMean)

title('Mean velovity (while moving)')
xlabel('Session')
ylabel('Velocity (cm/s)')

savefig('velocity.fig')
