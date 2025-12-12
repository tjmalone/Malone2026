%% trainingLaps
% Generates a time course of training data for an individual mouse over
% multiple days. This is used to manually generate a matrices containing
% the number of laps, number of rewards, and environment types for all
% mice combined.
%

%% Set paramters

clear; close all; clc

% load all text files in current folder and in all sub folders
d = dir('**/*T*.txt');

% sort d by date
[~,idxSort] = sort([d.datenum]);
d = d(idxSort);

% remove duplicate files
[~,idxDup] = unique({d(:).name});
d = d(idxDup);

numDays = size(d,1);            % number of virmen logs
numRuns = zeros(numDays,1);     % number of runs per day
numRewards = zeros(numDays,1);  % number of rewards per day
tLength = zeros(numDays,1);     % track length (for env validation)
mDates = cell(numDays,1);       % dates (for validation)


%% Process virmen logs

% cycle through all virmen files
for f = 1:numDays
    % load virmen log
    logData = readLog([d(f).folder '\' d(f).name], 9);
    params = logParams(logData);
    
    % get reward number
    numRewards(f) = length(params.rewardIdx);
    
    % calculate number of runs
    [~,numRuns(f)] = identifyLaps(params.y);
    
    % calculate date
    mFile = d(f).name;
    mDates{f} = [num2str(mFile([5 6])) '/' num2str(mFile([7 8]))...
        '/' num2str(mFile([3 4]))];
    
    % calculate track length
    tLength(f) = round(max(params.y));
    
end

% set nans to 0
numRuns(isnan(numRuns)) = 0;


%% Save raw training data

trData = struct();
trData.numRuns = numRuns;
trData.numRewards = numRewards;
trData.tLength = tLength;
trData.mDates = mDates;

save('expLaps.mat','trData')


%% Filter duplicate days

clear

switchDay = '10/02/23';

load('expLaps.mat','trData')
numRuns = trData.numRuns;
numRewards = trData.numRewards;
mDates = trData.mDates;
tLength = trData.tLength;

ii = 1;

while ii<length(numRuns)
    if ~strcmp(mDates{ii},switchDay) && strcmp(mDates{ii},mDates{ii+1}) && tLength(ii)==tLength(ii+1)
        numRuns(ii) = sum(numRuns(ii:ii+1));
        numRuns(ii+1) = [];
        numRewards(ii) = sum(numRewards(ii:ii+1));
        numRewards(ii+1) = [];
        mDates(ii+1) = [];
        tLength(ii+1) = [];
    end
    ii = ii+1;
end


%% Save filtered training data

trData = struct();
trData.numRuns = numRuns;
trData.numRewards = numRewards;
trData.tLength = tLength;
trData.mDates = mDates;

save('expLapsT.mat','trData')

