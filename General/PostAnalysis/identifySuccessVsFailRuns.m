%% identifySuccessVsFailRuns

% reclear all variables
clear; close all; clc

% reset directory
cd('D:\AD_Project\imagingData\data')
p1 = pwd;

folds = dir('*\*\loc*');

% reward threshold
rewardThresh=0.8;


%% Find lap indices

for ii = 1:length(folds)
    curFold = [folds(ii).folder '\' folds(ii).name];
    cd(curFold)
    
    if ~isfile('m.mat') || isfile('rewardInfo.mat')
        continue
    end

    disp(ii)
    % load m.mat
    load('m.mat')
    
    % get lap indices
    [runIdx,lapIdx,rewardIdx] = getLapIdx(m,rewardThresh);

    save('rewardInfo.mat','runIdx','rewardIdx','lapIdx')
end

cd(p1)
