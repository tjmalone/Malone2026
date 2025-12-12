%% appendTrueDays
% generate a true days variable to foldersLearning.m to mark the true days
% that the learning days corespond to (accounts for skipped activity days)

clear; close all; clc

% load data
load('foldersLearning.mat')

nDays = 11;
nFOV = length(foldersLearning);
trueDays = zeros(nFOV,nDays);

for ii = 1:nFOV

    idx = 1;
    for jj = 1:nDays
        % account for skip days
        if jj==4 && ismember(ii,[11 12])
            % all FOV with skip day 3
            idx = idx+1;
        elseif jj==7 && ismember(ii,13:18)
            % all FOV with skip day 6
            idx = idx+1;
        end

        % save true day index
        trueDays(ii,jj) = idx;
        idx = idx+1;
    end
end

% re-save data
save('foldersLearning.mat','alignsLearning','foldersLearning','trueDays')