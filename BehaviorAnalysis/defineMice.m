%% defineMice.m
% Saves indices defining which mice belong to which group (WT vs. AD) or
% sex (female vs. male) for behavior analyses
%

clear; clc; close all;

p = 'D:\AD_Project\Behavior';
cd(p)


%% Define genotype groups

groupIDs = {'WT','AD'};
sexIDs = {'all','female','male'};

% Manually defined groups for AD imaging batch 2-5 (all imaging mice)
groups = {[1:3 5:7 12 16 19 21],[4 8:11 13:15 17:18 20]};
sexes = {1:21,[1 3:4 6 10 14 16:17 19 21],[2 5 7:9 11:13 15 18 20]};

% if sum(cellfun(@length,groups)) ~= nMice
%     error('Group mismatch')
% end
% 
% if sum(cellfun(@length,sexes)) ~= nMice
%     error('Sex mismatch')
% end

save('groupIDs.mat','groupIDs','groups','sexIDs','sexes')
