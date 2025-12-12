%% prepCurveFitting.m
% Prepares behavior data for output to prism. Used with curve fitting
% analysis. Outputs both individual and groups behavior for the selected
% data field.
%

clear; close all; clc


%% Set inputs

% name of data matrix (allMiceDataB3 is current as of 12/20/2023)
dataMatrix = 'data\allMiceDataB4.mat';

load('groupIDs.mat')
% define genotypes (current as of 12/20/2023)
% groups = {[1 2 3 5 6 7 12],[4 8 9 10 11 13 14]}; % {[WT],[AD]}

% field to analyze
fieldName = 'perSuccess';
% fieldName = 'perSuccessTry';

useSex = 3;

%% Extract individual and group data

% load data
load(dataMatrix)

% maximum number of sessions
maxSessions = max(cellfun(@length,{B(:).(fieldName)}));

% initialize out arrays
curData = zeros(length(B),maxSessions);
data = cell(1,2);

% extract data for current field
for ii = 1:length(B)
    curD = B(ii).(fieldName);
    curD(end+1:maxSessions) = nan;
    curData(ii,:) = curD;
end

% sort data by group
for g = 1:2
    data{g} = curData(intersect(sexes{useSex},groups{g}),:)';
end

% output individual data
dataInd = cat(2,data{:});

% output groups means
dataMean = cellfun(@(x) nanmean(x,2),data,'UniformOutput',0);
dataMean = cat(2,dataMean{:});

