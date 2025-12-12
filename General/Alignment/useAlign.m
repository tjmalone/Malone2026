%% Generate results.mat

% Run once to generate results.mat. Takes a while to run, so is commented
% to avoid repeating.

alignmentPost(1:24,'results_all')


%% Generate alignment array

% Run to generate preliminary alignments. This output matrix can be copied
% into excel to save/process
clear; clc
alignData1 = uniqueAligns(1:24,'results_all');
alignData2 = peakAligns(alignData1);

%%

% Run to compare all potential alignments. Manually input whether each
% alignment is good or bad. useFOV can be copied to excel to process
[useFOV,projection_sum,SFC] = alignmentPostCheck([1 12],num2cell(alignData(1:end,2:end-1)),1,1);

%%

% Run to compare specific alignments using same reference frame
useFOV = alignmentPostCheck(1:5,num2cell(alignData([167 193 218 309 310 377],2:end-1)),3,0,projection_sum,SFC);

%%

% Run to compare specific alignments using different reference frame
useFOV = alignmentPostCheck([2 9 16 18 20]+3,num2cell(alignData([167 193 218 309 310 377],2:end-1)),2+3);
% useFOV = alignmentPostCheck([2 9 16 18 20]+1,num2cell(alignData([401:405],2:end-1)),3);


%%

alignsFinal = saveAligns(alignData,'results_recall');

