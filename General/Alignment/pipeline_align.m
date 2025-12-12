%% Automatically align cells
% Run alignAutoSuite2p (which calls alignmentSuite2p) or run
% alignmentSuite2p directly. directories for alignment must be pre-defined
% in .mat files (one file per FOV). The file structure can be a cell array
% of drectories or a struct. If a struct is used, a field of interest
% (fileType) must be defined, and this field must be a cell array of
% directories.


%% Check automated alignments
% Run alignCheck to manually evaluate the quality of all automated
% alignments.


%% Manually align cells
% Use runManAlignment to manually align session pairs as necessary. Run the
% first to generate .tifs to align. Input alignment information based on
% alignment in illustrator and run the second section. If proposed manually
% alignment is correct, run the final section to manually align.


%% Initialize input variables

clear; clc

% useFolders = 1:11;
% saveFolder = 'results_learning';

useFolders = [1 12];
saveFolder = 'results_recall';


%% Generate results.mat
% Run once to generate results.mat. Takes a while to run, so is commented
% to avoid repeating.

alignmentPost(useFolders,saveFolder)


%% Generate alignment array

% Run to generate preliminary alignments. This output matrix can be copied
% into excel to save/process

% alignData = uniqueAligns(useFolders,saveFolder);
alignData = uniqueAligns(length(useFolders),saveFolder);


%% Automatically check alignments

alignData = peakAligns(alignData,saveFolder);


%% Manually check alignments

% Run to compare all potential alignments. Manually input whether each
% alignment is good or bad. useFOV can be copied to excel to process
% [useFOV,projection_sum,SFC] = alignmentPostCheck([1 12],num2cell(alignData(1:end,2:end-1)),1,1);
[useFOV,projection_sum,SFC] = alignmentPostCheck(useFolders,num2cell(alignData(1:end,2:end-1)),6,1);


%%

% Run to compare specific alignments using same reference frame
% useFOV = alignmentPostCheck(1:5,num2cell(alignData([15],2:end-1)),1,0,projection_sum,SFC);
useFOV = alignmentPostCheck(useFolders,num2cell(alignData([13 16],2:end-1)),6,0,projection_sum,SFC);


%%

% Run to compare specific alignments using different reference frame
useFOV = alignmentPostCheck(useFolders,num2cell(alignData([7 10],2:end-1)),1,0);
% useFOV = alignmentPostCheck([2 9 16 18 20]+1,num2cell(alignData([401:405],2:end-1)),3);


%%

alignsFinal = saveAligns(alignData,saveFolder);

