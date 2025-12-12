%% Generate recall alignments

clear; close all; clc

useFolders = [1 12];
saveFolder = 'results_recall';


%% Generate results.mat
% Run once to generate results.mat. Takes a while to run, so is commented
% to avoid repeating.

alignmentPost(useFolders,saveFolder)


%% Generate alignment array

% Run to generate alignments. For recall alignments, all alignments from
% the original pairwise alignment are checkewd without manual or automated
% checking.

alignData = uniqueAligns(length(useFolders),saveFolder);

alignsFinal = saveAligns(alignData,saveFolder);

