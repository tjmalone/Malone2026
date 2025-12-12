%% runAD_timecourse_Leica

clear; close all; clc

run('quantTauLayersKC_AllTimecourses.m')
% Extract tau intensity
% Determine intensity threshold for each batch
    % For each of the young AD mice (2 and 3 months) in a batch, find the 99th percentile value intensity value in layer 2
    % Average this 99th percentile value for the young AD mice in each batch
    % makes 3 folders
        % backgroundSub: subtracting threshold from raw intensity
        % perArea: finds percent area above the threshold
        % raw: no thresholding
    % main output: tauDataGrouped.mat
        % has the tau intensity for each mouse in layers 1, 2 and 3
% Makes plots
% output folder: Z:\SNM\labMembers\KC\AD_Project\Histology\AD_timecourse_Leica\FinalVersion

run('plotTauTimecourseForFigs.m')
%Make final plots
    % by sex
    % young vs old

