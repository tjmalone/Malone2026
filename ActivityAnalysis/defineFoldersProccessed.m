%% defineFoldersProccessed
% Defines suite2p folders based on alignment folders for learning days and
% saves as nested cell array (FOV{day}). Also saves common cell alignment,
% and defines which session belongs to which group (WT vs. AD) or sex
% (female vs. male) for imaging analyses.

clear; clc; close all;

p = 'D:\AD_Project\imagingData\data';
cd(p)

% set replace paths
repPaths = {'D:\NewMice','D:\Imaging_data\AD_project\ADBatch6'};


%% Generate folder array

% find all alignments of desired type
resultType = 'results_learning';
foldMice = findSubF(resultType,3,[],0);
nMice = length(foldMice);

% initialize output array
foldersLearning = cell(nMice,1);
alignsLearning = cell(nMice,1);

% cycle through FOV
for ii = 1:nMice
    % change to alignment directory
    cd(foldMice{ii})
    
    % load alignment and folder indices
    load('alignsFinal.mat')
    load([resultType '.mat'])
    useSess = results.useIdx;
    
    % load data folders
    cd ..
    load('folders.mat')
    save('folders.mat','folders')

    % check folder validity
    try
        cd(folders{1})
    catch
        disp(ii)
        for jj = 1:length(repPaths)
            for kk = 1:length(folders)
                folders{kk} = replace(folders{kk},repPaths{jj},p);
            end
        end

        save('folders.mat','folders')
    end

    % save alignments and FOV folders
    foldersLearning{ii} = folders(useSess);
    alignsLearning{ii} = alignsFinal;
end

cd(p)

% save output array
save('foldersLearning.mat','foldersLearning','alignsLearning');


%% Define genotype groups

groupIDs = {'WT','AD'};
sexIDs = {'all','female','male'};

% Batch 2-4
% groups = {[1:6 9:14 23:24],[7:8 15:22 25:28]};

% Batch 2-4 + some Batch 5
% groups = {[1:6 9:14 23:24 33:36],[7:8 15:22 25:32]};

% Batch 2-5
groups = {[1:6 9:14 23:24 31:32 37:38 41:42],[7:8 15:22 25:28 29:30 33:36 39:40]};
sexes = {1:42,[1:2 5:8 11:12 19:20 27:28 31:34 37:38 41:42],[3:4 9:10 13:18 21:26 29:30 35:36 39:40]};

if sum(cellfun(@length,groups)) ~= nMice
    error('Group mismatch')
end

if sum(cellfun(@length,sexes(2:3))) ~= nMice
    error('Sex mismatch')
end

save('groupIDs.mat','groupIDs','groups','sexIDs','sexes')
