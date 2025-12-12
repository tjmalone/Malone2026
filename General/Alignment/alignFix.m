%% alignFix.m
%
% Cycles through FOV alignments to check for proper alignment. Allows user
% to manually align poorly aligned FOV pairs.
%

clear; close all; clc

start = [1,2];
cats = 1:2;


%%

clearvars -except cats start

p = pwd;

load('folders')
n = length(folders);

checkFile = 'Figures\Stage 2 - pre vs post alignment.fig';

% checkFiles = findSubF(checkName,3,p,0);

for ii = start(1):n-1
    for jj = start(2):n
        start = [ii,jj];
        
        curFolder = [num2str(ii) '\' num2str(ii) '_' num2str(jj)];
        cd(curFolder)
        
        F = openfig(fullfile(p,curFolder,checkFile));
        
        keep = input('Good alignment?: ');
        close(F)
        
        
        if keep==1
            cd(p)
            continue
        else
            return
        end
    end
    start(2) = ii+2;
end

cd(p)


%%

cd(p)

load allFile_names.mat

f = start(1);
k = start(2);

nFold = length(allFile_names);

clear results

% load reference day
load(allFile_names{f});
NRefCells=size(data,1);         % reference roi number

% initialize results struct
results.refDay = f; % reference day number
results.allFOVs = {}; % combined assignment data with other days

% [reference index, common cell index for other days]
results.allRefOthers(:,1) = 1:1:NRefCells;
results.perctCommonCells = []; % percentage of overlapping cells

% load current non-ref data
load(allFile_names{k});
NCells = size(data,1);

% make directory for current alignment (ref day - comparison day)
foldername = sprintf('%s_%s',num2str(f),num2str(k));
useIdx = [f k]; % indices of comparison days

% set inputs for alignment code
useFile_names = allFile_names(useIdx);
results_directory = fullfile(pwd,num2str(f),foldername);
save([results_directory '/useFile_names.mat'],'useFile_names');

% run alignment and save data
clear cell_registered_struct;

copyF = 'commonCellIdentificationManual.m';
if ~isfile(copyF)
    srcF = which(copyF);
    copyfile(srcF);
end

d=dir([results_directory '\cellRegistered*']);
if ~isempty(d)
    for ii = 1:length(d)
        delete([results_directory '\' d(ii).name])
    end
end
                
run(fullfile(pwd,'commonCellIdentificationManual.m'))

close all

