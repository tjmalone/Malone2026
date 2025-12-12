%% alignAuto.m
% Peforms alignments for analysis on a full set of days for one mouse.
%

clear all; close all; clc

warning('off')


%% Identify analysis folders

dataFolder = 'D:\AD_Project\imagingData';
cd(dataFolder)

saveFolder = 'Alignment11';
fileType = 'plus';

% find all pre-defined folder sets within mouse data
folds = findSubF('imageFiles11.mat',2,[],0);


%% Run alignment for all comparison sets

for ii = 1:length(folds)
    cd(erase(folds{ii},'\imageFiles11.mat'));
    
    mkdir(saveFolder)
    cd(saveFolder)
    
    % check claim
    if exist('claimedA2.mat','file')~=0
        disp(0)
        pause(0.1)
        continue
    else
        claimed = 1;
        save('claimedA2.mat','claimed')
        disp(1)
    end
    
    % load current comparison set
    load(folds{ii})
    folders = files.(fileType);
    save('folders.mat','folders')
    
    % run alignment
    copyfile('D:\AD_Project\imagingData\alignmentByIdx.m','alignmentByIdx.m')
    try
        alignmentByIdx()
    catch
        delete('claimedA2.mat');
    end
    
end

cd(dataFolder)
