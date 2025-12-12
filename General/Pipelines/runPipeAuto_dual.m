%% runPipeAuto
%

clear all; clc

warning('off')

% Establish inputs
baseFolder = 'D:\VirusMice\TM240611-V1\';

% Select mouse raw data folder
cd(baseFolder)

top = pwd;

% Find all days
d = dir();
isub = [d(:).isdir]; %# returns logical vector
nameDays = {d(isub).name}';
nameDays(ismember(nameDays,{'.','..'})) = [];

% Find all locs
allLocs = {};

for ff = 1:length(nameDays)
    cd(nameDays{ff})
    
    d2 = dir('loc*');
    isub = [d2(:).isdir]; %# returns logical vector
    subLocs = {d2(isub).name}';
    
    for gg = 1:length(subLocs)
        allLocs{end+1} = [nameDays{ff} '\' subLocs{gg}];
    end
    
    cd(top)
end


%% Run pipeline

for ff = 1:length(allLocs)
    cd([baseFolder allLocs{ff}])
    
    % check claim
    pause(rand*3+1)
    if exist('claimed.mat','file')~=0
        cd(top)
        disp(0)
        continue
    else
        save('claimed.mat','baseFolder')
        disp(1)
    end
    

    %% run pipeline
    
    run('runPipeline_dual.m')
    close all
    
    %% Reestablish Program
    
    % Establish inputs
    load('claimed.mat')
    
    % Select mouse raw data folder
    cd(baseFolder)
    
    top = pwd;
    
    % Find all days
    d = dir();
    isub = [d(:).isdir]; %# returns logical vector
    nameDays = {d(isub).name}';
    nameDays(ismember(nameDays,{'.','..'})) = [];
    
    % Find all locs
    allLocs = {};
    
    for qq = 1:length(nameDays)
        cd(nameDays{qq})
        
        d2 = dir('loc*');
        isub = [d2(:).isdir]; %# returns logical vector
        subLocs = {d2(isub).name}';
        
        for gg = 1:length(subLocs)
            allLocs{end+1} = [nameDays{qq} '\' subLocs{gg}];
        end
        
        cd(top)
    end
    %%
    
    cd(top)
end


cd(top)









