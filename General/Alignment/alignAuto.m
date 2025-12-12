%% alignAuto.m
% Peforms alignments for analysis on a full set of days for one mouse.
%

clear all; close all; clc

warning('off')


%% Identify analysis folders

p = pwd;

saveFolder = 'Alignments\';
dataFolder = [p '\Mice'];

% find all pre-defined folder sets within mouse data
folds = findSubF('folders5.mat',1,dataFolder,0);

%define analysis type/number correspondence
prior = 1;
same = 2;
% txtType = {'prior','same'};
% txtType = {'env3','ltm'};
% txtType = {'ltmS','same'};
% txtType = {'mix','same'};
txtType = {'man','NaN'};


%% Run alignment for all comparison sets

for i = 1:length(folds)
    
    % load current comparison set
    load(folds{i})
    
    % run alignment for prior/same analysis types
    for j = 1:2 % prior:same
        if isempty(folders(j).use) || folders(j).use==0
            continue
        end
        
        nLoc = length(folders(j).env1);
        
        % run alignment for all locs
        for k = 1:nLoc
            
            % extract mouse ID and loc info
            curID = regexp(folders(j).env1{k},'ID.*?\\','match','once');
            curLoc = regexp(folders(j).env1{k},'loc.*?\\','match','once');
            
            % generate alignment save folder
            curSF = [saveFolder curID curLoc txtType{j}];
            
            % move to save folder
            mkdir(curSF)
            cd(curSF)
            
            % run alignment
            copyfile('D:\AnalysisCode\Alignment\alignmentFF.m','alignmentFF.m')
            alignmentFF(folders(j).env1{k},folders(j).env2{k})
            
            %%%%%%
            %             novelCells = findNovel();
            %             save('novelCells.mat','novelCells')
            %%%%%%
            
            % return to base folder
            cd(p)
            
        end
    end
end
