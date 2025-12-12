% find and apply cue threshold

clear; clc; close all;

p1 = 'D:\AD_Project\imagingData';
cd(p1)

svFolder = 'data_uniThresh\cueLearning';

% load folders to analyze
load('foldersLearning.mat')

% set environment types
envGroups = {1,2:11};


%% Collect shuffles

sideName = {'Left','Right','All'};
suffName = {'_dfof','_sig'};
nSide = length(sideName);
nSuff = length(suffName);

% shuffle data
allShuffles = struct();
for gg = 1:length(envGroups)
    for ii = 1:nSide
        for jj = 1:nSuff
            allShuffles(gg).([sideName{ii} suffName{jj}]) = [];
        end
    end
end

% lags
allLags = [];

for gg = 1:length(envGroups)
    % current use sessions
    curSet = envGroups{gg};
    
    % cycle through FOV
    for n=1:length(foldersLearning)
        
        % initialize sessions folders
        dfolders = foldersLearning{n}(curSet);
        
        % cycle through sessions
        for nn=1:length(dfolders)
            % move to data directory
            disp([num2str(gg) '\' num2str(n) '-' num2str(nn)])
            cd(dfolders{nn});
            
            % load rois
            load('allROIs.mat');
            
            % cycle sides
            for ii = 1:nSide
                %cycle suffixes
                for jj = 1:nSuff
                    % load cue data
                    load(['cueAnalysis' suffName{jj} '\newScoreShuffleTemplate\'...
                        sideName{ii} '\cueCells.mat']);
                    
                    % add cue data to full data
                    A=cueCells.shuffleScores;
                    A=reshape(A,size(A,1)*size(A,2),1);
                    allShuffles(gg).([sideName{ii} suffName{jj}])(end+1:end+length(A)) = A;
                    allLags = [allLags; cueCells.realLags'];
                end
            end
        end
    end
end

cd(p1);

save([svFolder '\allLags.mat'],'allLags')
save([svFolder '\allShuffles.mat'],'allShuffles')



%% Process shuffles

% threshold values
thresholdVal = struct();

if ~isvarname('allShuffles')
    load([svFolder '\allShuffles.mat'])
end
prctileCutoff = 90;

% cycle enviornment groups
for gg = 1:length(envGroups)
    % cycle sides
    for ii = 1:nSide
        %cycle suffixes
        for jj = 1:nSuff
            % identify threshold
            thresholdVal(gg).([sideName{ii} suffName{jj}]) =...
                prctile(allShuffles(gg).([sideName{ii} suffName{jj}]),prctileCutoff);
        end
    end
end

% clear allShuffles


%% Apply thresholds

% cue cell percentage
perCue = struct();
for gg = 1:length(envGroups)
    for ii = 1:nSide
        for jj = 1:nSuff
            perCue(gg).([sideName{ii} suffName{jj}]) = [];
        end
    end
end

% cue cell percentage
scoreCue = struct();
for gg = 1:length(envGroups)
    for ii = 1:nSide
        for jj = 1:nSuff
            scoreCue(gg).([sideName{ii} suffName{jj}]) = [];
        end
    end
end

% keep fields
keepFields = 1:9;

for gg = 1:length(envGroups)
    % current use sessions
    curSet = envGroups{gg};
    
    % cycle through FOV
    for n=1:length(foldersLearning)
        
        % initialize sessions folders
        dfolders = foldersLearning{n}(curSet);
        
        % cycle through sessions
        for nn=1:length(dfolders)
            % move to data director
            disp([num2str(gg) '\' num2str(n) '-' num2str(nn)])
            cd(dfolders{nn});
            
            % load rois
            load('allROIs.mat');
            N=size(roi,3);
            
            p2 = pwd;
            
            % cycle sides
            for ii = 1:nSide
                % cycle suffixes
                for jj = 1:nSuff
                    % change directory
                    cd(['cueAnalysis' suffName{jj} '\newScoreShuffleTemplate\' sideName{ii}]);
                    
                    % current threshold
                    curThresh = thresholdVal(gg).([sideName{ii} suffName{jj}]);
                    
                    % load cue cells
                    load('cueCells.mat');
                    fNames = fieldnames(cueCells);
                    
                    % cue cell information structure
                    cueCellsUniThresh = struct;
                    
                    % copy unthresholded info
                    for kk = 1:length(keepFields)
                        cueCellsUniThresh.(fNames{keepFields(kk)}) =...
                            cueCells.(fNames{keepFields(kk)});
                    end
                       
                    % add thresholded info
                    cueCellsUniThresh.thresh = curThresh;
                    cueCellsUniThresh.cueCellInUseIdx = find(cueCells.realScores>curThresh);
                    cueCellsUniThresh.cueCellRealIdx = cueCells.useIdx(cueCells.realScores>curThresh);
                    cueCellsUniThresh.cueCellScores = cueCells.realScores(cueCells.realScores>curThresh);
                    cueCellsUniThresh.cueCellLags = cueCells.realLags(cueCells.realScores>curThresh);
                    cueCellsUniThresh.cueCellDfofAvg = cueCells.useDfofaverage(:,cueCells.realScores>curThresh);
                    
                    % save cue cell structure
                    save('cueCellsUniThresh.mat','cueCellsUniThresh');
                    
                    % update cue cell percent structure
                    perCue(gg).([sideName{ii} suffName{jj}])(end+1,1) =...
                        length(cueCellsUniThresh.cueCellRealIdx)/N;
                    
                    % update cue cell percent structure
                    scoreCue(gg).([sideName{ii} suffName{jj}])(end+1:end+length(cueCells.realScores),1) =...
                        cueCells.realScores;
                    
                    cd(p2)
                end
            end
        end
    end
end

cd(p1)

% calculate mean and SEM for cur percentages
perCueMean = struct();
perCueSEM = struct();
for gg = 1:length(envGroups)
    for ii = 1:nSide
        for jj = 1:nSuff
            perCueMean(gg).([sideName{ii} suffName{jj}]) =...
                mean(perCue(gg).([sideName{ii} suffName{jj}]),'omitnan');
            perCueSEM(gg).([sideName{ii} suffName{jj}]) =...
                nansem(perCue(gg).([sideName{ii} suffName{jj}]),1);
        end
    end
end

% save threshold information
save([svFolder '\cueInfo.mat'],'thresholdVal','perCue','perCueMean','perCueSEM','scoreCue')


%% Generate combined cue structure

clear; clc; close all;

p1 = 'D:\AD_Project\imagingData';
cd(p1)

% load alignment data
load('foldersLearning.mat')

% fields to save
keepFields  = {'cueCellInUseIdx','cueCellRealIdx','cueCellScores','cueCellLags'};

sideName = {'Left','Right','All'};
suffName = {'_dfof','_sig'};
nSide = length(sideName);
nSuff = length(suffName);

% cycle through FOV
for n=1:length(foldersLearning)
    
    % initialize sessions folders
    dfolders = foldersLearning{n};
    
    % cycle through sessions
    for nn=1:length(dfolders)
        % move to data directory
        disp([num2str(n) '-' num2str(nn)])
        cd(dfolders{nn});
        
        for ii = 1:nSuff
            % move to suffix directory
            p2 = pwd;
            cd(['cueAnalysis' suffName{ii} '\newScoreShuffleTemplate'])
            
            % load cue data
            for jj = 1:nSide
                load([sideName{jj} '\cueCellsUniThresh.mat'])
                if jj==1
                    nCells = max(cueCellsUniThresh.useIdx);
                    cueIdxs = nan(nCells,nSide);
                    curCueData = cueCellsUniThresh;
                else
                    curCueData(jj) = cueCellsUniThresh;
                end
                cueIdxs(cueCellsUniThresh.cueCellRealIdx,jj) = cueCellsUniThresh.cueCellScores;
            end

            % find type with max cue score when cue cells overlap
            [~,mxType] = max(cueIdxs,[],2,'omitnan');
            delIdx = cell(1,nSide);
            for jj = 1:nSide
                delLogic = find(mxType~=jj & ~isnan(cueIdxs(:,jj)));
                
                delIdx{jj} = ismember(curCueData(jj).cueCellRealIdx,delLogic);
            end
            
            % generate combined cue cell structure
            cueCellsAll = struct;
            for jj = 1:nSide
                cueCellsAll.(sideName{jj}).thresh = curCueData(jj).thresh;
                for kk = 1:length(keepFields)
                    cueCellsAll.(sideName{jj}).(keepFields{kk}) = curCueData(jj).(keepFields{kk});
                    cueCellsAll.(sideName{jj}).(keepFields{kk})(delIdx{jj}) = [];
                end
            end
            
            save('cueCellsAll.mat','cueCellsAll')
            cd(p2)
        end
    end
end

cd(p1)

