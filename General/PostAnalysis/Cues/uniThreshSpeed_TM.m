% find and apply speed threshold

clear; clc; close all;

p1 = 'D:\AD_Project\imagingData';
cd(p1)

svFolder = 'data_uniThresh\speedLearning';

% load folders to analyze
load('foldersLearning.mat')

% set environment types
envGroups = {1:11};


%% Collect shuffles

suffName = {'dfof','sig'};
suffName2 = {'dfof','dfof_sig'};
nSuff = length(suffName);

% shuffle data
allShuffles = struct();
for gg = 1:length(envGroups)
    for jj = 1:nSuff
        allShuffles(gg).(suffName{jj}) = [];
    end

end

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

            %cycle suffixes
            for jj = 1:nSuff
                % load speed data
                load(['speed_' suffName2{jj} '\speed.mat']);

                % add speed data to full data
                A=speed.shuffleScores;
                A=reshape(A,[],1);
                allShuffles(gg).(suffName{jj})(end+1:end+length(A)) = A;
            end
        end
    end
end

cd(p1);

save([svFolder '\allShuffles.mat'],'allShuffles')



%% Process shuffles

% threshold values
thresholdVal = struct();

if ~isvarname('allShuffles')
    load([svFolder '\allShuffles.mat'])
end
prctileCutoff = [1 99];

% cycle enviornment groups
for gg = 1:length(envGroups)
    %cycle suffixes
    for jj = 1:nSuff
        % identify threshold
        thresholdVal(gg).(suffName{jj}) =...
            prctile(allShuffles(gg).(suffName{jj}),prctileCutoff);
    end
end

% clear allShuffles


%% Apply thresholds

% speed cell percentage
perSpeed = struct();
for gg = 1:length(envGroups)
    for jj = 1:nSuff
        perSpeed(gg).(suffName{jj}) = [];
    end
end

% speed cell percentage
scoreSpeed = struct();
for gg = 1:length(envGroups)
    for jj = 1:nSuff
        scoreSpeed(gg).(suffName{jj}) = [];
    end
end

% keep fields
keepFields = 1;

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

            % cycle suffixes
            for jj = 1:nSuff
                % change directory
                cd(['speed_' suffName2{jj}]);

                % current threshold
                curThresh = thresholdVal(gg).(suffName{jj});

                % load speed cells
                load('speed.mat');
                fNames = fieldnames(speed);

                % speed cell information structure
                speedCellsUniThresh = struct;

                % copy unthresholded info
                for kk = 1:length(keepFields)
                    speedCellsUniThresh.(fNames{keepFields(kk)}) =...
                        speed.(fNames{keepFields(kk)});
                end

                % add thresholded info
                speedCellsUniThresh.percentile = prctileCutoff;
                speedCellsUniThresh.thresh = curThresh;
                speedCellsUniThresh.speedCellPost = find(speed.scores>curThresh(2));
                speedCellsUniThresh.speedCellNegt = find(speed.scores<curThresh(1));
                

                % save speed cell structure
                save('speedCellsUniThresh.mat','speedCellsUniThresh');

                % update speed cell percent structure
                perSpeed(gg).(suffName{jj})(end+1,1) = (length(speedCellsUniThresh.speedCellPost)+...
                    length(speedCellsUniThresh.speedCellNegt))/N;

                % update speed cell percent structure
                scoreSpeed(gg).(suffName{jj})(end+1:end+length(speed.scores),1) =...
                    speed.scores;

                cd(p2)
            end
        end
    end
end

cd(p1)

% calculate mean and SEM for speed percentages
perSpeedMean = struct();
perSpeedSEM = struct();
for gg = 1:length(envGroups)
    for jj = 1:nSuff
        perSpeedMean(gg).(suffName{jj}) = mean(perSpeed(gg).(suffName{jj}),'omitnan');
        perSpeedSEM(gg).(suffName{jj}) = nansem(perSpeed(gg).(suffName{jj}),1);
    end
end

% save threshold information
save([svFolder '\speedInfo.mat'],'thresholdVal','perSpeed','perSpeedMean','perSpeedSEM','scoreSpeed')

