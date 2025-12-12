%%  plotBehaviorForFigs
% Based on behaviorAD_all, which should be run first to generate the proper
% precursor data. Plots the graphs to be used for the AD paper formatted
% specifically for the paper. Statistics and legends should be added after
% the fact based on the original plots from behaviorAD_all (to avoid issues
% with combining by sex). Calculates and plots speed match plots
%


%% Load behavior data
% To skip recalculation, start running from here. All previously calculated
% data is loaded in allMiceData*.mat. This file name must match the current
% save output above.
%

% clear all variables
clear; close all; clc

% reset directory
cd('D:\AD_Project\Behavior')
p = pwd;

% load processed behavioral data
load('data\allMiceDataB5.mat','B')

% load genotypes and sexes
load('groupIDs.mat')

% define legend
legs = groupIDs;
nGroups = length(groups);
nSexes = length(sexes);

% define sex groups
sGroups = {1,[2 3]};
sSvLabels = {'all','bySex'};
nSGroups = length(sGroups);

% define plot colors
colors = {[0 0 1],[1 0 0]};
colors4 = {[102 105 255]/255,[255 102 102]/255;[0 0 175]/255,[175 0 0]/255};

% sefine color shift by sex
colShift = [1 0.5];

% set day categories (FE,NE,FER,NER)
dayCats = {1:4,5:14,15,16};
nCats = length(dayCats);

% set save folder name
svFile = [p '\Figures\manuscriptFigures\speedMatch'];

if ~isfolder(svFile)
    mkdir(svFile);
    for ss = 1:nSGroups
        mkdir(svFile)
        mkdir([svFile '\data'])
    end
end

% close all figures
clAll = 0;


%% Define valid speed indices by lap

% define types
useTypes = [1 2];
nTypes = length(useTypes);

% define mice
nMice = length(B);
envTypes = [1 2*ones(1,10)];
nDays = 11;
dayNames = ['FE' string(1:10)];
dayCats = {1,2:11};

% define matching threshold
maxDist = 100;

% initialize constrain function
constrainRange2D = @(x) max(min(mean(x,1),1-1/(size(x,1)*2)),1/(size(x,1)*2));

% define d' indices
idxDGoal = cell(1,2);
idxDInclude = cell(1,2);
for ii = 1:2
    curIdx = B(1).type{useTypes(ii)}(1);
    idxDGoal{ii} = B(1).idxGoal(curIdx);
    idxDInclude{ii} = B(1).idxInclude{curIdx};
end

% define variables and labels
fldNames = {'vel','success','globalDP'};
nFields = length(fldNames);
yLabs = {'Mean velocity (cm/s)','Success rate (% of runs)','Global d'''};
ttls = {'Mean Velocity','Success Rate','Global d'''};

% initialize pre data struct
preData = struct();
for ff = 1:nFields
    preData.(fldNames{ff}) = cell(nMice,nDays);
end

for mm = 1:nMice
    % define enviornment indices
    idx = cat(1,B(mm).type{useTypes});

    % get velocity info
    curVel = B(mm).avgVelMovingRawPerLap(idx);

    % get success info
    curSuccess = B(mm).successPerLap(idx);
    curSuccessRoll = B(mm).successRoll(idx);

    % store current data
    preData.vel(mm,:) = curVel;
    preData.success(mm,:) = curSuccess;
    preData.successRoll(mm,:) = curSuccessRoll;
end

% concatenate all lap velocities
velCat = cat(1,preData.vel{:});

% initialize figure
figure;
tiledlayout(nSexes,nFields)
sgtitle('Velocity matched behavior (by lap)')

% calculate matching for each sex
for ss = 1:nSexes
    % initialize behavior struct
    matchData = struct();
    matchData.vel = cell(nDays,2);
    matchData.success = cell(nDays,2);
    matchData.globalDP = cell(nDays,2);

    % get group indices
    gpIdx = cell(1,2);
    gpMouseMap = cell(nDays,2);
    gpRunMap = cell(nDays,2);
    gpVelCat = cell(nDays,2);
    matchIdx = cell(nDays,2);
    for gg = 1:2
        gpIdx{gg} = intersect(groups{gg},sexes{ss});

        for dd = 1:nDays
            % create mouse id map
            gpMouseID = arrayfun(@(ii) repmat(ii, length(preData.vel{gpIdx{gg}(ii),dd}),1), ...
                (1:numel(gpIdx{gg}))','UniformOutput', false);
            gpMouseMap{dd,gg} = cat(1, gpMouseID{:});

            % create mouse run map
            gpRunID = arrayfun(@(ii) 1:length(preData.vel{gpIdx{gg}(ii),dd}), ...
                (1:numel(gpIdx{gg}))','UniformOutput', false);
            gpRunMap{dd,gg} = cat(2,gpRunID{:})';

            % concatenate group velocity
            gpVelCat{dd,gg} = cat(1,preData.vel{gpIdx{gg},dd});
        end
    end

    % perform daily velocity matching
    for dd = 1:nDays
        % define d' indices
        curDGoal = idxDGoal{envTypes(dd)};
        curDInclude = idxDInclude{envTypes(dd)};

        % calculate matching pairs
        costMatrix = pdist2(gpVelCat{dd,1},gpVelCat{dd,2},'cityblock');

        % apply distance threshold to costs beyond threshold to Inf
        costMatrix(costMatrix > maxDist) = Inf;

        % perform matching under constraint
        pairs = matchpairs(costMatrix,1e6);

        % apply distance threshold
        pairDists = abs(gpVelCat{dd,1}(pairs(:,1))-gpVelCat{dd,2}(pairs(:,2)));
        validPairs = pairs(pairDists<=maxDist,:);

        % recreate cell array
        for gg = 1:2
            % define remap mice
            validIdx = validPairs(:,gg);
            validRuns = gpRunMap{dd,gg}(validIdx);
            validMap = gpMouseMap{dd,gg}(validIdx);

            % match indices to mouse
            useIdx = cell(length(gpIdx{gg}),1);
            for m = unique(validMap)'
                useIdx{m} = validRuns(validMap==m);
            end
            matchIdx{dd,gg} = useIdx;

            % calculate new success rate by mouse
            curSuccess = cellfun(@(x,y) x(y), preData.success(gpIdx{gg},dd),useIdx,...
                'UniformOutput',false);
            matchData.success{dd,gg} = cellfun(@(x) mean(x,'omitnan'),curSuccess);

            % calculate new velocity by mouse
            curVel = cellfun(@(x,y) x(y), preData.vel(gpIdx{gg},dd),useIdx,...
                'UniformOutput',false);
            matchData.vel{dd,gg} = cellfun(@(x) mean(x,'omitnan'),curVel);

            % calculate new global d'
            curSuccessRoll = cellfun(@(x,y) x(y,:), preData.successRoll(gpIdx{gg},dd),useIdx,...
                'UniformOutput',false);
            conSuccessRoll = cellfun(@(x) constrainRange2D(x),curSuccessRoll,'UniformOutput',false);

            % calculate d-prime
            dprimeRoll = cellfun(@(x) norminv(x(curDGoal))-norminv(x(curDInclude)),...
                conSuccessRoll,'UniformOutput',false);
            matchData.globalDP{dd,gg} = cellfun(@(x) mean(x,'omitnan'), dprimeRoll);
        end
    end

    %% Plot updated behavior

    for ff = 1:nFields
        % get current data
        plotData = matchData.(fldNames{ff})';

        % plot figure
        nexttile; hold on
        [pAnova,pMC] = errorSig(dayNames,plotData,colors,groupIDs,dayCats,0);

        % set labels
        title([ttls{ff} ': ' sexIDs{ss} ' mice'])
        xlim([0.5 11.5])
        xlabel('Session')
        ylabel(yLabs{ff})
    end

end

savefig([svFile '\velocityMatchbyLap.fig'])


%% Define valid speed indices by session

% define types
useTypes = [1 2];

% define mice
nMice = length(B);
envTypes = [1 2*ones(1,10)];
nDays = 11;
dayNames = ['FE' string(1:10)];
dayCats = {1,2:11};

% define matching threshold
maxDist = 5;

% define variables and labels
fldNames = {'vel','success','globalDP'};
nFields = length(fldNames);
yLabs = {'Mean velocity (cm/s)','Success rate (% of runs)','Global d'''};
yLims = [10 40;0 110;-1 4.5];
ttls = {'Mean Velocity','Success Rate','Global d'''};

% initialize pre data struct
preData = struct();
for ff = 1:nFields
    preData.(fldNames{ff}) = zeros(nMice,nDays);
end

for mm = 1:nMice
    % define enviornment indices
    idx = cat(1,B(mm).type{useTypes});

    % get velocity info
    curVel = B(mm).avgVelMovingRaw(idx);

    % get success info
    curSuccess = B(mm).perSuccess(idx);
    curDPrime = B(mm).dprimeGlobal(idx);

    % store current data
    preData.vel(mm,:) = curVel;
    preData.success(mm,:) = curSuccess;
    preData.globalDP(mm,:) = curDPrime;
end

% concatenate all lap velocities
velCat = preData.vel(:);

% initialize figure
figure;
tiledlayout(nSexes,nFields)
sgtitle('Velocity matched behavior (by mouse)')

% initialize matched mouse store
validPairsAll = cell(nSexes,nDays);
storeData = struct();

% calculate matching for each sex
for ss = 1:nSexes
    % initialize behavior struct
    matchData = struct();
    for ff = 1:nFields
        matchData.(fldNames{ff}) = cell(nDays,2);
    end

    % get group indices
    gpIdx = cell(1,2);
    gpVelCat = cell(nDays,2);
    matchIdx = cell(nDays,2);
    for gg = 1:2
        gpIdx{gg} = intersect(groups{gg},sexes{ss});

        for dd = 1:nDays
            % concatenate group velocity
            gpVelCat{dd,gg} = preData.vel(gpIdx{gg},dd);
        end
    end

    % perform daily velocity matching
    for dd = 1:nDays
        % calculate matching pairs
        costMatrix = pdist2(gpVelCat{dd,1}, gpVelCat{dd,2}, 'cityblock');

        % apply distance threshold to costs beyond threshold to Inf
        costMatrix(costMatrix > maxDist) = Inf;

        % perform matching under constraint
        pairs = matchpairs(costMatrix,1e6);

        % apply distance threshold
        pairDists = abs(gpVelCat{dd,1}(pairs(:,1))-gpVelCat{dd,2}(pairs(:,2)));
        validPairs = pairs(pairDists<=maxDist,:);
        validPairsAll{ss,dd} = validPairs;

        % recreate cell array
        for gg = 1:2
            % calculate new velocity by mouse
            matchData.vel{dd,gg} = preData.vel(gpIdx{gg}(validPairs(:,gg)),dd);

            % calculate new success rate by mouse
            matchData.success{dd,gg} = preData.success(gpIdx{gg}(validPairs(:,gg)),dd);

            % calculate new global d'
            matchData.globalDP{dd,gg} = preData.globalDP(gpIdx{gg}(validPairs(:,gg)),dd);
        end
    end


    %% Plot updated behavior
    
    % store final matched data
    storeData.(sexIDs{ss}) = matchData;

    for ff = 1:nFields
        % get current data
        plotData = matchData.(fldNames{ff})';

        % plot figure
        nexttile; hold on
        [pAnova,pMC] = errorSig(dayNames,plotData,colors,groupIDs,dayCats,0);

        % set labels
        title([ttls{ff} ': ' sexIDs{ss} ' mice'])
        xlim([0.5 11.5])
        ylim(yLims(ff,:))
        xlabel('Session')
        ylabel(yLabs{ff})
        legend('off')
    end


end

save([svFile '\data\velocityMatchbySession.mat'],'validPairsAll','storeData')
savefig([svFile '\velocityMatchbySession.fig'])

