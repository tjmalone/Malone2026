function licking = predLickActiveAlt(params,envParams,useLaps,varargin)
%% predLickActive
% Calculates predictive licking from active training in virmen using an
% alternate calculation (post-reward exclusion is a fixed length rather
% than entire reward zone.
%
% Inputs:
%       params - virmen log data (output of logParams.m)
%       envParams - environment parameters (reward zone, false alarm zone,
%       and track length)
%       useLaps - range of laps to be analzyed
%       varargin - variable input. Any left empty will be set to default.
%           postRewBuf - post reward buffer (post reward exclusion)
%           minLickDist - minimum lick distance (default of 0 means no
%               exclusion)
%           figOn - whether to plot figures
%
% Outputs:
%       licking - predictive licking data struct
%


%% Process inputs

% read environment parameters
rewStart = envParams.rewStart;
rewStop = envParams.rewStop;
rewWin = rewStop-rewStart;
trackLength = envParams.trackLength;

% read use laps
if isempty(useLaps) || isequal(useLaps,0)
    startLap = 1;
    stopLap = 0;
else
    startLap = useLaps(1);
    stopLap = useLaps(end);
end

% varargin: preRewBuff (cm), postRewBuff(cm),minLickDist(cm),figOn
Defaults = {20,30,0,0};
idx = ~cellfun('isempty',varargin);
Defaults(idx) = varargin(idx);

% set variable values
preRewBuff = Defaults{1};
postRewBuf = Defaults{2};
minLickDist = Defaults{3};
figOn = Defaults{4};



%% Calculate number of runs

y = params.y;

% calculate number of runs
[runIdx,nRuns] = identifyLaps(y,startLap,stopLap);


%% Find reward/lick locations

rewardIdx = params.rewardIdx;
lickIdx = params.lickIdx;

% remove licks after reward delivery
for ii = 1:length(rewardIdx)
    % find the index at end of post reward buffer for each reward
    idx = find(y(rewardIdx(ii):end)>(y(rewardIdx(ii))+postRewBuf),1) + rewardIdx(ii);
    
    % if no points are present in the same lap, find the first index in the
    % next lap (subtract 2 to only select same lap). Only relevant if the
    % end of the reward zone is within buffer distance of the end of the
    % track. Therefore, it is never used in current data
    if isempty(idx)
        idx = find(y(rewardIdx(ii):end)<postRewBuf,1) + rewardIdx(ii)-2;
    end
    
    % if all future points are in the buffer zone, find the last y position
    if isempty(idx)
        idx = length(y);
    end
    
    % find all lick indices between reward and end of post reward buffer
    idxDel = lickIdx>rewardIdx(ii) & lickIdx<idx;
    
    % delete reward related licks
    lickIdx(idxDel) = [];
end


%% Determine licks by lap

% initialize variables
rewIdx = nan(nRuns,1);                  % reward index by lap
licksAll = cell(nRuns,1);               % all lick locations by lap
licksPre = cell(nRuns,1);               % pre reward licks by lap
licksOth = cell(nRuns,1);               % other licks by lap
licksPreLoc = zeros(nRuns,trackLength); % distance limited reward licks by lap
licksOthLoc = zeros(nRuns,trackLength); % distance limited other licks by lap
licksPrePer = zeros(nRuns,1);           % reward licks per distance
licksOthPer = zeros(nRuns,1);           % other licks per distance

for ii = 1:nRuns
    % find reward index in current lap
    curRewIdx = rewardIdx(rewardIdx>=runIdx(ii,1) & rewardIdx<=runIdx(ii,2));
    curRew = y(curRewIdx);
    
    % set predictive lick zone
    if isempty(curRew)
        % if no reward recieved, predictive window is one reward zone
        % length before start of reward zone through end of reward zone
        curRew = rewStop;
        curZone = [rewStop-preRewBuff,rewStop];
    else
        % if reward recieved, predictive window starts at the reward
        % location minus the pre-ward buffer
        curZone = [curRew-preRewBuff,curRew];
        rewIdx(ii) = curRewIdx;
    end
    
    % licks in current lap
    curLick = lickIdx(lickIdx>=runIdx(ii,1) & lickIdx<=runIdx(ii,2));
    licksAll{ii} = y(curLick);
    
    % calculate predictive and other licks
    curLickPre = curLick(y(curLick)>=curZone(1) & y(curLick)<=curZone(2));
    curLickOth = curLick(y(curLick)<curZone(1) | y(curLick)>curZone(2));
    
    % individual licks must be separated by set value (minLickDist)
    if minLickDist~=0
        % predictive and other licks
        for jj = 1:2
            
            % sort y positions
            if jj==1
                [L,I] = sort(y(curLickPre));
            elseif jj==2
                [L,I] = sort(y(curLickOth));
            end
            
            % iteratively remove next lick if is within minLickDist
            kk = 1;
            while kk<length(L)
                if L(kk+1)-L(kk)>minLickDist
                    kk = kk+1;
                else
                    L(kk+1) = [];
                    I(kk+1) = [];
                end
            end
            
            % keep only included licks
            if jj==1
                curLickPre = curLickPre(I);
            elseif jj==2
                curLickOth = curLickOth(I);
            end
        end
    end
    
    % calculate zone lengths
    preLength = diff(curZone);
    othLength = trackLength - diff([curZone(1) curRew+postRewBuf]);
    
    % caluclate licks and y locations
    licksPre{ii} = curLickPre;
    licksOth{ii} = curLickOth;
    licksPreLoc(ii,max(1,round(y(curLickPre)))) = 1;
    licksOthLoc(ii,max(1,round(y(curLickOth)))) = 1;
    
    % calculate licks per distance
    if minLickDist~=0
        licksPrePer(ii) = length(curLickPre)/(preLength/minLickDist)*100;
        licksOthPer(ii) = length(curLickOth)/(othLength/minLickDist)*100;
    else
        licksPrePer(ii) = length(curLickPre)/(preLength)*100;
        licksOthPer(ii) = length(curLickOth)/(othLength)*100;
    end
end

% calculate lick distributions
meanLickPreLoc = mean(licksPreLoc,1);
meanLickOthLoc = mean(licksOthLoc,1);

% calculate predctive licking (percent of laps with more licks per distance
% in predictive zone than in other zone)
licksDiff = licksPrePer-licksOthPer;
predLickByLap = licksDiff>0;
predLickPer = mean(predLickByLap)*100;

% calculate global predicitve licking (number of total predictive licks
% over total number of licks)
sumLickPre = sum(cellfun(@length,licksPre));
sumLickOth = sum(cellfun(@length,licksOth));
perPredInAll= sumLickPre/(sumLickPre+sumLickOth)*100;

% calculate combined lick locations
licksAllComb = cat(1,licksAll{:});

%% Save results

licking.licksPre = licksPre;
licking.licksOth = licksOth;
licking.licksPrePer = licksPrePer;
licking.licksOthPer = licksOthPer;
licking.licksPreLoc = licksPreLoc;
licking.licksOthLoc = licksOthLoc;
licking.meanLickPreLoc = meanLickPreLoc;
licking.meanLickOthLoc = meanLickOthLoc;
licking.licksDiff = licksDiff;
licking.predLickByLap = predLickByLap;
licking.predLickPer = predLickPer;
licking.licksAll = licksAllComb;
licking.perPredInAll = perPredInAll;


%% Plot results

if figOn
    %% Plot lick distribution
    figure; hold on
    
    % plot lick distributions
    plot(smooth(meanLickPreLoc,5),'r-')
    plot(smooth(meanLickOthLoc,5),':b')
    
    
    %% Plot per lap predicitve licking
    
    figure; hold on
    
    % plot y, rewards, and licks
    plot(y,'k');
    plot(rewardIdx,y(rewardIdx),'g.','MarkerSize',10)
    plot(cell2mat(licksPre),y(cell2mat(licksPre)),'r.','MarkerSize',10);
    plot(cell2mat(licksOth),y(cell2mat(licksOth)),'b.','MarkerSize',10);
    
    % mark laps with more predictive licking
    tempRew = rewIdx;
    
    % calculate rewards with predictive licking
    for ii = 1:length(rewIdx)
        if isnan(rewIdx(ii)) && predLickByLap(ii)
            tempRew(ii) = licksPre{ii}(1);
        end
    end
    predRewardIdx = tempRew(predLickByLap);
    
    % plot predictive rewards
    plot(predRewardIdx,y(predRewardIdx),'mo','MarkerSize',10)
    
    % label figure
    xlim([runIdx(1,1) runIdx(end,2)])
    title(['Predictive licking: ' num2str(predLickPer,3),'%']);
    
end

end

