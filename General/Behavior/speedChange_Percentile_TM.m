function [SpeedChangesPercentile] = speedChange_Percentile_TM(logDirectory,d,distAfterReward,numberOfRuns,start_track,stop_track,movemean_value)
%this code is used to analyze the percentile of the speed changes before
%reward delivery throughout the whole track


%input
%(1) logDirectory: the directory of the txt file: for example, this one:
%'C:\Users\guy5\Documents\Lab\systemSetup\analyzingMouseBehaviorOnNewRig\20201117T115335.txt'
%(2) d: the distance before reward delivery for comparison; e.g., 100 for 100 cm
%before reward deliver; then it will compute the mean speed changes within
%every 100 cm and then the percentile of the speed changes of 100 cm before
%reward delivery. for 10 m new and old track, it's better to set d = 150;
%for shorter 4 m track, it's better to set d ~ 50. Don't use d > 200, especially for shorter track. 
%(3) distAfterReward:the distance after reward delivery which will be
%excluded from speed calculation. e.g., 30 for 30 cm. So 30 cm after reward
%delivery will not be used to calcaulte speed changes
% Notice: if the distance after reward delivery + distAfterReward is
% shorter than d, then the velocity after reward delivery will not be
% included to compute. For example, in the new 10 m track, the reward
% delivery is ~890, if distAfterReward = 30, 890+30 = 920, so the rest of
% the track is 1000-920 = 80. If d is >=80, then the script will ignore
% 920-1000 for calculation.
%(4) numberOfRuns: this is to determine how many runs (counted from the
%beginning) are used for the analysis. The reason for doing this is that
%sometimes the animal's behavior does not look good after a while so only
%use the first couple of runs to evaluate the behavior. If numberOfRuns IS
%0, all the runs will be used. 0.5 means use half. Otherwise, we will just use the indicated
%number runs.
%(5) start_track: the start position of the track you want to include to
%calculate the speed percentile. See (6) below for example
%(6) stop_track: the end position of the track for calculating the speed
%percentile. For example,for the old 10 m environment, the script will only
%take the second reward delivery to calculate. So you may want to exclude
%the first 50 cm of the track for the calculation as the first reward
%delivery is 35 cm. For this case, the start = 50, the stop = 1000; For the
%new 10 m environment, start = 0, stop = 1000.
%(7) movemean_value: the average moving window for smoothen the velocity; 0
%means no average is applied
% Example: 
% new 10 m track: speedChange_Percentile('20201117T115335.txt',150,30,0,0,1000,0) 
% old 10 m track: speedChange_Percentile('20201117T115335.txt',150,30,0,50,700,0)
% 4 m track: speedChange_Percentile('20201117T115335.txt',70,30,0,0,400,50)
% Notes: for every new environment, the parameteres may need to modify to
% give reasonable results. One way to fine tune the parameteres is to run
% multiple different combinations of bin size, start and stop locations
% with multiple mice. The idea output should give a wide range between ~40
% and ~80. If most of the mice have extremely low values (e.g., ~30) or
% high values, it may indicate that the parameteres need to be changed.
% To determine the start and stop locations, it's better to first plot the
% velocity changes along the track from many mice to get an idea where
% should not be included for calculation. For tracks with more than one
% reward location (e.g., old 10 m track), extra cares are required for the
% start and stop locations as this script only considers the last reward location.
% For example, several mice started slowing down at 700 cm in the old 10 m
% track for the reward at 35 cm. In this case, it's better to have
% stop_track = 700.
%
% 
%output
% SpeedChangesPercentile.percentile: the percentile of the speed changes
% before reward delivery for each run
% SpeedChangesPercentile.meanCentile: the average percentile of the speed
% changes before reward delivery
% SpeedChangesPercentile.speedChanges: speed changes throughout the whole
% track for each run
% SpeedChangesPercentile.deceleration: the speed changes before reward
% location

%% load the log file
data = readLog(logDirectory);
position = data(:,3);
reward = data(:,8);
time = (data(:,1)-data(1,1))*24*60;

% % truncate each run
% C=diff(position);
% C=[C;0];
% A=C<=-max(position)*0.2;%0 and 1. Indices of the teleportation point will be 1; in the analysis,
% 
% start_contiguous = contiguous(A,[1]);
% 
% begin = [1; start_contiguous{2}(:,1)];% begin is the start for each run, right before the mouse starts move

% truncate each run
start = find(position==0);
start = [start; start(end)+1000];
start_pause = diff(start); % when the mouse pauses at 0 position
real_start = find(start_pause>500);

begin = start(real_start);% begin is the start for each run, right before the mouse starts move


%%

if numberOfRuns == 0
    if length(begin)>1
        runs = length(begin)-1;
    else
        runs = length(begin);
    end
elseif numberOfRuns == 0.5
    runs = round(length(begin)/2);
else
    runs = numberOfRuns;
end

for j = 1:runs % did not use the last run in case the mouse did not complete
    if runs >1
        run_position = position(begin(j):begin(j+1)-1);
        time_segment = time(begin(j):begin(j+1)-1);
    else
        index = find(position==max(position));
        run_position = position(begin(j):index);
        time_segment = time(begin(j):index);
    end
%     velocity = diff(run_position');
%     velocity(velocity<-60) = 0;
%     velocity = velocity*60;
%     run_position = movmean(run_position,10);
    displacement = diff(run_position');
    displacement(displacement<-60) = 0;
    time_interval = diff(time_segment);
    velocity = (displacement./time_interval')/60;
    
    %****************apply moving averaging****************
    if movemean_value~=0
        velocity = movmean(velocity,movemean_value);        
    end
    %*******************************************************
    
    % find the indices before actual reward location
    if runs>1
        reward_run = reward(begin(j):begin(j+1)-1);
    else
        reward_run = reward(begin(j):index);
    end
    reward_run_index = find(reward_run==1);
    
    if reward_run_index    
        reward_run_index = reward_run_index(end);
        run_position_beforeReward = run_position(1:reward_run_index);
        velocity_beforeReward = velocity(1:reward_run_index);
        
        % find the speed after actual reward location, excluding moving
        % backwards
        reward_zone = run_position(reward_run_index) + distAfterReward;
        reward_zone = reward_zone(end);
        actual_reward_zone_idx = find(run_position - reward_zone>0,1,'first'); % find the actual reward zone
        
        run_position_afterReward = run_position(actual_reward_zone_idx:end);
        velocity_afterReward = velocity(actual_reward_zone_idx:end);
        
        before = start_track:1:run_position(reward_run_index)-90; % 90 cm before 
        % reward delivery will not used to compare to the speed changes before reward delivery
        after = reward_zone:1:stop_track;
        
        % calculate how many samples virmen collected within the specified
        % distance
        distance = d;
        samples(j) = sum(run_position_beforeReward >= run_position_beforeReward(end)-distance);
        
        % calculate the velocity changes within the specified 
        y = velocity_beforeReward(end-samples(j)+1:end);
        SpeedChangesPercentile.deceleration(j) = nanmean(diff(y));
        
        index = 1;
        for i = start_track:before(end)-distance+1
            start_position = i;
            stop_position = i+distance-1;
            start_idx = find(run_position_beforeReward>=start_position,1);
            stop_idx = find(run_position_beforeReward>=stop_position,1);
            if stop_idx-start_idx+1 >=1 % check if there are at least two points 
                y = velocity_beforeReward(start_idx:stop_idx);
                BeforeSpeedChanges(j,index) = nanmean(diff(y));
                index = index+1;
            else
                BeforeSpeedChanges(j,index) = nan;
                index = index+1;
            end
        end
        
        index = 1;
        if stop_track-reward_zone > d
            for i = reward_zone:after(end)-distance+1
                start_position = i;
                stop_position = i+distance-1;
                start_idx = find(run_position_afterReward>=start_position,1);
                stop_idx = find(run_position_afterReward>=stop_position,1);
                if stop_idx-start_idx+1 >=1 % check if there are at least two points 
                    y = velocity_afterReward(start_idx:stop_idx);
                    AfterSpeedChanges(j,index) = nanmean(diff(y));
                    index = index+1;
                else
                    AfterSpeedChanges(j,index) = nan;
                    index = index+1;
                end
            end
        else
            AfterSpeedChanges(j,index) = nan;
            index = index+1;            
        end   
        
    else
        BeforeSpeedChanges(j,:) = nan;
        AfterSpeedChanges(j,:) = nan;
        SpeedChangesPercentile.deceleration(j) = nan;
    end
end
for j = 1:runs
    % combine before and after deceleration in time, then calculate
    % the percentile    
    if isnan(BeforeSpeedChanges(j,1))
        BeforeSpeedChanges(j,:) = nan;
    end
    try
        if stop_track-reward_zone > d
            SpeedChangesPercentile.speedChanges(j,:) = [BeforeSpeedChanges(j,:) AfterSpeedChanges(j,:)]*(-1);
        else
            SpeedChangesPercentile.speedChanges(j,:) = BeforeSpeedChanges(j,:)*(-1);
        end
        nless = sum(SpeedChangesPercentile.speedChanges(j,:) < SpeedChangesPercentile.deceleration(j)*(-1));
        nequal = sum(SpeedChangesPercentile.speedChanges(j,:) == SpeedChangesPercentile.deceleration(j)*(-1));
        SpeedChangesPercentile.percentile(j) = 100 * (nless + 0.5*nequal) / length(SpeedChangesPercentile.speedChanges(j,:));
    catch
        SpeedChangesPercentile.percentile(j) = nan;
        continue
    end
    
end
try
    SpeedChangesPercentile.meanCentile = nanmean(SpeedChangesPercentile.percentile);
    SpeedChangesPercentile.runs = runs;
%     disp(['percentile of speed changes before reward delivery = ',num2str(SpeedChangesPercentile.meanCentile)])
catch
%     disp('less than 2 runs')
    SpeedChangesPercentile.meanCentile = nan;
    SpeedChangesPercentile.runs = nan;
end


end