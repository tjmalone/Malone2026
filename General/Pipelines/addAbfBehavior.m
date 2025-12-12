%% addAbfBehavior
% adds behavior information to abf files


%% Define paramters

clear; close all; clc

% define base folder
p1 = 'D:\AD_Project\imagingData\data';
cd(p1)

% identify analysis folders
folds = findSubF('suite2p',4,[],0);


%% Generate abf

for ff = 1:length(folds)
    %% Plot behavior
    %5th column is lick
    %6th column is reward

    cd(folds{ff})
    disp(ff)

    % load abf and behavior file
    load('abf.mat','abf');
    try
        abf.rewardFrame;
        continue
    catch
    end

    load('..\..\m.mat')

    % identify synchronized points
    thresh=4.8;
    a=m(:,2)>thresh;
    b=contiguous(a,1);
    c=b{1,2};
    middleFrameIdx=ceil(mean(c,2));

    % calculate all imaging indices
    cImaging=b{1,2};

    % calculate all lick indices
    lickThresh=0.8;
    lick=m(:,5)>lickThresh;
    lickAllIdx=find(lick);

    % calculate all reward indices
    rewardThresh=0.8;
    reward=m(:,6)>rewardThresh;
    b=contiguous(reward,1);
    c=b{1,2};
    if ~isempty(c)
        rewardAllIdx = c(:,1);
    else
        rewardAllIdx = c;
    end

    nOriginalFrames = length(abf.imageIndex);

    ii = str2double(folds{ff}(end-9));

    % calculate current imaging indices
    curIdx = middleFrameIdx(ii:2:end);

    cImagingFrame = [];
    cImagingIdx = [];
    for n = ii:2:nOriginalFrames*2
        cImagingFrame = [cImagingFrame; (cImaging(n,1):cImaging(n,2))'];
        cImagingIdx = [cImagingIdx; n*ones(diff(cImaging(n,:))+1,1)];
    end


    %% Calculate licks
    % Calculate the overlap betweeen the lick trace and imaging trace. Uses
    % all lick indices (not just those that coincide with peak lick signal).
    % lickIdx - all lick indices matched to the nearest imaging frame.
    % lickFrame - lick indices that overlap an imaging frame

    % calculate nearest image frame to licks
    closestIndex = zeros(length(lickAllIdx),1);
    for n = 1:length(lickAllIdx)
        [~,closestIndex(n)] = min(abs(lickAllIdx(n)-curIdx));
    end

    % calculate unique lick indices up to max frame number
    closestIndexUnq = unique(closestIndex);
    closestIndexUnq(closestIndexUnq>=nOriginalFrames) = [];
    closestIndexUnq(closestIndexUnq<1) = [];

    % save lick indices (all lick indices matched to the nearest frame)
    lickIdx = false(nOriginalFrames,1);
    lickIdx(closestIndexUnq) = 1;
    abf.lickIdx = lickIdx;

    % calculate lick overlaps
    [~,~,lickInt] = intersect(lickAllIdx,cImagingFrame);
    lickFrameIdx = unique(cImagingIdx(lickInt));

    % save lick frames (lick indices that overlap an imaging frame)
    lickFrame = zeros(nOriginalFrames,1);
    lickFrame(lickFrameIdx) = 1;
    abf.lickFrame = lickFrame;


    %% Calculate rewards
    % Calculate the overlap betweeen the reward trace and imaging trace.
    % Uses all reward indices (not just those that coincide with start of
    % reward signal).
    % rewardIdx - all reward indices matched to the nearest imaging frame.
    % rewardFrame - reward indices that overlap an imaging frame

    % calculate nearest image frame to rewards
    closestIndex = zeros(length(rewardAllIdx),1);
    for n = 1:length(rewardAllIdx)
        [~,closestIndex(n)] = min(abs(rewardAllIdx(n)-curIdx));
    end

    % calculate unique reward indices up to max frame number
    closestIndexUnq = unique(closestIndex);
    closestIndexUnq(closestIndexUnq>=nOriginalFrames) = [];
    closestIndexUnq(closestIndexUnq<1) = [];

    % save reward indices (all reward indices matched to the nearest frame)
    rewardIdx = false(nOriginalFrames,1);
    rewardIdx(closestIndexUnq) = 1;
    abf.rewardIdx = rewardIdx;

    % calculate reward overlaps
    [~,~,rewardInt] = intersect(rewardAllIdx,cImagingFrame);
    rewardFrameIdx = unique(cImagingIdx(rewardInt));

    % save reward frames (reward indices that overlap an imaging frame)
    rewardFrame = zeros(nOriginalFrames,1);
    rewardFrame(rewardFrameIdx) = 1;
    abf.rewardFrame = rewardFrame;


    %% Plot behavior

    % calculate y positions
    ylick = abf.y(lickIdx);
    yreward = abf.y(rewardIdx);

    % calculate times
    t = abf.t;
    tlick = abf.t(lickIdx);
    treward = abf.t(rewardIdx);

    scaleX = 1;
    figure; hold on
    plot(t(1:scaleX:end,1),abf.y(1:scaleX:end),'k');
    plot(tlick,ylick,'r.','MarkerSize',15);
    plot(treward,yreward,'g.','MarkerSize',15);

    saveas(gcf,'behavior.fig');
    close

    %% Save

    save('abf.mat','abf')

    abfFake = abf;
    save('abfFake.mat','abfFake')
end

cd(p1)
