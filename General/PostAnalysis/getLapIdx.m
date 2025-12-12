function [runIdx,lapIdx,rewardIdx] = getLapIdx(m,rewardThresh)
%% getLapIdx
% use a .m or .md behavior sync file to identify lap indices. Outputs
% include the lap numbers (lapIdx), the indices of each lap (runIdx), and
% the laps with rewards (rewardIdx).

if nargin<2 || isempty(rewardThresh)
    rewardThresh = 0.8;
end

% set behavior variables
laps = m(:,4)*100;
y = m(:,3)*100;
rewards = m(:,6)>rewardThresh;

% find lap ends
N = contiguous(round(laps));

% correct and concatenate lap indices
nRuns = size(N,1);
NN = zeros(nRuns,2);
for jj = 1:size(N,1)
    curN = N{jj,2};

    if size(curN,1)==1
        curCorrect = curN;
    else
        curCorrect = [];
        for kk = 1:size(curN,1)
            if diff(curN(kk,:))>100
                if isempty(curCorrect)
                    curCorrect = curN(kk,:);
                else
                    curCorrect(2) = curN(kk,2);
                end
            end
        end
    end

    NN(jj,:) = curCorrect;
end

% correct for teleportation
for jj = 1:nRuns
    C = diff(y(NN(jj,1):NN(jj,2)));

    % remove teleportation at start of lap
    if C(1)<=-max(y)*0.2
        idxFirst = find(C>-max(y)*0.2,1,'first');
        NN(jj,1) = NN(jj,1) + idxFirst - 1;
    end

    % remove teleportation at end of lap
    if C(end)<=-max(y)*0.2
        idxLast = find(C>-max(y)*0.2,1,'last');
        NN(jj,2) = NN(jj,2)+idxLast-diff(NN(jj,:));
    end
end

% save output
lapIdx = unique(round(laps))+1;
runIdx = NN;
if length(lapIdx)~=size(runIdx,1)
    error('Size mismatch')
end

% get reward laps
rewardIdx = zeros(nRuns,1);
for jj = 1:nRuns
    if any(rewards(runIdx(jj,1):runIdx(jj,2)))
        rewardIdx(jj) = 1;
    end
end

end

