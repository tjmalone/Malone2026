function [params] = logParams(logData, varargin)
%this code extract parameters from virmen logs on linear track: on my new
%rig

%input:
%logDirectory: the directory of the txt file: for example, this one:
%'C:\Users\guy5\Documents\Lab\\20201117T115335.txt'
%
%output:
%params.y:y position
%params.lickIdx=licking indices
%params.rewardIdx=the indices when the reward is delivered



%%%%%%%
% 1/4 - change parameter -> logDirectory = logData to be able to remove
% first minutes of training

% a=readLog(logDirectory);
a = logData;
t = (logData(:,1)-logData(1,1))*24*60*60; % time in seconds
y=a(:,3);%%%%%%%%%%%%%%%%%%%y position
b=a(:,9);%licking
r=a(:,8);%reward

b(b<0.7)=0;
b(b>=0.7)=1;%turning licking signal to 1 and 0

if any(b>0.7)
    n=contiguous(b,1); %find contiguous 1
    nn=n{1,2};
    lickIdx=nn(:,1);%%%%%%%%%%%%%%%%%%%licking indices
    yLick=y(lickIdx);%y position of licking
else
    lickIdx=nan;
    yLick=nan;
end

rewardIdx=find(r);%%%%%%%%%%%%%%%%%reward indices
yReward=y(rewardIdx);%y position of reward

params.t=t;
params.y=y;
params.lickIdx=lickIdx;
params.yLick=yLick;
params.rewardIdx=rewardIdx;
params.yReward=yReward;

if ~isempty(varargin)
    if varargin{1}
        figure,plot(y,'k');
        hold on
        plot(rewardIdx,yReward,'g.','MarkerSize',10)
        hold on
        plot(lickIdx,yLick,'r.')
    end
end


end

%%
function runs = contiguous(A,varargin)
%   RUNS = CONTIGUOUS(A,NUM) returns the start and stop indices for contiguous
%   runs of the elements NUM within vector A.  A and NUM can be vectors of
%   integers or characters.  Output RUNS is a 2-column cell array where the ith
%   row of the first column contains the ith value from vector NUM and the ith
%   row of the second column contains a matrix of start and stop indices for runs
%   of the ith value from vector NUM.    These matrices have the following form:
%
%   [startRun1  stopRun1]
%   [startRun2  stopRun2]
%   [   ...        ...  ]
%   [startRunN  stopRunN]
%
%   Example:  Find the runs of '0' and '2' in vector A, where
%             A = [0 0 0 1 1 2 2 2 0 2 2 1 0 0];
%
%   runs = contiguous(A,[0 2])
%   runs =
%           [0]    [3x2 double]
%           [2]    [2x2 double]
%
%   The start/stop indices for the runs of '0' are given by runs{1,2}:
%
%           1     3
%           9     9
%          13    14
%
%   RUNS = CONTIGUOUS(A) with only one input returns the start and stop
%   indices for runs of all unique elements contained in A.
%
%   CONTIGUOUS is intended for use with vectors of integers or characters, and
%   is probably not appropriate for floating point values.  You decide.
if prod(size(A)) ~= length(A)   %#okay
    error('A must be a vector.')
end

if isempty(varargin)
    num = unique(A);
else
    num = varargin{1};
    if prod(size(num)) ~= length(num)   %#okay
        error('NUM must be a scalar or vector.')
    end
end

for numCount = 1:length(num)
    
    indexVect = find(A(:) == num(numCount));
    shiftVect = [indexVect(2:end);indexVect(end)];
    diffVect = shiftVect - indexVect;
    
    % The location of a non-one is the last element of the run:
    transitions = (find(diffVect ~= 1));
    
    runEnd = indexVect(transitions);
    runStart = [indexVect(1);indexVect(transitions(1:end-1)+1)];
    
    runs{numCount,1} = num(numCount);   %#okay
    runs{numCount,2} = [runStart runEnd];   %#okay
end
end

