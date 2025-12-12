function quantData = distrQuant2(dist1,dist2,rawZones,absZones,combZones,selfZones)
%% distrQuant
% Performs quantification of two arbitrary spatial distributions.
% Specifically, calculates various comparison statistics between the the
% two distributions and some individual distribution statistics: raw
% differences, absolute differences, comparison correlation, and self
% correlation. Outputs data as a single struct.
%
% Note: Primary testing was only performed using numerical array inputs, so
% additional testing is necessary to confirm functionality with
% cells/structs.
%
% Inputs:
%   dist1/dist2 - input distribution data. Can be a numerical or cell
%       arrays of numerical arrays or structS of cell or numerical arrays.
%       At lowest level, numerical arrays must be identical in size along
%       second dimension (spatial dimension). The two distributions must
%       have full identical sizes
%   rawZones - indices at which to compare the raw difference between the
%       distributions. Can be numerical or cell array or numerical arrays.
%       Default is to not analyze.
%   absZones - indices at which to compare the absolute difference between
%       the distributions. Can be numerical or cell array or numerical
%       arrays. Default is to not analyze.
%   combZones - indices at which to calculate the correlation between the
%       distributions. Can be numerical or cell array or numerical arrays.
%       Default is to not analyze.
%   selfZones - indices at which to calculate the correlation between two
%       segments of each distribution. Must be an n-by-2 cell array with
%       the two zones in the second dimension. The number of elements for
%       each comparison must be identical. Default is to not analyze.
%
% Outputs:
%   quantData - output structure. Each field represents one calculation type
%       and is the same data type/size as the input distribution data.
%


%% Process inputs

% process rawZones
if nargin<3 || isempty(rawZones)
    rawCalc = 0;
    rawN = 0;
else
    rawCalc = 1;

    % check rawZones type
    if isnumeric(rawZones)
        % check rawZones size
        szCur = size(rawZones);
        if length(szCur)~=2 || sum(szCur>1)>1
            error('Invalid sizes')
        end
        rawZones = {rawZones};

    elseif iscell(rawZones)
        % check rawZones size
        szCur = size(rawZones);
        if length(szCur)~=2 || sum(szCur>1)>1
            error('Invalid sizes')
        elseif szCur(2)>=1
            rawZones = rawZones';
        end
    else
        error('Invalid type')
    end
    rawN = length(rawZones);
end

% process absZones
if nargin<4 || isempty(absZones)
    absCalc = 0;
    absN = 0;
else
    absCalc = 1;

    % check absZones type
    if isnumeric(absZones)
        % check absZones size
        szCur = size(absZones);
        if length(szCur)~=2 || sum(szCur>1)>1
            error('Invalid sizes')
        end
        absZones = {absZones};

    elseif iscell(absZones)
        % check absZones size
        szCur = size(absZones);
        if length(szCur)~=2 || sum(szCur>1)>1
            error('Invalid sizes')
        elseif szCur(2)>=1
            absZones = absZones';
        end
    else
        error('Invalid type')
    end
    absN = length(absZones);
end

% process combZones
if nargin<5 || isempty(combZones)
    combCalc = 0;
    combN = 0;
else
    combCalc = 1;

    % check combZones type
    if isnumeric(combZones)
        % check combZones size
        szCur = size(combZones);
        if length(szCur)~=2 || sum(szCur>1)>1
            error('Invalid sizes')
        elseif szCur(2)>=1
            combZones = combZones';
        end
        combZones = {combZones};

    elseif iscell(combZones)
        % check combZones size
        szCur = size(combZones);
        if length(szCur)~=2 || sum(szCur>1)>1
            error('Invalid sizes')
            
        end
    else
        error('Invalid type')
    end
    combN = length(combZones);
end

% process selfZones
if nargin<6 || isempty(selfZones)
    selfCalc = 0;
    selfN = 0;
else
    selfCalc = 1;

    % check selfZones type
    if iscell(selfZones)
        % check selfZones size
        szCur = size(selfZones);
        if length(szCur)~=2 || szCur(2)~=2
            error('Invalid sizes')
        end
    else
        error('Invalid type')
    end
    selfN = size(selfZones,1);
end


%%

% clear; clc
% 
% dist1 = [1:240;randi(150,1,240)];
% dist2 = [randi(35,1,240);3:242];

% distrRaw = rand(100,120);
% 
% A = {'a','b'};
% distrRaw = struct();
% for jj = 1:2
%     distrRaw.(A{jj}) = cell(2,3);
%     for ii = 1:6
%         cellNum = randi(240)+10;
%         distrRaw.(A{jj}){ii} = rand(cellNum,120);
%     end
% end
% 

% rawZones = {40:60};
% absZones = {10:40;60:70};
% combZones = {121:200};
% selfZones = {10:30,100:120};
% rawCalc = 1;
% absCalc = 1;
% combCalc = 1;
% selfCalc = 1;
% rawN = size(rawZones,1);
% absN = size(absZones,1);
% combN = size(combZones,1);
% selfN = size(selfZones,1);


% determine structure info
if isstruct(dist1)
    structBool = 1;
    fields = fieldnames(dist1);
else
    structBool = 0;
    fields = 1;
end
nFields = length(fields);

% combine final analysis parameters
allParams = struct();
allParams.calc = [rawCalc, absCalc, combCalc, selfCalc];
allParams.zones = {rawZones,absZones,combZones,selfZones};
allParams.n = [rawN,absN,combN,selfN];
allParams.names = {'Raw','Abs','Comb','Self'};
allParams.fun = {@rawDiff,@absDiff,@combCorr,@selfCorr};
nSets = length(allParams.names);


%% Quantify distribution

quantData = struct();

% loop through structs
for ss = 1:nFields
    % load current field
    if structBool
        curStrData1 = dist1.(fields{ss});
        curStrData2 = dist2.(fields{ss});
    else
        curStrData1 = dist1;
        curStrData2 = dist2;
    end

    % determine cell info
    if iscell(curStrData1)
        cellBool = 1;
        nCells = numel(curStrData1);
    else
        cellBool = 0;
        nCells = 1;
    end

    % loop through index sets
    for curSet = 1:nSets
        % define zone calculation
        if allParams.calc(curSet)==0; continue; end
        curZone = allParams.zones{curSet};
        curZoneN = allParams.n(curSet);
        curZoneName = allParams.names{curSet};
        curZoneFun = allParams.fun{curSet};

        if cellBool
            outData = cell(size(curStrData1));
        end

        % loop through cells
        for cc = 1:nCells
            %% Load base level distribution data

            % load current data
            if cellBool
                curData1 = curStrData1{cc};
                curData2 = curStrData2{cc};
            else
                curData1 = curStrData1;
                curData2 = curStrData2;
            end
            curData = cat(3,curData1,curData2);


            %% Perform calculations

            % loop through sub-calculations
            curOut = cell(curZoneN,1);
            for curSub = 1:curZoneN
                % define sub-calculation
                subZone = curZone(curSub,:);

                % perform calcuation
                subOut = curZoneFun(curData,subZone);
                curOut{curSub} = subOut;
            end

            % set output data
            if cellBool
                outData{cc} = curOut;
            else
                outData = curOut;
            end
        end

        % store output data in output struct
        if structBool
            quantData.(curZoneName).(fields{ss}) = outData;
        else
            quantData.(curZoneName) = outData;
        end
    end
end

end


%% Define analysis functions

function out = rawDiff(curData,subZone)
% raw difference calculation
out = mean(diff(curData(:,subZone{1},:),[],3),2,'omitnan');
end

% absolute difference calculation
function out = absDiff(curData,subZone)
out = mean(abs(diff(curData(:,subZone{1},:),[],3)),2,'omitnan');
end

% combination correlation calculation
function out = combCorr(curData,subZone)
A = curData(:,subZone{1},:);
out = diag(corr(A(:,:,1)',A(:,:,2)'));
end

% self correlation calculation
function out = selfCorr(curData,subZone)
A = cat(2,curData(:,:,1)',curData(:,:,2)'); 
B = diag(corr(A(subZone{1},:),A(subZone{2},:)));
out = reshape(B,[],2);
end

