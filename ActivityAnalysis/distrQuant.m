function barData = distrQuant(distrRaw,cueTemp,rewLocIn,cueDil,cueExlBin,useCalcs)
%% distrQuant
% Performs quantification of an arbitrary spatial distribution.
% Specifically, calculates various single number characteristics of the
% distrbution: in-cue, out-cue, pre-reward, in-reward, post-reward,
% ridge-background ratio, no-lag correlation, lag to peak correlation.
% Outputs data as a single struct.
%
% Inputs:
%   distrRaw - input distribution data. Can be a numerical or cell array of
%       numerical arrays or struct of cell or numerical arrays. At lowest
%       level, numerical arrays must be identical in size along second
%       dimension (spatial dimension).
%   cueTemp - cue template. Must be numerical array the same size as base
%       distribution.
%   rewLocIn - reward window indices. Pre- and post-reward zones will be
%       the same size as in-reward zone.
%   cueDil - cue dilation factor
%   cueExlBin - defines exluded bins. If length is one, all bins after will
%       be excluded. If length is greater than one, all bins in array will
%       be excluded.
%
% Outputs:
%   barData - output structure. Each field represents one calculation type
%       and is the same data type/size as the input distribution data.
%

%% Process inputs

% set ridge background settings
rdgRange = 1:3;
bkgRange = 7:9;
maxDistance = 16;

% set correlation settings
maxShift = 20;

if nargin<3
    error('Distribution quantification requires 3 inputs')
end

if nargin<4 || isempty(cueDil)
    cueDil = 1;
end

% correct cue template dimension
szCue = size(cueTemp);
if length(szCue)>2 || sum(szCue>1)~=1
    error('Cue template must be a 1D vector')
end
if szCue(1)~=1
    cueTemp = cueTemp';
end

if nargin<5 || isempty(cueExlBin)
    cueExlBin = length(cueTemp);
end

% correct reward zone
if length(rewLocIn)==2
    rewLocIn = rewLocIn(1):rewLocIn(2);
end

% determine structure info
if isstruct(distrRaw)
    structBool = 1;
    fields = fieldnames(distrRaw);
else
    structBool = 0;
    fields = 1;
end
nFields = length(fields);

% find the start and end indices of each cluster
cueLen = length(cueTemp);
allStarts = strfind([0, cueTemp], [0 1]);
allEnds = strfind([cueTemp, 0], [1 0]);
allCenters = round((allStarts + allEnds+1) / 2);
nCues = length(allStarts);

% find expanded cue bins
cueBinsIn = [];
for ii = 1:nCues
    cuehalfwidth = (cueDil-1)*round(0.5*length(allStarts(ii):(allEnds(ii))));
    cueBinsIn  = [cueBinsIn, (allStarts(ii)-cuehalfwidth):(allEnds(ii)+cuehalfwidth)];
end
cueBinsIn(cueBinsIn<1 | cueBinsIn>cueLen) = [];
cueBinsOut = setdiff(1:cueLen,cueBinsIn);
binsAll = 1:cueLen;

% correct cue bins
if length(cueExlBin)==1
    cueBinsIn(cueBinsIn>cueExlBin) = [];
    cueBinsOut(cueBinsOut>cueExlBin) = [];
    binsAll(binsAll>cueExlBin) = [];
else
    cueBinsIn(ismember(cueBinsIn,cueExlBin)) = [];
    cueBinsOut(ismember(cueBinsOut,cueExlBin)) = [];
    binsAll(ismember(binsAll,cueExlBin)) = [];
end

% calculate the distance of each index to the nearest center
distancesToCenters = zeros(1,cueLen);
for ii = 1:cueLen
    distancesToCenters(ii) = min(abs(allCenters - ii));
end

% find reward zones
if any(rewLocIn<2) || any(rewLocIn>cueLen-1)
    error('Pre- and post-reward must cover at least 1 bin')
end
rewLen = length(rewLocIn);
rewLocPre = max(1,rewLocIn(1)-rewLen):rewLocIn(1)-1;
rewLocPost = rewLocIn(end)+1:min(cueLen,rewLocIn(end)+rewLen);

% set final list of indices
setIdx = {cueBinsIn,cueBinsOut,binsAll,rewLocPre,rewLocIn,rewLocPost};
setLabel = {'CueIn','CueOut','binsAll','RewPre','RewIn','RewPost','RdgBkg','CorrNoLag','LagToPeakCorr','Rdg','Bkg'};
nSets = length(setLabel);

if nargin<6 || isempty(useCalcs)
    useCalcs = 1:nSets;
end


%% Quantify distribution

barData = struct();

% loop through structs
for ss = 1:nFields
    % load current field
    if structBool
        curStrData = distrRaw.(fields{ss});
    else
        curStrData = distrRaw;
    end

    % determine cell info
    if iscell(curStrData)
        cellBool = 1;
        nCells = numel(curStrData);
    else
        cellBool = 0;
        nCells = 1;
    end

    % loop through index sets
    for curSet = useCalcs
        if cellBool
            outData = cell(size(curStrData));
        end

        % loop through cells
        for cc = 1:nCells
            % load current field
            if cellBool
                curData = curStrData{cc};
            else
                curData = curStrData;
            end
            cellNum = size(curData,1);

            if curSet<7
                % calculate data mean
                curDataOut = mean(curData(:,setIdx{curSet}),2,'omitnan');

            elseif curSet==7 || ismember(curSet,[10 11])
                % % calculate mean as a function of distance to cue
                % dataToCenters = zeros(cellNum,maxDistance);
                % for ii = 1:maxDistance
                %     if length(cueExlBin)==1
                %         dataToCenters(:,ii) = mean(curData(:,distancesToCenters==ii-1 &...
                %             1:length(distancesToCenters)<=cueExlBin),2,'omitnan');
                %     else
                %         dataToCenters(:,ii) = mean(curData(:,distancesToCenters==ii-1 &...
                %             ~ismember(1:length(distancesToCenters),cueExlBin)),2,'omitnan');
                %     end
                % end
                % 
                % % normalize curve to cue centers
                % % dataRB = dataToCenters;
                % dataRB = dataToCenters./mean(dataToCenters(:,1:2),2,'omitnan');
                % % dataRB = 0.5*log((1+dataToCenters)./(1-dataToCenters));
                % 
                % % calculate ridge background ratio
                % rdgM = mean(dataRB(:,rdgRange),2,'omitnan');
                % bkgM = mean(dataRB(:,bkgRange),2,'omitnan');
                % curDataOut = rdgM./bkgM;
                % % curDataOut = (rdgM-bkgM)./(rdgM+bkgM);
                % % curDataOut = (rdgM-bkgM)./(abs(rdgM)+abs(bkgM));
                % % curDataOut = log(rdgM./bkgM);

                if ~exist('rdgM','var') && ~exist('bkgM','var')
                    % calculate mean as a function of distance to cue
                    dataToCenters = cell(1,maxDistance);
                    % distrNorm = curData./mean(distrRaw,2,'omitnan');
                    for ii = 1:maxDistance
                        dataToCenters{:,ii} = curData(:,distancesToCenters==ii-1);
                    end

                    % calculate ridge background ratio
                    rdgM = mean(cat(2,dataToCenters{rdgRange}),2,'omitnan');
                    bkgM = mean(cat(2,dataToCenters{bkgRange}),2,'omitnan');
                end

                % store correct output
                if curSet==7
                    curDataOut = rdgM./bkgM;
                    curDataOut(curDataOut==Inf) = NaN;
                elseif curSet==10
                    curDataOut = rdgM;
                elseif curSet==11
                    curDataOut = bkgM;
                end

            elseif curSet==8
                % calculate no-lag correlation
                if ~isempty(curData)
                    curDataOut = corr(curData',cueTemp');
                else
                    curDataOut = NaN;
                end

            elseif curSet==9
                curDataOut = zeros(cellNum,1);
                for ii = 1:cellNum
                    % calculate cross correlation
                    [curCorrs,curLags] = xcorr(curData(ii,:),cueTemp,maxShift,'coeff');

                    % get lag to peak correlation
                    [~,idxMaxCorr] = max(curCorrs);
                    curDataOut(ii) = abs(curLags(idxMaxCorr));
                end
            end

            % set output data
            if cellBool
                outData{cc} = curDataOut;
            else
                outData = curDataOut;
            end
        end

        % store output data in output struct
        if structBool
            barData.(setLabel{curSet}).(fields{ss}) = outData;
        else
            barData.(setLabel{curSet}) = outData;
        end
    end
end

end

