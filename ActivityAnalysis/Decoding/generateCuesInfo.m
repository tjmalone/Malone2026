%% generateCuesInfo
% Generates cue and anti-cue information for FE and NE environments. Must
% be run before cue discrimination and decoding analyses. Hard coded
% generation of .mat files containing the position; side and identity
% information for cues in the 6m active learning environments.
%

clear; close all; clc

p1 = 'D:\AD_Project\imagingData\analysis_Decoding';
cd(p1)


%% Generate FE cue info

centers = [75; 100; 215; 315; 505; 540];
radii = [5; 5; 10; 10; 5; 10];
idents = [1; 2; 1; 3; 2; 3];
sides = [1; -1; -1; 1; -1; 1];
shapes = [1; 1; 2; 2; 1; 2];

rewLoc = [240 290];

cueInfo = [centers,radii,idents,sides,shapes];
key = {'centers','radii','idents','sides','shapes'};

save('cueInfo_FE.mat','cueInfo','key','rewLoc')


%% Generate NE cue info

centers = [75; 135; 280; 370; 515; 550];
radii = [5; 5; 10; 5; 10; 5];
idents = [1; 2; 1; 3; 1; 2];
sides = [-1; 1; 1; -1; -1; 1];
shapes = [1; 1; 2; 1; 2; 1];

rewLoc = [510 560];

cueInfo = [centers,radii,idents,sides,shapes];
key = {'centers','radii','idents','sides','shapes'};

save('cueInfo_NE.mat','cueInfo','key','rewLoc')


%% Generate NE anti-cue info

% define cue/reward info
cueFolder = 'D:\AnalysisCode\PostAnalysis\Cues\6mEnv2\';

binWidth = 5;

% load cue templates
load([cueFolder 'tempL.mat'])
load([cueFolder 'tempR.mat'])
temps = {tempL,tempR};

allStarts = [];
allEnds = [];
for ii = 1:2
    % Find the start and end indices of each cluster
    curStart = strfind([0, temps{ii}'], [0 1]);
    curEnd = strfind([temps{ii}', 0], [1 0]);

    allStarts = [allStarts,curStart];
    allEnds = [allEnds,curEnd];
end
nCues = length(allStarts);

% find cue bins
cueBins = [];
for ii = 1:nCues
    cuehalfwidth = round(0.5*length(allStarts(ii):(allEnds(ii))));
    cueBins  = [cueBins, (allStarts(ii)-cuehalfwidth):(allEnds(ii)+cuehalfwidth)];
end
cueBins = unique(cueBins);
antiCueIdx = setdiff(1:length(tempR),cueBins);
% antiCueIdx = 1:length(tempR);
% set random seed for reproducibility
rng(42)

% set anti-cue parameters
antiWidths = [7,11];
antiProbs = [2/3 1/3];

nCuePairs = 200;
antiCuePairs = cell(nCuePairs,2);
antiCuePairInfo = zeros(nCuePairs,4);

plotON = 1;
if plotON
    figure; hold on
    colors = 'mg';
end


% figure; hold on
for ii = 1:nCuePairs
    q = 0;
    while q<1000
        q = q+1;

        % select current cue range
        curStart = randsample(antiCueIdx,2);
        curWidth = randsample(antiWidths,2,true,antiProbs);

        curRange = cell(1,2);
        flag = 0;
        for jj = 1:2
            curRange{jj} = curStart(jj):curStart(jj)+curWidth(jj)-1;

            % check for border/cue intersects
            if any(~ismember(curRange{jj},antiCueIdx))
                flag = 1;
                continue
            end
        end

        % check for intersects
        if flag==1; continue; end
        if any(ismember(curRange{1},curRange{2})); continue; end

        % store anti-cue info
        for jj = 1:2
            antiCuePairs{ii,jj} = curRange{jj};
        end
        antiBCenters = cellfun(@mean,curRange)+0.5;
        antiBDistance = abs(diff(antiBCenters));
        antiCmDistance = antiBDistance*binWidth;
        antiCuePairInfo(ii,:) = [antiBDistance antiCmDistance antiBCenters];

        break
    end

    if plotON
        for jj = 1:2
            plot(curRange{jj},ii*ones(size(curRange{jj})),colors(jj));
        end
    end


    if q>1000
        error('No convergence')
    end
end

if plotON
    % define cue/reward info
    colorsCue = [0 0 0;0.5 0.5 0.5];
    colorsRew = [0.75 1 1];
    cueFolder = 'D:\AnalysisCode\PostAnalysis\Cues\6mEnv2\';
    rewLoc = [510 560];

    % load cue templates
    load([cueFolder 'tempL.mat'])
    load([cueFolder 'tempR.mat'])
    cueTemp = [tempL,tempR];
    cueX = 1:length(cueTemp);

    % plot cues
    cueLvl = nCuePairs;
    for cc = 1:size(cueTemp,2)
        plotCues(cueX,cueTemp(:,cc),cueLvl,colorsCue(cc,:));
    end

    xlim([0 120])
    ylim([0.5 nCuePairs+0.5])
end

antiKey = {'bin distance','cm distance','binCenter1','binCeneter2'};
save('antiCuePairInfo_NE.mat','antiCuePairInfo','antiCuePairs','antiKey','rewLoc')


%% Generate pairwise cue info

sfx = {'FE','NE'};
cueFolders = {'D:\AnalysisCode\PostAnalysis\Cues\6mEnv1\',...
    'D:\AnalysisCode\PostAnalysis\Cues\6mEnv2\'};
binWidth = 5;

for ff = 1:2
    %% Calculate cue pair info

    load(['cueInfo_' sfx{ff} '.mat'],'cueInfo')

    % define cue pairs
    nCues = size(cueInfo,1);
    [x,y] = ndgrid(1:nCues,1:nCues);
    cuePairIdxs = [x(:) y(:)];

    % initialize data matrix
    cuePairInfo = zeros(size(cuePairIdxs,1),9);
    cuePairs = cell(size(cuePairIdxs,1),2);
    delPairs = [];

    % load cue templates
    load([cueFolders{ff} 'tempL.mat'])
    load([cueFolders{ff} 'tempR.mat'])
    temps = {tempL,tempR};

    % define cue bins
    allStarts = [];
    allEnds = [];
    for ii = 1:2
        % find the start and end indices of each cluster
        curStart = strfind([0, temps{ii}'], [0 1]);
        curEnd = strfind([temps{ii}', 0], [1 0]);

        % store start and end indices
        allStarts = [allStarts,curStart];
        allEnds = [allEnds,curEnd];
    end


    for ii = 1:size(cuePairIdxs,1)
        % identify cues
        cue1 = cuePairIdxs(ii,1);
        cue2 = cuePairIdxs(ii,2);

        % skip same cues and duplicate cue pairs
        if cue1>=cue2
            delPairs(end+1) = ii;
            continue
        end

        % calculate pair info
        cmDistance = abs(diff(cueInfo([cue1 cue2],1)));
        bDistance = cmDistance/binWidth;
        bCenter1 = binCenter(cueInfo(cue1,1)-cueInfo(cue1,2),cueInfo(cue1,1)+cueInfo(cue1,2),binWidth);
        bCenter2 = binCenter(cueInfo(cue2,1)-cueInfo(cue2,2),cueInfo(cue2,1)+cueInfo(cue2,2),binWidth);
        sameIdent = cueInfo(cue1,3)==cueInfo(cue2,3);
        sameSide = cueInfo(cue1,4)==cueInfo(cue2,4);
        sameShape = cueInfo(cue1,5)==cueInfo(cue2,5);

        % store pair info
        cuePairInfo(ii,:) = [cue1, cue2, bDistance, cmDistance, bCenter1, bCenter2, sameIdent, sameSide, sameShape];

        % calculate cue bins
        cuehalfwidth = round(0.5*length(allStarts(cue1):(allEnds(cue1))));
        cuePairs{ii,1}  = (allStarts(cue1)-cuehalfwidth):(allEnds(cue1)+cuehalfwidth);
        cuehalfwidth = round(0.5*length(allStarts(cue2):(allEnds(cue2))));
        cuePairs{ii,2}  = (allStarts(cue2)-cuehalfwidth):(allEnds(cue2)+cuehalfwidth);

    end

    % delete same cues
    cuePairInfo(delPairs,:) = [];
    cuePairs(delPairs,:) = [];


    %% Save cue pair info

    % save cue pair info
    key = {'cue1','cue2','bin distance','cm distance','binCenter1','binCeneter2','sameIdent','sameSide','sameShape'};
    save(['cuePairInfo_'  sfx{ff} '.mat'],'cuePairInfo','cuePairs','key')
end

