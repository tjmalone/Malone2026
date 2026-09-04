%% Load folders/file infomation

clear; close all; clc

% set directory
foldMan = 'D:\AD_Project\Revision\PyrSteAnalysis\';
cd(foldMan);

% load projections
load('data_PyrSte_rev.mat','dataPS');
nFOV = length(dataPS.allThresh);

% threshold for morphology range
threshRange = 10;

% pixel resolution
pr = 750/512;

% distance radii (D)
D = 5:5:100;
nDists = length(D);

% neighbors (N)
N = 40:41;
nNeigh = length(N);

combX = {D,N};
combXLab = {'Distance (um)','Nearest neighbors'};

% shuffles
rng(42)
nShuff = 10;
shuffThresh = [5 50 95];
shuffColor = [0.5 0.5 0.5];

% random cell selection reps
rng(42)
nReps = 10;

% cell identities
IDs = [-1 1];
nID = length(IDs);

% set cell type info
cellTypes = {'calbindin','reelin'};


%% Calculate centroids

dataCluster = struct();
dataCluster.centroids = cell(nFOV,1);
dataCluster.ident = cell(nFOV,1);
dataCluster.IDDist = cell(nFOV,1);
dataCluster.IDNeigh = cell(nFOV,1);

plotDataGlob = cell(2,nID,nFOV);

% loop through FOV
for ff = 1:nFOV
    %% Calculate cell identity
    disp(ff)

    % threshold for morphology center
    threshCenter = dataPS.allThresh(ff);
    if isnan(threshCenter)
        continue
    end

    % current morphology edges
    thLow = threshCenter-threshRange;
    thHigh = threshCenter+threshRange;

    % calculate morphology classes
    ident = zeros(size(dataPS.allArea{ff},1),1);
    ident(dataPS.allArea{ff}<=thLow) = -1;
    ident(dataPS.allArea{ff}>=thHigh) = 1;


    %% Calculate centroids

    % get current manual footprints
    curROIsManual = dataPS.allRoiMan{ff};
    nROI = size(curROIsManual,3);

    % calculate properties of manual rois
    centroidMan = zeros(size(curROIsManual,3),2);
    for kk = 1:size(curROIsManual,3)
        rProps = regionprops(curROIsManual(:,:,kk),'centroid','PixelIdxList');
        centroidMan(kk,:) = rProps(1).Centroid;
    end


    %% Calculate distance

    % shuffle identity
    identShuff = zeros(nROI,nShuff);
    for ss = 1:nShuff
        identShuff(:,ss) = randsample(ident,nROI);
    end

    % calculate all distances
    [distsAll,distsIdx] = pdist2(centroidMan,centroidMan,...
        'euclidean','Smallest',nROI);

    % precompute shuffle permutation
    repIdx = cell(nROI,1);
    for n = 1:nROI
        [~,repIdx{n}] = sort(rand(n,nReps),1);
    end

    % initialize identity matrics
    distNeigh = zeros(nROI,nNeigh);
    distNeighShuff = zeros(nROI,nNeigh,nShuff);

    for rr = 1:nROI
        % find nearest same type neighbors
        curMatchDists = nan(nROI,1);
        curMatchIdx = find(ident(distsIdx(:,rr))==ident(rr));
        curMatchDists(1:length(curMatchIdx)-1) = distsAll(curMatchIdx(2:end),rr);

        % find nearest same type neighbors for shuffle
        curShuffDists = nan(nROI,nShuff);
        for ss = 1:nShuff
            curShuffIdx = find(identShuff(distsIdx(:,rr),ss)==identShuff(rr,ss));
            curShuffDists(1:length(curShuffIdx)-1,ss) = distsAll(curShuffIdx(2:end),rr);
        end

        % generate random cell order
        nRReal = find(~isnan(curMatchDists),1,'last');

        % calculate average distance of same type identity
        for nn = 1:nNeigh
            % find cells within N neighbors
            curN = N(nn);
            curRepIdxReal = repIdx{nRReal}(1:min(curN,nRReal),:);

            % calculate true average
            distNeigh(rr,nn) = mean(curMatchDists(curRepIdxReal),[1 2],'omitnan');

            for ss = 1:nShuff
                nRShuff = find(~isnan(curShuffDists(:,ss)),1,'last');

                curRepIdxShuff = repIdx{nRShuff}(1:min(curN,nRShuff),:);

                distNeighShuff(rr,nn,ss) = mean(curShuffDists(curRepIdxShuff,ss),'all','omitnan');
            end
        end
    end

    %% Calculate identity averages

    figure; tiledlayout(1,nID)

    for ii = 1:nID
        % find cells with current identity
        curCells = ident==IDs(ii);
        curCellsShuff = identShuff==(IDs(ii));

        % find mean identity
        catReal = distNeigh(curCells,:);
        meanReal = mean(catReal,1,'omitnan')';

        % find mean shuffle
        catShuff = zeros(sum(curCells),nNeigh,nShuff);
        for ss = 1:nShuff
            catShuff(:,:,ss) = distNeighShuff(curCellsShuff(:,ss),:,ss);
        end
        meanShuff = squeeze(mean(catShuff,1,'omitnan'));

        % find shuffle percentiles
        shuffTile = prctile(meanShuff,shuffThresh,2);

        % plot identity
        nexttile(ii); hold on
        plot(shuffTile,'Color',shuffColor)
        plot(meanReal,'k','LineWidth',2)
        % ylim([0 100])
        title(cellTypes{ii})

        % store data
        dataCluster.meanReal{ff,ii} = meanReal;
        dataCluster.meanShuff{ff,ii} = meanShuff;
    end

    close


    %% Store results

    dataCluster.centroids{ff} = centroidMan;
    dataCluster.ident{ff} = ident;
    dataCluster.distNeigh{ff} = distNeigh;
    dataCluster.distNeighShuff{ff} = distNeighShuff;
end

save('Z:\labMembers\TM\AD_Project\Revision\Histology\histologyValidation\ManClusteringType.mat','dataCluster')


%% Plot combined data

clear; close all; clc
load('Z:\labMembers\TM\AD_Project\Revision\Histology\histologyValidation\ManClusteringType.mat','dataCluster')
nFOV = sum(~cellfun(@isempty,dataCluster.ident));

% define use data
N = 40:41;
nNeigh = length(N);
useIdx = 2;
useN = N(useIdx);

% define plot settings
kX = 125:275;
ksWidth = 7;
nX = length(kX);
trueColors = {'b','r'};
yPos = 0.035;

% cell identities
IDs = [-1 1];
nID = length(IDs);
legs = {'Pyr','Ste'};

figure; hold on
indData = cell(1,nID);

h = zeros(1,nID);
shuffTileGlob = zeros(1,nID);
for ii = 1:nID
    % get current data
    plotReal = cat(2,dataCluster.meanReal{:,ii});
    plotReal = plotReal(useIdx,:)';
    indData{ii} = plotReal;

    % get current shuffle
    plotShuff = cat(3,dataCluster.meanShuff{:,ii});
    plotShuff = squeeze(plotShuff(useIdx,:,:))';

    % plot shuffles
    kYShuff = zeros(nFOV,nX);
    for ff = 1:nFOV
        kYShuff(ff,:) = ksdensity(plotShuff(ff,:),kX,'Width',ksWidth);
    end
    semshade(kYShuff,0.3,trueColors{ii},kX);

    % plot real data mean +/- SEM
    M = mean(plotReal,'omitnan');
    E = nansem(plotReal,1);
    h(ii) = errorbar(M,yPos,E,'horizontal','Color',trueColors{ii},...
        'LineWidth',1.5);
    xline(mean(plotReal),'Color',trueColors{ii},'LineWidth',1.5)

    % find global percentile
    meanShuff = mean(kYShuff,1);
    shuffTileGlob(ii) = sum(meanShuff(kX<M))*100;
end

% add labels
title(['Inter-Cell Distance: (' num2str(useN) ' cells)']);
xlabel('Distance (\mum)')
ylabel('Shuffle percentile')
legend(h,legs)
ylim([0 0.04])
xlim([kX(1) kX(end)])
% savefig('Results\Clustering\clusteringManual.fig')
% save('Results\Clustering\clusteringMan.mat','dataCluster','plotDataGlob','indData')


%% Calculate statistics

statData = cell(1,2);
for ii = 1:nID
    statData{ii} = cat(2,dataCluster.meanReal{:,ii})';
end

% find percentiles per FOV
nFOVAll = size(dataCluster.meanShuff,1);
shuffTile = zeros(nFOVAll,ii);
for ff = 1:nFOVAll
    for ii = 1:nID
        if ~isempty(dataCluster.meanReal{ff,ii})
            curShuff = dataCluster.meanShuff{ff,ii}(useIdx,:)';
            curReal = dataCluster.meanReal{ff,ii}(useIdx);
            shuffTile(ff,ii) = prctileInv(curShuff,curReal);
        else
            shuffTile(ff,ii) = NaN;
        end
    end
end
shuffTileMean = mean(shuffTile,1,'omitnan');

% set test parameters
nUnits = ' FOVs';
testName = 'two-tailed paired Students t-test';
testPair = 1;
testMC = 0;
testLimitP = 0;

% perfrom ttest for all compairsons
outStats = {};
outStats(end+1,:) = ['' ttestEffectSize(statData{1},statData{2},testName,nUnits,testPair,testMC,testLimitP)];

% perfrom ttest for use compairson
outStats(end+1,:) = ['' ttestEffectSize(statData{1}(:,useIdx),statData{2}(:,useIdx),testName,nUnits,testPair,testMC,testLimitP)];


%%

szs = zeros(nFOV,nID);
for ff = 1:nFOV
    curIdents = dataCluster.ident{ff};
    for ii = 1:nID
        szs(ff,ii) = sum(curIdents==IDs(ii));
    end
end

