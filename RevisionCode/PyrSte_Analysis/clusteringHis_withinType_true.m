%% Set inputs

clear; close all; clc

% set base folder
p1 = 'Z:\labMembers\TM\AD_Project\Revision\Histology\revisionHistology';
cd(p1)

% find all validation selections
foldsHV = findSubF('HV',2,[],0);

% define use FOV
useFOV = [3:6 10:12 16];

nFOV = length(useFOV);

% load prior cluster data
load('D:\AD_Project\Revision\PyrSteAnalysis\dataCluster.mat','dataCluster')
oldDataCluster = dataCluster;
oldFOV = [11:12 14:17];
oldFOVMap = [2:3 5:8];

% neighbors (N)
N = 1:25;
nNeigh = length(N);

% shuffles
rng(42)
nShuff = 1000;
shuffThresh = [5 50 95];
shuffColor = [0.5 0.5 0.5];

pixelSize = 850.19/1024;    % conversion from pixels to distance

% random cell selection reps
nReps = 100;

% cell identities
IDs = [-1 1];
nID = length(IDs);

% set cell type info
cellTypes = {'calbindin','reelin'};


%% Proccess ROIs

dataCluster = struct();
dataCluster.centroids = cell(nFOV,1);
dataCluster.ident = cell(nFOV,1);
dataCluster.IDDist = cell(nFOV,1);
dataCluster.IDNeigh = cell(nFOV,1);

plotDataGlob = cell(nID,nFOV);

for ff = 1:nFOV
    %% Load and process cells

    % load current roi data
    curF = useFOV(ff);
    load([foldsHV{curF} '\roiData.mat'],'roiData')

    if ismember(curF,oldFOV)
        % update with pre-saved info
        useOld = oldFOVMap(oldFOV==curF);
        centroidsAll = oldDataCluster.centroids{useOld}*pixelSize;
        identsAll = oldDataCluster.ident{useOld};
    else
        % concatenate data
        identsAll = roiData.identity;
        centroidsAll  = roiData.centroid*pixelSize;
    end
    nROI = size(identsAll,1);


    %% Calculate distance

    % shuffle identity
    identShuff = zeros(nROI,nShuff);
    for ss = 1:nShuff
        identShuff(:,ss) = randsample(identsAll,nROI);
    end

    % calculate all distances
    [distsAll,distsIdx] = pdist2(centroidsAll,centroidsAll,...
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
        curMatchIdx = find(identsAll(distsIdx(:,rr))==identsAll(rr));
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
        curCells = identsAll==IDs(ii);
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

    dataCluster.centroids{ff} = centroidsAll;
    dataCluster.ident{ff} = identsAll;
    dataCluster.distNeigh{ff} = distNeigh;
    dataCluster.distNeighShuff{ff} = distNeighShuff;
end

save('Z:\labMembers\TM\AD_Project\Revision\Histology\histologyValidation\histClusteringType_trueB.mat','dataCluster')


%% Plot combined data

clear; close all; clc
load('Z:\labMembers\TM\AD_Project\Revision\Histology\histologyValidation\histClusteringType_trueB.mat','dataCluster')
nFOV = length(dataCluster.ident);

% define use data
N = 1:25;
nNeigh = length(N);
useIdx = 25;
useN = N(useIdx);

% define plot settings
kX = 350:650;
nX = length(kX);
trueColors = {'b','r'};
yPos = 0.025;

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
        kYShuff(ff,:) = ksdensity(plotShuff(ff,:),kX);
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
ylim([0 0.03])


%% Calculate statistics
    
statData = cell(1,2);
for ii = 1:nID
    statData{ii} = cat(2,dataCluster.meanReal{:,ii})';
end

% find percentiles per FOV
shuffTile = cell(1,nID);
for ii = 1:nID
    shuffTile{ii} = zeros(nFOV,nNeigh);
    for ff = 1:nFOV
        for nn = 1:nNeigh
            curShuff = dataCluster.meanShuff{ff,ii}(N(nn),:)';
            curReal = dataCluster.meanReal{ff,ii}(N(nn));
            shuffTile{ii}(ff,nn) = prctileInv(curShuff,curReal);
        end
    end
end
% statData = shuffTile;
meanTile = cellfun(@(x) mean(x(:,useIdx)),shuffTile);

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
