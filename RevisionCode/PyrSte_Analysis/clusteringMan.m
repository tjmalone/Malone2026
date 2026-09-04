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
N = 1:20;
nNeigh = length(N);

combX = {D,N};
combXLab = {'Distance (um)','Nearest neighbors'};

% shuffles
rng(42)
nShuff = 1000;
shuffThresh = [5 50 95];
shuffColor = [0.5 0.5 0.5];

% cell identities
IDs = -1:1;
nID = length(IDs);


%% Calculate centroids

dataCluster = struct();
dataCluster.centroids = cell(nFOV,1);
dataCluster.ident = cell(nFOV,1);
dataCluster.IDDist = cell(nFOV,1);
dataCluster.IDNeigh = cell(nFOV,1);

plotDataGlob = cell(2,nID,nFOV);

% loop through FOV
for ff = 8:nFOV
    %% Calculate cell identity

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

    % initialize identity matrics
    IDDist = zeros(nROI,nDists);
    IDDistShuff = zeros(nROI,nDists,nShuff);
    IDNeigh = zeros(nROI,nNeigh);
    IDNeighShuff = zeros(nROI,nNeigh,nShuff);

    for rr = 1:nROI
        % calculate average identity by distance
        for dd = 1:nDists
            % find cells within distance D
            curD = D(dd);
            curIdx = find(distsAll(:,rr)<=curD);

            % get average cell identity
            if length(curIdx)>1
                IDDist(rr,dd) = mean(ident(distsIdx(curIdx(2:end),rr)));
                IDDistShuff(rr,dd,:) = mean(identShuff(distsIdx(curIdx(2:end),rr),:),1);
            else
                IDDist(rr,dd) = NaN;
                IDDistShuff(rr,dd,:) = NaN;
            end
        end

        % calculate average identity by neighbors
        for nn = 1:nNeigh
            % find cells within N neighbors
            curN = N(nn);
            IDNeigh(rr,nn) = mean(ident(distsIdx(2:curN+1,rr)));
            IDNeighShuff(rr,nn,:) = mean(identShuff(distsIdx(2:curN+1,rr),:),1);
        end
    end

    % store data for plot
    storeReal = {IDDist,IDNeigh};
    storeShuff = {IDDistShuff,IDNeighShuff};


    %% Calculate identity averages

    figure; tiledlayout(2,nID)

    for tt = 1:2
        curReal = storeReal{tt};
        curShuff = storeShuff{tt};
        for ii = 1:nID
            % find cells with current identity
            curCells = ident==IDs(ii);
            curCellsShuff = identShuff==(IDs(ii));

            % find mean identity
            meanReal = mean(curReal(curCells,:),1,'omitnan')';
            meanRealShuff = zeros(nDists,nShuff);
            for ss = 1:nShuff
                meanRealShuff(:,ss) = mean(curShuff(curCellsShuff(:,ss),:,ss),1,'omitnan');
            end

            % find shuffle percentiles
            shuffTile = prctile(meanRealShuff,shuffThresh,2);

            % plot identity
            nexttile((tt-1)*nID+ii); hold on
            plot(shuffTile,'Color',shuffColor)
            plot(meanReal,'k','LineWidth',2)
            ylim([-1 1])

            % calculate and store inverse percentile
            tileInv = invPercentile(meanReal',meanRealShuff');
            plotDataGlob{tt,ii,ff} = tileInv;
        end
    end

    close


    %% Store results

    dataCluster.centroids{ff} = centroidMan;
    dataCluster.ident{ff} = ident;
    dataCluster.IDDist{ff} = IDDist;
    dataCluster.IDNeigh{ff} = IDNeigh;

end


%% Plot combined data

figure; tiledlayout(2,nID)
indData = cell(2,nID);

ttSupps = {'Distance','N Neighbors'};
iiSupps = {'Pyr','Other','Ste'};

for tt = 1:2
    for ii = 1:nID
        % get current data
        plotData = squeeze(plotDataGlob(tt,ii,:));
        plotData(cellfun(@isempty,plotData)) = [];

        % concatenate data
        plotDataCat = squeeze(cat(3,plotData{:}));
        indData{tt,ii} = plotDataCat';

        % plot identity
        nexttile((tt-1)*nID+ii); hold on
        title(['Average Identity of Neighbors by ' ttSupps{tt} ' (' iiSupps{ii} ' cells)'])

        M = mean(plotDataCat,2);
        E = nansem(plotDataCat,2);
        errorbar(combX{tt},M,E,'k','LineWidth',1)

        % add labels
        xlabel(combXLab{tt})
        ylabel('Shuffle percentile')
        ylim([0 100])
        for yy = 1:length(shuffThresh)
            yline(shuffThresh(yy),'Color',[0.5 0.5 0.5])
        end
    end
end

savefig('Results\Clustering\clusteringManual.fig')
save('Results\Clustering\clusteringMan.mat','dataCluster','plotDataGlob','indData')
