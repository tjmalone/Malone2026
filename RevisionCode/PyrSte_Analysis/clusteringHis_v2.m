%% Set inputs

clear; close all; clc

% set base folder
p1 = 'Z:\labMembers\TM\AD_Project\Revision\Histology\revisionHistology';
cd(p1)

% find all validation selections
foldsHV = findSubF('HV',2,[],0);

% define use FOV
useFOV = [3:6 11:12 16];

nFOV = length(useFOV);

% load prior cluster data
load('D:\AD_Project\Revision\PyrSteAnalysis\dataCluster.mat','dataCluster')
oldDataCluster = dataCluster;
oldFOV = [11:12 14:17];
oldFOVMap = [2:3 5:8];

% load putative cell identities
load('Z:\labMembers\TM\AD_Project\Revision\Histology\histologyValidation\histValAccuracy.mat','dataAcc')
putIdents = dataAcc.ind.ident;

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
IDs = [-1 1];
nID = length(IDs);

% set cell type info
cellTypes = {'calbindin','reelin'};


%% Proccess ROIs

identTypes = {'True','Putative'};

% loop through identity types
for jj = 1:2
    dataCluster = struct();
    dataCluster.centroids = cell(nFOV,1);
    dataCluster.ident = cell(nFOV,1);
    dataCluster.IDDist = cell(nFOV,1);
    dataCluster.IDNeigh = cell(nFOV,1);

    plotDataGlob = cell(2,nID,nFOV);

    for ff = 1:nFOV
        %% Load and process cells

        % load current roi data
        curF = useFOV(ff);
        load([foldsHV{curF} '\roiData.mat'],'roiData')

        if ismember(curF,oldFOV)
            % update with pre-saved info
            useOld = oldFOVMap(oldFOV==curF);
            centroidsAll = oldDataCluster.centroids{useOld};
            identsAll = oldDataCluster.ident{useOld};
        else
            % concatenate data
            identsAll = roiData.identity;
            centroidsAll  = roiData.centroid;
        end
        nROI = size(identsAll,1);

        % check roi match
        if length(identsAll)~=length(putIdents{ff})
            error('Cell number mismatch')
        end

        % set current identity type
        if jj==1
            curIdents = identsAll;
        else
            curIdents = putIdents{ff};
        end


        %% Calculate distance

        % shuffle identity
        identShuff = zeros(nROI,nShuff);
        for ss = 1:nShuff
            identShuff(:,ss) = randsample(curIdents,nROI);
        end

        % calculate all distances
        [distsAll,distsIdx] = pdist2(centroidsAll,centroidsAll,...
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
                    IDDist(rr,dd) = mean(curIdents(distsIdx(curIdx(2:end),rr)));
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
                IDNeigh(rr,nn) = mean(curIdents(distsIdx(2:curN+1,rr)));
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
                curCells = curIdents==IDs(ii);
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


        %% Store results

        dataCluster.centroids{ff} = centroidsAll;
        dataCluster.ident{ff} = curIdents;
        dataCluster.IDDist{ff} = IDDist;
        dataCluster.IDNeigh{ff} = IDNeigh;

        close
    end


    %% Plot combined data

    figure; tiledlayout(2,nID)
    sgtitle([identTypes{jj} ' Cell Identity'])
    indData = cell(2,nID);

    ttSupps = {'Distance','N Neighbors'};
    iiSupps = {'Pyr','Ste'};

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

    % cd(foldSave)

    % savefig('Results\Clustering\clusteringHis.fig')
    % save('Results\Clustering\clusteringHis.mat','dataCluster','plotDataGlob','indData')
end

