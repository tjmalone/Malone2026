%% speedDecoding
% Calculate speed decoding for all FOV and plot data separated by funtional
% cell type, morphology, and sex. Initial implementation uses all cells,
% not just common cells, of all types.
%

clear; clc; close all;

p1 = 'D:\AD_Project\imagingData';
cd(p1)

% load alignment data
load('foldersLearning.mat')
nFOV = length(alignsLearning);
useDays = 2:size(alignsLearning{1},2);
nDays = length(useDays);

% set decoding settings (bin size must be divisible by speed range)
speedMin = 5;
speedMax = 50;
speedStep = 5;
nBins = (speedMax-speedMin)/speedStep;
binErrorThresh = 1; % correct decoding bin threshold
fileBase = ['decodingDataSp_' num2str(speedMin) '-' num2str(speedStep) '-'...
    num2str(speedMax) '_' num2str(binErrorThresh) '_'];

trackLength = 600;
nSubCell = 30;  % number of cels per iteration
nItr = 100;     % number of iterations

% load genotype and sex information
load('groupIDs.mat')
useSex = 1:3;
nSex = length(useSex);
nGroups = length(groups);

% load cell type information
load('data\globalCellTypes.mat','globalCellTypeLogical')
fldNamesT = fieldnames(globalCellTypeLogical);
useFieldsT = 6:10;
nType = length(useFieldsT);

svFolder = [p1 '\Figures\Figures' char(datetime('today','format','yyyyMMdd')) '\decoding'];
if ~isfolder(svFolder)
    mkdir(svFolder)
end

% set cell type categories
% cellTypes = {'all','common','grid','cue'};
cellTypes = {'all','common'};
nCellType = length(cellTypes);

% set random seed
rng(42)

% whether to rerun activity extraction
newRun = 1;


%% Extract activity data

% loop through cell types
for cType = 1:nCellType
    curCellType = cellTypes{cType};


    if isfile(['data\decoding\' fileBase curCellType '.mat']) && newRun==0
        load(['data\decoding\' fileBase curCellType '.mat'])
    else
        decodingDataSp = struct();
        decodingDataSp.meanPosError = zeros(nFOV,nDays);
        decodingDataSp.meanPerCorrect = zeros(nFOV,nDays);
        decodingDataSp.meanPosErrorSpeed = zeros(nFOV,nDays,nBins);
        decodingDataSp.meanPerCorrectSpeed = zeros(nFOV,nDays,nBins);
        decodingDataSp.complete = false(nFOV,nDays);
    end

    % loop through FOV
    for ii = 1:nFOV

        % loop through days
        for jj = 1:nDays
            %% Calculate cue identity discrimination

            if decodingDataSp.complete(ii,jj)==1; continue; end

            cd(foldersLearning{ii}{useDays(jj)})
            disp([num2str(ii) '-' num2str(jj)])

            % load activity information
            load('abf.mat')
            d = dir('dfof_*.mat');
            load(d(1).name,'dfof')

            % select cells
            load('localCellTypes_sig.mat','localTypeMat')
            nCells = size(localTypeMat,1);

            % define use cells
            if strcmp(curCellType,'all')
                useCells = 1:nCells;
            elseif strcmp(curCellType,'common')
                useCells = alignsLearning{ii}(:,useDays(jj));
            elseif strcmp(curCellType,'grid')
                useCells = find(localTypeMat(:,end)==1);
            elseif strcmp(curCellType,'cue')
                useCells = find(max(localTypeMat(:,1:end-1),[],2)==1);
            end

            % check cell number
            flag = false;
            if length(useCells)<nSubCell
                flag = true;
            end

            if ~flag
                % update dfof
                dfofUse = dfof(:,useCells);

                % calculate decoding with all y values
                 decdPos = decodingBySpeed_TM(dfofUse,abf,...
                    speedMin,speedMax,speedStep,binErrorThresh,nItr,nSubCell);
                meanPosErrorAll = decdPos.meanPosError;
                perCorrect = decdPos.perCorrect;
                errorPosSpeed = decdPos.errorPosSpeed;
                perCorrectSpeed = decdPos.perCorrectSpeed;

                % calculat data means
                meanPosError = mean(meanPosErrorAll,'omitnan');
                meanPerCorrect = mean(perCorrect,'omitnan');
                meanPosErrorSpeed = mean(errorPosSpeed,1,'omitnan');
                meanPerCorrectSpeed = mean(perCorrectSpeed,1,'omitnan');
            else
                meanPosError = NaN;
                meanPerCorrect = NaN;
                meanPosErrorSpeed = nan(1,1,nBins);
                meanPerCorrectSpeed = nan(1,1,nBins);
            end

            % % add histogram or position erros
            % figure
            % posErrors = decdPos.errorBin;
            % histogram(posErrors)
            % title([num2str(ii) '-' num2str(jj)])

            % identify true day index
            trueIdx = trueDays(ii,useDays(jj));
            if trueIdx>useDays(jj)
                % set skipped day data to nan
                if trueDays(ii,useDays(jj)-1)==trueIdx-2
                    decodingDataSp.meanPosError(ii,jj) = NaN;
                    decodingDataSp.meanPerCorrect(ii,jj) = NaN;
                    decodingDataSp.meanPosErrorSpeed(ii,jj,:) = nan(1,1,nBins);
                    decodingDataSp.meanPerCorrectSpeed(ii,jj,:) = nan(1,1,nBins);
                end

                % skip days past day limits
                if trueIdx>max(useDays)
                    decodingDataSp.complete(ii,jj) = true;
                    continue
                end
            end

            % store daily activity matrices
            storeIdx = trueIdx-useDays(1)+1;
            decodingDataSp.meanPosError(ii,storeIdx) = meanPosError;
            decodingDataSp.meanPerCorrect(ii,storeIdx) = meanPerCorrect;
            decodingDataSp.meanPosErrorSpeed(ii,storeIdx,:) = meanPosErrorSpeed;
            decodingDataSp.meanPerCorrectSpeed(ii,storeIdx,:) = meanPerCorrectSpeed;

            % save current data
            decodingDataSp.complete(ii,jj) = true;
            save([p1 '\data\' fileBase curCellType '.mat'],'decodingDataSp')

        end
    end

    cd(p1)

    save(['data\decoding\' fileBase curCellType '.mat'],'decodingDataSp')


    %% Concatenate activity data

    decodingDataGroupSp = struct();
    fldNames = fieldnames(decodingDataSp);

    % group decoding information by genotype group
    for ff = 1:length(fldNames)-1
        for ss = 1:length(sexes)
            decodingDataGroupSp.(fldNames{ff}).(sexIDs{ss}) = cell(nGroups,nDays);
            for gg = 1:nGroups
                curMice = intersect(groups{gg},sexes{ss});
                for jj = 1:nDays
                    decodingDataGroupSp.(fldNames{ff}).(sexIDs{ss}){gg,jj} =...
                        squeeze(decodingDataSp.(fldNames{ff})(curMice,jj,:));
                end
            end
        end
    end

    save(['data\decoding\' fileBase curCellType '.mat'],'decodingDataSp','decodingDataGroupSp')


    %% Plot overall decoding

    close all

    % define pair details
    colors = {[0 0 1],[1 0 0]};
    legs = {'WT','AD'};
    useDayCats = {1,7:10};
    nCats = length(useDayCats);
    Xlab = {'Pre-Learning','Post-Learning'};

    fldPlot = {'meanPerCorrect','meanPosError'};
    fldPlotSpeed = {'meanPerCorrectSpeed','meanPosErrorSpeed'};
    fldTitles = {'Percent Correct Speed Bin','Speed Decoding Error'};
    yLabs = {'% correct bins','cm/s'};

    ttlSuffix = [': ' num2str(speedMin) 'cm/s - ' num2str(speedMax) 'cm/s'];
    xLimits = [0.5 nDays+0.5];
    % yLimits = [0.1 0.2];
    nFields = length(fldPlot);

    useSexes = 1;

    for ss = 1:length(useSexes)
        curSex = useSexes(ss);
        curSexID = sexIDs{curSex};

        h1 = figure;
        h2 = figure;
        for ff = 1:nFields
            %% Plot by day

            % current data
            curDataCat = cell(nGroups,1);
            curData = decodingDataGroupSp.(fldPlot{ff}).(curSexID);
            for gg = 1:nGroups
                curDataCat{gg} = cat(2,curData{gg,:});
            end

            % calculate anova p values with multiple comparisons
            [pACur,pMC] = anovaRM2W_full_BH(curDataCat{1},curDataCat{2},1);
            pAnova = pACur([1 3])'.*ones(nDays,2);

            % plot data
            figure(h1)
            subplot(1,nFields,ff); hold on
            plotErrorSig(1:nDays,curDataCat{1},curDataCat{2},legs,pAnova,pMC,colors)

            % set plot labels
            xlabel('NE day')
            ylabel(yLabs{ff})
            xlim(xLimits)
            % ylim(yLimits)
            title(fldTitles{ff})


            %% Plot before versus after learning

            % calculate statistics for category days
            catData = cell(nCats,nGroups);
            for cc = 1:nCats
                for gg = 1:nGroups
                    catData{cc,gg} = mean(curDataCat{gg}(:,useDayCats{cc}),2,'omitnan');
                end
            end

            % plot bar graph
            figure(h2)
            subplot(1,nFields,ff); hold on
            pShow = [1 2;1 3;2 4;3 4];
            barGroup(Xlab,catData,colors,pShow);

            % set plot labels
            title(fldTitles{ff})
            ylabel(yLabs{ff})
            % ylim([0 yLimits(2)])
            set(gca,'FontSize',12)

        end
        figure(h1)
        sgtitle(['Overall Speed Decoding: ' curCellType ' cells, ' curSexID ' mice (nSub=30)'])
        savefig([svFolder '\' fileBase '_All_' curCellType '_' curSexID])

        figure(h2)
        sgtitle(['Overall Speed Decoding: ' curCellType ' cells, ' curSexID ' mice (nSub=30)'])
        savefig([svFolder '\' fileBase '_All-bar_' curCellType '_' curSexID])


        %% Plot speed decoding

        figure; hold on
        useDayCats = {1,7:10};
        speedX = speedMin+speedStep/2:speedStep:speedMax-speedStep/2;

        for ff = 1:nFields
            for dd = 1:length(useDayCats)
                subplot(length(useDayCats),nFields,(dd-1)*nFields + ff); hold on
                curDays = useDayCats{dd};

                curData = cell(nGroups,1);
                for gg = 1:nGroups
                    tempData = cat(3,decodingDataGroupSp.(fldPlotSpeed{ff}).(curSexID){gg,curDays});
                    curData{gg} = mean(tempData,3,'omitnan');
                end

                % remove nan FOV and columns
                curData = cellfun(@(x) x(~all(isnan(x),2),:),curData,'UniformOutput',false);
                nc = cellfun(@(x) any(isnan(x),1),curData,'UniformOutput',false);
                nanCol = find(~any(cat(1,nc{:}),1));

                % calculate anova p values with multiple comparisons
                pAnova = nan(nBins,2);
                pMC = nan(nBins,1);
                [pACur,pMCCur] = anovaRM2W_full_BH(curData{1}(:,nanCol),curData{2}(:,nanCol),0);
                pAnova(nanCol,:) = pACur([1 3])'.*ones(length(nanCol),2);
                pMC(nanCol) = pMCCur;
                pAnova = pAnova(nanCol);
                pMC = pMC(nanCol);

                % plot data
                plotErrorSig(speedX(nanCol),curData{1}(:,nanCol),curData{2}(:,nanCol),legs,pAnova,pMC,colors)

                % set labels
                xlabel('Speed (cm/s)')
                ylabel(yLabs{ff})
                title([fldTitles{ff} ': Day ' num2str(curDays)])
                set(gca,'fontsize',12)
            end

            sgtitle(['Speed Decoding: ' curCellType ' cells (' curSexID ' mice)' ])
            savefig([svFolder '\' fileBase '_Speed_' curCellType '_' curSexID])
        end
    end
end
