%% spatialDecodingDiscrimination
% Calculate decoding, spatial decoding, and cue identity decoding for all
% FOV and plot data separated by funtional cell type, morphology, and sex.
% Initial implementation uses all cells, not just common cells, of all
% types.

clear; clc; close all;

p1 = 'D:\AD_Project\imagingData';
cd(p1)

% load alignment data
load('foldersLearning.mat')
nFOV = length(alignsLearning);
useDays = 2:size(alignsLearning{1},2);
nDays = length(useDays);

% set decoding settings
trackLength = 600;
binWidth = 5;
nBins = trackLength/binWidth;
nSubCell = 30;  % number of cels per iteration
nItr = 100;     % number of iterations
binErrorThresh = 4; % correct decoding bin threshold
ySubset = [150 450];

% load cue pair information
load('D:\AD_Project\imagingData\analysis_Decoding\cuePairInfo_NE.mat','cuePairInfo','cuePairs')
load('D:\AD_Project\imagingData\analysis_Decoding\antiCuePairInfo_NE.mat','antiCuePairs')
nPairs = size(cuePairInfo,1);
nAntiPairs = size(antiCuePairs,1);

% load genotype and sex information
load('groupIDs.mat')
useSex = 1:3;
nSex = length(useSex);
nGroups = length(groups);

svFolder = [p1 '\Figures\Figures' char(datetime('today','format','yyyyMMdd'))  '\decoding'];
if ~isfolder(svFolder)
    mkdir(svFolder)
end

% set cell type categories
cellTypes = {'all','common','grid','nongrid','cue','other'};
cellTypes = {'all'};
nCellType = length(cellTypes);

% set random seed
rng(42)

% whether to rerun activity extraction
newRun = 0;


%% Extract activity data

% loop through cell types
for cType = 1:nCellType
    curCellType = cellTypes{cType};


    if isfile(['data\decoding\decodingData_' curCellType '.mat']) && newRun==0
        load(['data\decoding\decodingData_' curCellType '.mat'])
    else
        decodingData = struct();
        decodingData.meanPosError = zeros(nFOV,nDays);
        decodingData.meanPerCorrect = zeros(nFOV,nDays);
        decodingData.meanPosErrorSpatial = zeros(nFOV,nDays,nBins);
        decodingData.meanPerCorrectSpatial = zeros(nFOV,nDays,nBins);
        decodingData.meanPerCorrectCues = zeros(nFOV,nDays,nPairs);
        decodingData.meanPerCorrectAntiCues = zeros(nFOV,nDays,nAntiPairs);
        decodingData.meanPosErrorYSub = zeros(nFOV,nDays);
        decodingData.meanPerCorrectYSub = zeros(nFOV,nDays);
        decodingData.meanPosErrorSpatialYSub = zeros(nFOV,nDays,nBins);
        decodingData.meanPerCorrectSpatialYSub = zeros(nFOV,nDays,nBins);
        decodingData.complete = false(nFOV,nDays);
    end

    % loop through FOV
    for ii = 1:nFOV

        % loop through days
        for jj = 1:nDays
            %% Calculate cue identity discrimination

            if decodingData.complete(ii,jj)==1; continue; end

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
            elseif strcmp(curCellType,'nongrid')
                useCells = find(localTypeMat(:,end)==0);
            elseif strcmp(curCellType,'cue')
                useCells = find(max(localTypeMat(:,1:end-1),[],2)==1);
            elseif strcmp(curCellType,'other')
                useCells = find(~any(localTypeMat,2));
            else
                error('Invalid cell type')
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
                [decdPos,decdCue] = decodingByLocation_TM(dfofUse,abf,...
                    trackLength,binWidth,binErrorThresh,nItr,nSubCell,{cuePairs,antiCuePairs});
                meanPosErrorAll = decdPos.meanPosError;
                perCorrect = decdPos.perCorrect;
                errorPosSpatial = decdPos.errorPosSpatial;
                perCorrectSpatial = decdPos.perCorrectSpatial;
                perCorrectCues = decdCue(1).perCorrectAllCue;
                perCorrectAntiCues = decdCue(2).perCorrectAllCue;

                % calculate decoding without ends of track
                [decdPos,decdCue] = decodingByLocation_TM(dfofUse,abf,...
                    trackLength,binWidth,binErrorThresh,nItr,nSubCell,[],ySubset);
                meanPosErrorYSubAll = decdPos.meanPosError;
                perCorrectYSub = decdPos.perCorrect;
                errorPosSpatialYSub = decdPos.errorPosSpatial;
                perCorrectSpatialYSub = decdPos.perCorrectSpatial;

                % calculat data means
                meanPosError = mean(meanPosErrorAll,'omitnan');
                meanPerCorrect = mean(perCorrect,'omitnan');
                meanPosErrorSpatial = mean(errorPosSpatial,1,'omitnan');
                meanPerCorrectSpatial = mean(perCorrectSpatial,1,'omitnan');
                meanPerCorrectCues = mean(perCorrectCues,1,'omitnan');
                meanPerCorrectAntiCues = mean(perCorrectAntiCues,1,'omitnan');
                meanPosErrorYSub = mean(meanPosErrorYSubAll,'omitnan');
                meanPerCorrectYSub = mean(perCorrectYSub,'omitnan');
                meanPosErrorSpatialYSub = mean(errorPosSpatialYSub,1,'omitnan');
                meanPerCorrectSpatialYSub = mean(perCorrectSpatialYSub,1,'omitnan');
            else
                meanPosError = NaN;
                meanPerCorrect = NaN;
                meanPosErrorSpatial = nan(1,1,nBins);
                meanPerCorrectSpatial = nan(1,1,nBins);
                meanPerCorrectCues = nan(1,1,nPairs);
                meanPerCorrectAntiCues = nan(1,1,nAntiPairs);
                meanPosErrorYSub = NaN;
                meanPerCorrectYSub = NaN;
                meanPosErrorSpatialYSub = nan(1,1,nBins);
                meanPerCorrectSpatialYSub = nan(1,1,nBins);
            end

            % identify true day index
            trueIdx = trueDays(ii,useDays(jj));
            if trueIdx>useDays(jj)
                % set skipped day data to nan
                if trueDays(ii,useDays(jj)-1)==trueIdx-2
                    decodingData.meanPosError(ii,jj) = NaN;
                    decodingData.meanPerCorrect(ii,jj) = NaN;
                    decodingData.meanPosErrorSpatial(ii,jj,:) = nan(1,1,nBins);
                    decodingData.meanPerCorrectSpatial(ii,jj,:) = nan(1,1,nBins);
                    decodingData.meanPerCorrectCues(ii,jj,:) = nan(1,1,nPairs);
                    decodingData.meanPerCorrectAntiCues(ii,jj,:) = nan(1,1,nAntiPairs);
                    decodingData.meanPosErrorYSub(ii,jj) = NaN;
                    decodingData.meanPerCorrectYSub(ii,jj) = NaN;
                    decodingData.meanPosErrorSpatialYSub(ii,jj,:) = nan(1,1,nBins);
                    decodingData.meanPerCorrectSpatialYSub(ii,jj,:) = nan(1,1,nBins);
                end

                % skip days past day limits
                if trueIdx>max(useDays)
                    decodingData.complete(ii,jj) = true;
                    continue
                end
            end

            % store daily activity matrices
            storeIdx = trueIdx-useDays(1)+1;
            decodingData.meanPosError(ii,storeIdx) = meanPosError;
            decodingData.meanPerCorrect(ii,storeIdx) = meanPerCorrect;
            decodingData.meanPosErrorSpatial(ii,storeIdx,:) = meanPosErrorSpatial;
            decodingData.meanPerCorrectSpatial(ii,storeIdx,:) = meanPerCorrectSpatial;
            decodingData.meanPerCorrectCues(ii,storeIdx,:) = meanPerCorrectCues;
            decodingData.meanPerCorrectAntiCues(ii,storeIdx,:) = meanPerCorrectAntiCues;
            decodingData.meanPosErrorYSub(ii,storeIdx) = meanPosErrorYSub;
            decodingData.meanPerCorrectYSub(ii,storeIdx) = meanPerCorrectYSub;
            decodingData.meanPosErrorSpatialYSub(ii,storeIdx,:) = meanPosErrorSpatialYSub;
            decodingData.meanPerCorrectSpatialYSub(ii,storeIdx,:) = meanPerCorrectSpatialYSub;

            % save current data
            decodingData.complete(ii,jj) = true;
            save([p1 '\data\decoding\decodingData_' curCellType '.mat'],'decodingData')

        end
    end

    cd(p1)

    save(['data\decoding\decodingData_' curCellType '.mat'],'decodingData')


    %% Concatenate activity data

    decodingDataGroup = struct();
    fldNames = fieldnames(decodingData);

    % group decoding information by genotype group
    for ff = 1:length(fldNames)-1
        for ss = 1:length(sexes)
            decodingDataGroup.(fldNames{ff}).(sexIDs{ss}) = cell(nGroups,nDays);
            for gg = 1:nGroups
                curMice = intersect(groups{gg},sexes{ss});
                for jj = 1:nDays
                    decodingDataGroup.(fldNames{ff}).(sexIDs{ss}){gg,jj} =...
                        squeeze(decodingData.(fldNames{ff})(curMice,jj,:));
                end
            end
        end
    end

    save(['data\decoding\decodingData_' curCellType '.mat'],'decodingData','decodingDataGroup')


    %% Plot overall decoding
load('D:\AD_Project\imagingData\data\decoding\decodingData4_all.mat')

    close all

    % define pair details
    colors = {[0 0 1],[1 0 0]};
    legs = {'WT','AD'};
    useDayCats = {1,7:10,1:10};
    nCats = length(useDayCats);
    Xlab = {'Pre-Learning','Post-Learning'};

    fldPlot = {'meanPerCorrect','meanPerCorrectYSub','meanPosError','meanPosErrorYSub'};
    fldTitles = {'% Correct: 5cm - 595cm','% Correct: 150cm - 450cm',...
        'Error: 5cm - 595cm', 'Error: 150cm - 450cm'};
    xLimits = [0.5 nDays+0.5];
    % yLimits = [0.2 0.65];
    nFields = length(fldPlot);
    yLabs = {'% correct bins','% correct bins','Decoding error (cm)','Decoding error (cm)'};

    useSexes = 1:3;

    for ss = 1:length(useSexes)
        curSex = useSexes(ss);
        curSexID = sexIDs{curSex};

        h1 = figure;
        h2 = figure;
        for ff = 1:nFields
            %% Plot by day

            % current data
            curDataCat = cell(nGroups,1);
            curData = decodingDataGroup.(fldPlot{ff}).(curSexID);
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
            pShow = generatePShow(nCats,[1 2]);
            barGroup(Xlab,catData,colors,pShow);

            % set plot labels
            title(fldTitles{ff})
            ylabel(yLabs{ff})
            set(gca,'FontSize',12)

        end

        figure(h1)
        sgtitle(['Mean Decoding Correct Rate: ' curCellType ' cells, ' curSexID ' mice (nSub=30)'])
        savefig([svFolder '\decoding_All_' curCellType '_' curSexID])

        figure(h2)
        sgtitle(['Mean Decoding Correct Rate: ' curCellType ' cells, ' curSexID ' mice (nSub=30)'])
        savefig([svFolder '\decoding_All-bar_' curCellType '_' curSexID])


        %% Plot spatial decoding


        fldSpatial = {'meanPosErrorSpatial','meanPerCorrectSpatial'};
        fldTtls = {'Spatial Decoding Correct Rate','Spatial Decoding Error'};
        yLabsSp = {'% correct bins','Decoding error (cm)'};

        nFldSpatial = length(fldSpatial);
        for ff = 1:nFldSpatial
            figure; hold on

            % define cue/reward info
            colorsCue = [0 0 0;0.5 0.5 0.5];
            colorsRew = [0.75 1 1];
            cueFolder = 'D:\AnalysisCode\PostAnalysis\Cues\6mEnv2\';
            rewLoc = [510 560];

            % load cue templates
            load([cueFolder 'tempL.mat'])
            load([cueFolder 'tempR.mat'])
            cueTemp = [tempL,tempR];
            cueX = 2.5:5:597.5;
            colors = {[0 0 1],[1 0 0]};

            for dd = 1:length(useDayCats)
                subplot(1,length(useDayCats),dd); hold on
                curDays = useDayCats{dd};

                curData = cell(nGroups,1);
                for gg = 1:nGroups
                    tempData = cat(3,decodingDataGroup.(fldSpatial{ff}).(curSexID){gg,curDays});
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

                % plot cues
                cueLvl = 1;
                for cc = 1:size(cueTemp,2)
                    plotCues(cueX,cueTemp(:,cc),cueLvl,colorsCue(cc,:));
                end

                % plot data
                plotErrorSig(cueX(nanCol),curData{1}(:,nanCol),curData{2}(:,nanCol),legs,pAnova,pMC,colors)

                % set labels
                xlabel('Track Location (cm)')
                ylabel(yLabsSp{ff})
                title(['Day ' num2str(curDays)])
                set(gca,'fontsize',12)
            end

            sgtitle([fldTtls{ff} ': ' curCellType ' cells (' curSexID ' mice)' ])
            savefig([svFolder '\decoding_' fldSpatial{ff} '_' curCellType '_' curSexID])
        end

        %% Plot cue decoding

        % % cue pair type
        % t1 = 1:size(cuePairInfo,1);
        % t2 = find(cuePairInfo(:,7)==0);
        % t3 = find(cuePairInfo(:,7)==1);
        % t4 = find(cuePairInfo(:,8)==0);
        % t5 = find(cuePairInfo(:,8)==1);
        % t6 = find(cuePairInfo(:,9)==0);
        % t7 = find(cuePairInfo(:,9)==1);
        % typePairs = {t1,t2,t3,t4,t5,t6,t7};
        % 
        % 
        % % only look at all pairs
        % for tt = 1%:length(typePairs)
        %     for cp = 1:2
        %         if cp==1
        %             curData = decodingDataGroup.meanPerCorrectCues.(curSexID);
        %             ttlSuff = 'Cue Decoding Correct Rate: ';
        %             svSuff = 'decoding_CueDist';
        %             tGroup = typePairs{tt};
        %         else
        %             curData = decodingDataGroup.meanPerCorrectAntiCues.(curSexID);
        %             ttlSuff = 'Anti Cue Decoding Correct Rate: ';
        %             svSuff = 'decoding_AntiCueDist';
        %             tGroup = 1:nAntiPairs;
        %         end
        % 
        %         %% Plot by day
        % 
        %         % define pair details
        %         colors = {[0 0 1],[1 0 0]};
        % 
        %         % current data
        %         curDataCat = cell(nGroups,1);
        %         curDataMean = cellfun(@(x) mean(x(:,tGroup),2,'omitnan'),curData,'UniformOutput',false);
        %         for gg = 1:nGroups
        %             curDataCat{gg} = cat(2,curDataMean{gg,:});
        %         end
        % 
        %         % calculate anova p values with multiple comparisons
        %         [pACur,pMC] = anovaRM2W_full_BH(curDataCat{1},curDataCat{2},1);
        %         pAnova = pACur([1 3])'.*ones(nDays,2);
        % 
        %         % plot data
        %         figure; hold on
        %         plotErrorSig(1:nDays,curDataCat{1},curDataCat{2},legs,pAnova,pMC,colors)
        % 
        %         % set labels
        %         xlabel('NE day')
        %         ylabel('correct cue %')
        %         title([ttlSuff curCellType  ' cells (' curSexID ' mice)'])
        %         xlim([0.5 nDays+0.5])
        % 
        %         savefig([svFolder '\' svSuff '_' curCellType '_' curSexID])
        % 
        % 
        %         %% Plot before versus after learning
        % 
        %         % calculate statistics for category days
        %         catData = cell(nCats,nGroups);
        %         for cc = 1:nCats
        %             for gg = 1:nGroups
        %                 catData{cc,gg} = mean(curDataCat{gg}(:,useDayCats{cc}),2,'omitnan');
        %             end
        %         end
        % 
        %         % plot bar graph
        %         figure; hold on
        %         pShow = [1 2;1 3;2 4;3 4];
        %         barGroup(Xlab,catData,colors,pShow);
        % 
        %         % set plot labels
        %         title([ttlSuff curCellType ' cells (' curSexID ' mice)'])
        %         ylabel('correct cue %')
        %         set(gca,'FontSize',12)
        % 
        %         savefig([svFolder '\' svSuff '-bar_' curCellType '_' curSexID])
        %     end
        % end
        % 
        % 
        % %% Plot cue decoding percentile within anti-cue
        % 
        % cueData = decodingDataGroup.meanPerCorrectCues.(curSexID);
        % cueDataMean = cellfun(@(x) mean(x,2,'omitnan'),cueData,'UniformOutput',false);
        % antiData = decodingDataGroup.meanPerCorrectAntiCues.(curSexID);
        % 
        % 
        % tileData = cell(nGroups,nDays);
        % for gg = 1:nGroups
        %     for jj = 1:nDays
        %         curCueData = cueDataMean{gg,jj};
        %         curAntiData = antiData{gg,jj};
        %         tileData{gg,jj} = zeros(size(curCueData,1),1);
        % 
        %         for ii = 1:size(curCueData,1)
        %             sortAnti = sort(curAntiData(ii,:));
        %             curLess = sum(sortAnti<curCueData(ii));
        %             curEqual = sum(sortAnti==curCueData(ii));
        % 
        %             curTile = (curLess+curEqual/2)/nAntiPairs*100;
        %             if curTile>100
        %                 return
        %             end
        %             tileData{gg,jj}(ii) = curTile;
        %         end
        %     end
        % end
        % 
        % % current data
        % curDataCat = cell(nGroups,1);
        % for gg = 1:nGroups
        %     curDataCat{gg} = cat(2,tileData{gg,:});
        % end
        % 
        % % calculate anova p values with multiple comparisons
        % [pACur,pMC] = anovaRM2W_full_BH(curDataCat{1},curDataCat{2},1);
        % pAnova = pACur([1 3])'.*ones(nDays,2);
        % 
        % % plot data
        % figure; hold on
        % plotErrorSig(1:nDays,curDataCat{1},curDataCat{2},legs,pAnova,pMC,colors)
        % 
        % % set labels
        % xlabel('NE day')
        % ylabel('Percentile')
        % title(['Cue Decoding Percentile: ' curCellType ' cells (' curSexID ' mice)'])
        % xlim([0.5 nDays+0.5])
        % 
        % savefig([svFolder '\decoding_CueTile_' curCellType '_' curSexID])
        % 
        % close all
    end
end