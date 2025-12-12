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
useDays = 1:size(alignsLearning{1},2);
nDays = length(useDays);

% set decoding settings
trackLength = 600;
binWidth = 5;
nBins = trackLength/binWidth;
nSubCell = 30;  % number of cels per iteration
nItr = 100;     % number of iterations

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

svFolder = [p1 '\Figures\Figures' char(datetime('today','format','yyyyMMdd'))  '\decodingMaster'];
if ~isfolder(svFolder)
    mkdir(svFolder)
end

% set cell type categories
cellTypes = {'all','common','grid','nongrid','cue','other'};
% cellTypes = {'all'};
nCellType = length(cellTypes);

% set random seed
rng(42)

% whether to rerun activity extraction
newRun = 0;


%% Extract activity data

% loop through cell types
for cType = 1:nCellType
    curCellType = cellTypes{cType};
    curFile = [p1 '\data\decoding\decodingDataMaster_' curCellType '.mat'];

    if isfile(curFile) && newRun==0
        load(curFile)
        continue
    else
        decodingData = struct();
        decodingData.errorBin = cell(nFOV,nDays);
        decodingData.testYBins = cell(nFOV,nDays);
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
                decdPos = decodingByLocationMaster_TM(dfofUse,abf,trackLength,...
                    binWidth,nItr,nSubCell);

                % store decoding data
                errorBin = decdPos.errorBin;
                testYBins = decdPos.testYBins;

            else
                errorBin = nan(nItr,1);
                testYBins = nan(nItr,1);
            end

            % identify true day index
            trueIdx = trueDays(ii,useDays(jj));
            if trueIdx>useDays(jj)
                % set skipped day data to nan
                if trueDays(ii,useDays(jj)-1)==trueIdx-2
                    decodingData.errorBin{ii,jj} = nan(nItr,1);
                    decodingData.testYBins{ii,jj} = nan(nItr,1);
                end

                % skip days past day limits
                if trueIdx>max(useDays)
                    decodingData.complete(ii,jj) = true;
                    continue
                end
            end

            % store daily activity matrices
            storeIdx = trueIdx-useDays(1)+1;
            decodingData.errorBin{ii,storeIdx} = errorBin;
            decodingData.testYBins{ii,storeIdx} = testYBins;

            % save current data
            % decodingData.complete(ii,jj) = true;
            % save(curFile,'decodingData')
        end
    end

    cd(p1)

    save(curFile,'decodingData')


    %% Process decoding

    % set decoding parameters
    binErrorThresh = 0:4;                   % criteria for correct decoding
    excludeBins = [0,1,2,4,6,8,10,15,25,50];  % number of bins to exclude at track start and end
    trackLoop = 1;                          % whether to loop the track for errors

    nThresh = length(binErrorThresh);
    nExclude = length(excludeBins);

    % load y positions
    testYBins = decodingData.testYBins;

    % perform track looping
    if trackLoop==1
        errorBinAdj = cellfun(@(x) min(x,nBins-x),decodingData.errorBin,'UniformOutput',false);
    else
        errorBinAdj = decodingData.errorBin;
    end

    % initialize data structure
    decodingDataItr = struct();
    itrN = 0;

    for th = 1:nThresh
        % set current correct threshold
        curThresh = binErrorThresh(th);

        % calculate spatial decoding
        meanPosErrorSpatial = zeros(nFOV,nDays,nBins);
        meanPerCorrectSpatial = zeros(nFOV,nDays,nBins);
        for binNumber = 1:nBins
            curErrors = cellfun(@(x,y) x(:,y==binNumber),errorBinAdj,testYBins,'UniformOutput',false);
            meanPosErrorSpatial(:,:,binNumber) = cellfun(@(x) mean(x,'all','omitnan')*binWidth,curErrors);
            meanPerCorrectSpatial(:,:,binNumber) = cellfun(@(x) mean(mean(x<=curThresh,2)),curErrors);
        end

        % calculate general decoding
        for ex = 1:nExclude
            % set current use bins
            curExclude = excludeBins(ex);
            binSubset = 1+curExclude:nBins-curExclude;

            % apply exclusion filter
            keepIdx = cellfun(@(x) ismember(x,binSubset),testYBins,'UniformOutput',false);

            % reload errors
            errorBinCur = cellfun(@(x,y) x(:,y),errorBinAdj,keepIdx,'UniformOutput',false);
            yBinCur = cellfun(@(x,y) x(y),testYBins,keepIdx,'UniformOutput',false);

            % calculate overall decoding statistics
            meanPosError = cellfun(@(x) mean(x,'all','omitnan')*binWidth,errorBinCur);
            meanPerCorrect = cellfun(@(x) mean(mean(x<=curThresh,2)),errorBinCur);

            % reset iteration index
            itrN = itrN+1;

            % store results
            decodingDataItr(itrN).thresh = curThresh;
            decodingDataItr(itrN).exclude = curExclude;
            decodingDataItr(itrN).meanPosError = meanPosError;
            decodingDataItr(itrN).meanPerCorrect = meanPerCorrect;
            decodingDataItr(itrN).meanPosErrorSpatial = meanPosErrorSpatial;
            decodingDataItr(itrN).meanPerCorrectSpatial = meanPerCorrectSpatial;
        end
    end

    save(curFile,'decodingData','decodingDataItr')


    %% Separate data by group

    decodingDataGroup = struct();
    fldNames = fieldnames(decodingDataItr);
    splitIdx = 3;

    for ii = 1:itrN
        % group decoding information by genotype group

        for ff = 1:length(fldNames)
            if ff<splitIdx
                decodingDataGroup(ii).(fldNames{ff}) = decodingDataItr(ii).(fldNames{ff});
            else
                % split mice by sex and genotype
                for ss = 1:length(sexes)
                    decodingDataGroup(ii).(fldNames{ff}).(sexIDs{ss}) = cell(nGroups,nDays);
                    for gg = 1:nGroups
                        curMice = intersect(groups{gg},sexes{ss});
                        for jj = 1:nDays
                            decodingDataGroup(ii).(fldNames{ff}).(sexIDs{ss}){gg,jj} =...
                                squeeze(decodingDataItr(ii).(fldNames{ff})(curMice,jj,:));
                        end
                    end
                end
            end
        end
    end

    save(curFile,'decodingData','decodingDataItr','decodingDataGroup')
end

return
%% Plot overall decoding

close all
load('\data\decoding\decodingDataMaster_all.mat');

% define pair details
colors = {[0 0 1],[1 0 0]};
legs = {'WT','AD'};

% define spatial/cue/reward info
spDays = 2:11;
colorsCue = [0 0 0;0.5 0.5 0.5];
colorsRew = [0.75 1 1];
cueFolder = 'D:\AnalysisCode\PostAnalysis\Cues\6mEnv2\';
rewLoc = [510 560];

% load cue templates
load([cueFolder 'tempL.mat'])
load([cueFolder 'tempR.mat'])
cueTemp = [tempL,tempR];
cueX = 2.5:5:597.5;

fldPlot = {'meanPerCorrect','meanPosError';'meanPerCorrectSpatial','meanPosErrorSpatial'};
fldTitles = {'% Correct Bins: ','Looped Decoding Error: '};
nFields = numel(fldPlot);
xLimits = [0.5 nDays+0.5];
% yLimits = [0.2 0.65];
yLabs = {'% correct bins','Decoding error (cm)'};

useSexes = 2:3;
nSexes = length(useSexes);

plotDays = 2:nDays;


for ii = 1:itrN
    curDecodeData = decodingDataGroup(ii);

    h1 = figure;
    tiledlayout(nSexes,nFields)

    for ss = 1:nSexes
        curSex = useSexes(ss);
        curSexID = sexIDs{curSex};

        for ff = 1:nFields/2
            %% Plot by day

            % current data
            curDataCat = cell(nGroups,1);
            curData = curDecodeData.(fldPlot{1,ff}).(curSexID);
            for gg = 1:nGroups
                curDataCat{gg} = cat(2,curData{gg,:});
            end

            % calculate anova p values with multiple comparisons
            pAnova = nan(nDays,2);
            pMC = nan(nDays,1);
            try
                [pACur,pMCCur] = anovaRM2W_full_BH(curDataCat{1}(:,spDays),curDataCat{2}(:,spDays),1);
                pAnova(spDays,:) = pACur([1 3])'.*ones(length(spDays),2);
                pMC(spDays) = pMCCur;
            catch
            end

            % plot data
            nexttile((ss-1)*nFields+ff); hold on
            plotErrorSig(1:length(plotDays),curDataCat{1}(:,plotDays),curDataCat{2}(:,plotDays),...
                legs,pAnova(plotDays,:),pMC(plotDays),colors)

            % set plot labels
            xlabel('Session')
            ylabel(yLabs{ff})
            xlim(xLimits)
            title([fldTitles{ff} curSexID])
            set(gca,'fontsize',12)


            %% Plot spatial decoding

            nexttile((ss-0.5)*nFields+ff); hold on

            curData = cell(nGroups,1);
            for gg = 1:nGroups
                tempData = cat(3,curDecodeData.(fldPlot{2,ff}).(curSexID){gg,spDays});
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
            plotErrorSig(cueX(nanCol),curData{1}(:,nanCol),curData{2}(:,nanCol),legs,pAnova,pMC,colors)

            % plot cues
            cueLvl = max(ylim);
            for cc = 1:size(cueTemp,2)
                plotCues(cueX,cueTemp(:,cc),cueLvl,colorsCue(cc,:));
            end

            % set labels
            xlabel('Track Location (cm)')
            ylabel(yLabs{ff})
            title(['Spatial ' fldTitles{ff} curSexID])
            set(gca,'fontsize',12)
        end
    end

    sgtitle(['Decoding: ' curCellType ' cells, thresh=' num2str(curDecodeData.thresh)...
        ', excludeBins=' num2str(curDecodeData.exclude)])
    savefig([svFolder '\decoding_th' num2str(curDecodeData.thresh) '_'...
        'exc' num2str(curDecodeData.exclude)])

    close
end


%% Plot spatial decoding

close all

% define spatial/cue/reward info
spDaysCat = {2:4,9:11};
nSP = length(spDaysCat);

fldPlot = {'meanPerCorrectSpatial','meanPosErrorSpatial'};
fldTitles = {'% Correct Bins: ','Looped Decoding Error: '};
nFields = length(fldPlot);
xLimits = [0.5 nDays+0.5];
% yLimits = [0.2 0.65];
yLabs = {'% correct bins','Decoding error (cm)'};

useSexes = 2:3;
nSexes = length(useSexes);

useItr = [21 30];

for ii = useItr
    curDecodeData = decodingDataGroup(ii);

    h1 = figure;
    tiledlayout(nSexes,nFields*nSP)

    for ss = 1:nSexes
        curSex = useSexes(ss);
        curSexID = sexIDs{curSex};

        for ff = 1:nFields
            for sp = 1:nSP
                %% Plot spatial decoding

                nexttile((ss-1)*nFields*nSP+2*ff+sp-2); hold on

                curData = cell(nGroups,1);
                for gg = 1:nGroups
                    tempData = cat(3,curDecodeData.(fldPlot{ff}).(curSexID){gg,spDaysCat{sp}});
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
                plotErrorSig(cueX(nanCol),curData{1}(:,nanCol),curData{2}(:,nanCol),legs,pAnova,pMC,colors)

                % plot cues
                cueLvl = max(ylim);
                for cc = 1:size(cueTemp,2)
                    plotCues(cueX,cueTemp(:,cc),cueLvl,colorsCue(cc,:));
                end

                % set labels
                xlabel('Track Location (cm)')
                ylabel(yLabs{ff})
                title(['Spatial ' fldTitles{ff} curSexID])
                set(gca,'fontsize',12)
            end
        end
    end

    sgtitle(['Spatial Decoding: ' curCellType ' cells, thresh=' num2str(curDecodeData.thresh)...
        ', excludeBins=' num2str(curDecodeData.exclude)])
    savefig([svFolder '\spatialDecoding_th' num2str(curDecodeData.thresh) '_'...
        'exc' num2str(curDecodeData.exclude)])

    close
end

