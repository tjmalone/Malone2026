%% cueIdentityDiscrimination
% Calculate cue identity discrimination for all FOV and plot data separated
% by funtional cell type, morphology, and sex. Initial implementation uses
% all cells, not just common cells, of all types. Uses learning data with
% cue lag shifted activity.
%

clear; clc; close all;

p1 = 'D:\AD_Project\imagingData';
cd(p1)

% load alignment data
load('foldersLearning.mat')
nFOV = length(alignsLearning);
useDays = 2:size(alignsLearning{1},2);
nDays = length(useDays);

% load cue pair info
load('analysis_Decoding\cuePairInfo_NE.mat')
load('analysis_Decoding\antiCuePairInfo_NE.mat')
nPairs = size(cuePairInfo,1);
nAntiPairs = size(antiCuePairs,1);

pairTypes = {'cue','anticue'};
nPT = length(pairTypes);

% number of cells per correlation group
nSubCell = 30;

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

% load morphology information
load('D:\AD_Project\imagingData\data\morph\morphIdxs.mat','morphIdxs')
useMorphs = morphIdxs.commonCells;
fldNamesM = fieldnames(useMorphs);
useFieldsM = 1:2;
nMorph = length(useFieldsM);

% set cell type categories
cellTypes = {'all','common','grid','cue','ste','pyr'};
nCellType = length(cellTypes);

svFolder = [p1 '\Figures\Figures' char(datetime('today','format','yyyyMMdd')) '\cueDiscrimination'];
if ~isfolder(svFolder)
    mkdir(svFolder)
end

% whether to rerun activity extraction
newRun = 0;


%% Extract activity data

% loop through cell types
for cType = 1:nCellType
    curCellType = cellTypes{cType};

    if isfile(['data\cueIdent\cueIdentityData_' curCellType '.mat']) & newRun==0
        load(['data\cueIdent\cueIdentityData_' curCellType '.mat'])
    else
        cueIdentityData = struct();
    end

    for pt = 1:nPT
        if newRun==1 || ~isfield(cueIdentityData,(pairTypes{pt}))
            cueIdentityData.(pairTypes{pt}).ampDiffMean = cell(nFOV,nDays);
            % cueIdentityData.(pairTypes{pt}).ampDiffPeak = cell(nFOV,nDays);
            cueIdentityData.(pairTypes{pt}).corrDiffMean = cell(nFOV,nDays);
            % cueIdentityData.(pairTypes{pt}).corrDiffPeak = cell(nFOV,nDays);
            cueIdentityData.(pairTypes{pt}).complete = false(nFOV,nDays);
        end

        % current pair number
        if pt==1
            nPairsCur = size(cuePairs,1);
        elseif pt==2
            nPairsCur = size(antiCuePairs,1);
        end

        % loop through FOV
        for ii = 1:nFOV
            if all(cueIdentityData.(pairTypes{pt}).complete(ii,:)==true)
                continue
            end

            % loop through days
            for jj = 1:nDays
                %% Calculate cue identity discrimination

                if cueIdentityData.(pairTypes{pt}).complete(ii,jj)==1; continue; end

                cd(foldersLearning{ii}{useDays(jj)})
                disp([num2str(ii) '-' num2str(jj)])

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
                elseif strcmp(curCellType,'ste')
                    useCells = useMorphs.ste{ii,useDays(jj)};
                elseif strcmp(curCellType,'pyr')
                    useCells = useMorphs.pyr{ii,useDays(jj)};
                end
                nUseCells = sum(useCells~=0);

                % check cell number
                flag = false;
                if nUseCells<nSubCell || cType>=5
                    flag = true;
                end

                if pt==1
                    [curAmpDiffMean,curAmpDifPeak] = ampDiff_TM({useCells},[],cuePairInfo);
                    if ~flag
                        [curCorrDiffMean,curCorrDifPeak] =...
                            corrAtCues_TM({useCells},[],cuePairInfo,nSubCell);
                    end
                else
                    [curAmpDiffMean,curAmpDifPeak] = ampDiff_TM({useCells},[],antiCuePairInfo,antiCuePairs);
                    if ~flag
                        [curCorrDiffMean,curCorrDifPeak] =...
                            corrAtCues_TM({useCells},[],antiCuePairInfo,nSubCell,[],antiCuePairs);
                    end
                end

                if flag
                    curCorrDiffMean = {nan(nUseCells,nPairsCur)};
                    curCorrDifPeak = {nan(nUseCells,nPairsCur)};
                end

                % identify true day index
                trueIdx = trueDays(ii,useDays(jj));
                if trueIdx>useDays(jj)
                    % set skipped day data to nan
                    if trueDays(ii,useDays(jj)-1)==trueIdx-2
                        cueIdentityData.(pairTypes{pt}).ampDiffMean{ii,jj} = nan(nUseCells,nPairsCur);
                        % cueIdentityData.(pairTypes{pt}).ampDiffPeak{ii,jj} = nan(nUseCells,nPairsCur);
                        cueIdentityData.(pairTypes{pt}).corrDiffMean{ii,jj} = nan(nUseCells,nPairsCur);
                        % cueIdentityData.(pairTypes{pt}).corrDiffPeak{ii,jj} = nan(nUseCells,nPairsCur);
                    end

                    % skip days past day limits
                    if trueIdx>max(useDays)
                        continue
                    end
                end

                % store daily activity matrices
                cueIdentityData.(pairTypes{pt}).ampDiffMean{ii,trueIdx-1} = curAmpDiffMean{1};
                % cueIdentityData.(pairTypes{pt}).ampDiffPeak{ii,trueIdx-1} = curAmpDifPeak{1};
                cueIdentityData.(pairTypes{pt}).corrDiffMean{ii,trueIdx-1} = curCorrDiffMean{1};
                % cueIdentityData.(pairTypes{pt}).corrDiffPeak{ii,trueIdx-1} = curCorrDifPeak{1};

            end
            % save current data
            cueIdentityData.(pairTypes{pt}).complete(ii,:) = true;
            save([p1 '\data\cueIdent\cueIdentityData_' curCellType '.mat'],'cueIdentityData','-v7.3')
        end

        cd(p1)
    end


    %% Concatenate activity data

    cueIdentityDataCat = struct();

    for pt = 1:nPT
        fldNames = fieldnames(cueIdentityData.(pairTypes{pt}));
        for ff = 1:length(fldNames)-1
            for ss = 1:length(sexes)
                cueIdentityDataCat.(fldNames{ff}).(sexIDs{ss}).(pairTypes{pt}) = cell(nGroups,nDays);

                for gg = 1:nGroups
                    curMice = intersect(groups{gg},sexes{ss});

                    for jj = 1:nDays
                        cueIdentityDataCat.(fldNames{ff}).(sexIDs{ss}).(pairTypes{pt}){gg,jj} =...
                            cat(1,cueIdentityData.(pairTypes{pt}).(fldNames{ff}){curMice,jj});
                    end
                end
            end
        end
    end

    save(['data\cueIdent\cueIdentityDataCat_' curCellType '.mat'],'cueIdentityDataCat','-v7.3')
    clear cueIdentityData


    %% Plot cue identity discrimination

    close all

    % define plot fields
    if cType<=4
        fldNames = {'ampDiffMean','corrDiffMean'};
    elseif cType
        fldNames = {'ampDiffMean'};
    end
    nFields = length(fldNames);

    % define legend
    legs = groupIDs;

    % define plot colors
    colors = {[0 0 1],[1 0 0]};

    % define bar plot criteria
    dayCats = {1,[7 8 9 10]};
    ffBar = [1,2];
    ptBar = 1;
    sttBar = 1;
    Xlab = {'Pre-Learning','Post-Learning'};

    t1 = 1:size(cuePairInfo,1);
    t2 = find(cuePairInfo(:,7)==0);
    t3 = find(cuePairInfo(:,7)==1);
    t4 = find(cuePairInfo(:,8)==0);
    t5 = find(cuePairInfo(:,8)==1);
    t6 = find(cuePairInfo(:,9)==0);
    t7 = find(cuePairInfo(:,9)==1);

    typePairs = {t1,t2,t3,t4,t5,t6,t7};
    tGroups = {1,2:3,4:5,6:7};
    typeNames = {'All','Pattern','Side','Shape'};
    typeSubNames = {'Different','Identical'};

    yLabs = {'Activity Variation (cell): high is discrimination',...
        'Rank Correlation (population): low is discrimination'};

    nTypes = length(tGroups);

    useSexes = 1:3;

    for ss = 1:length(useSexes)
        curSex = useSexes(ss);
        curSexID = sexIDs{curSex};

        % loop through cue match types
        for tt = 1:nTypes
            %% Plot figure

            % skip pair sub type anlaysis for sex subsets
            if ss>1 && tt>1; continue; end

            h1 = figure; hold on
            set(gca,'FontSize',12)
            nSubplotX = nFields*nPT;
            nSubPlotY = length(tGroups{tt});

            % loop through fields
            for ff = 1:nFields
                yLimits = nan(1,2);

                % loop through cue types (cue/anti-cue)
                for pt = 1:nPT

                    % select current data
                    dataCur = cueIdentityDataCat.(fldNames{ff}).(curSexID).(pairTypes{pt});

                    % loop through cue match sub types (different/idenitical)
                    for stt = 1:length(tGroups{tt})
                        %% Plot across days

                        if pt==1
                            useCols = typePairs{tGroups{tt}(stt)};
                        else
                            useCols = 1:size(antiCuePairs,1);
                        end
                        dataSub = cellfun(@(x) mean(x(:,useCols),2,'omitnan'),dataCur,'UniformOutput',false);

                        % calcualte p values
                        [pACur,pMC] = anovaO2W_BH(dataSub,1);
                        pAnova = pACur([1 3])'.*ones(nDays,2);

                        % plot line graphs
                        figure(h1)
                        subplot(nSubPlotY,nSubplotX,(stt-1)*nSubplotX+(ff-1)*nPT+pt); hold on
                        plotErrorSigCell(1:nDays,dataSub,legs,pAnova,pMC',colors)

                        % set plot labels
                        title([fldNames{ff} '-' pairTypes{pt} ': ' typeSubNames{stt}])
                        xlabel('NE Day')
                        ylabel(yLabs{ff})
                        curYLimits = ylim;
                        yLimits(1) = min(yLimits(1),curYLimits(1));
                        yLimits(2) = max(yLimits(2),curYLimits(2));
                        set(gca,'FontSize',12)


                        %% Plot pre- versus post- learning

                        if ismember(ff,ffBar) && ismember(pt,ptBar) && ismember(tGroups{tt}(stt),sttBar)
                            % calculate statistics for all days
                            subMeans = cellfun(@(x) mean(x,'omitnan'),dataSub);
                            subSTD = cellfun(@(x) std(x,'omitnan'),dataSub);
                            subSEM = cellfun(@(x) nansem(x,1),dataSub);
                            subNs = cellfun(@(x) sum(~isnan(x)),dataSub);

                            % calculate statistics for category days
                            catM = zeros(2,nGroups);
                            catSTD = zeros(2,nGroups);
                            catSEM = zeros(2,nGroups);
                            catN = zeros(2,nGroups);
                            for cc = 1:2
                                catM(cc,:) = mean(subMeans(:,dayCats{cc}),2,'omitnan');
                                catSTD(cc,:) = mean(subSTD(:,dayCats{cc}),2,'omitnan');
                                catSEM(cc,:) = mean(subSEM(:,dayCats{cc}),2,'omitnan');
                                catN(cc,:) = mean(subNs(:,dayCats{cc}),2,'omitnan');
                            end

                            % plot bar graph
                            figure; hold on
                            pShow = [1 2;1 3;2 4;3 4];
                            barMan(Xlab,catM,catSEM,catN,catSTD,colors,pShow);

                            % set figure labels
                            title([fldNames{ff} '-' pairTypes{pt} '(' curCellType ' cells, ' curSexID ' mice): Pre- versus Post- Learning'])
                            ylabel(yLabs{ff})
                            set(gca,'fontsize',12)
                            savefig([svFolder '\cueDiscrim_' typeNames{tt} '-' fldNames{ff} '_' curCellType '-bar_' curSexID '.fig'])
                        end
                    end
                end

                figure(h1)
                for pt = 1:nPT
                    for stt = 1:length(tGroups{tt})
                        subplot(nSubPlotY,nSubplotX,(stt-1)*nSubplotX+(ff-1)*nPT+pt)
                        ylim(yLimits)
                    end
                end
            end

            sgtitle([typeNames{tt} ': ' curCellType ' cells, ' curSexID ' mice'])

            savefig([svFolder '\cueDiscrim_' typeNames{tt} '_' curCellType '_' curSexID '.fig'])


            %% Plot cue percentile within anti-cue

            if tt==1
                figure
                for ff = 1:nFields
                    cueData = cueIdentityDataCat.(fldNames{ff}).(curSexID).cue;
                    cueDataMean = cellfun(@(x) mean(x,2,'omitnan'),cueData,'UniformOutput',false);
                    antiData = cueIdentityDataCat.(fldNames{ff}).(curSexID).anticue;

                    tileData = cell(nGroups,nDays);
                    for gg = 1:nGroups
                        for jj = 1:nDays
                            curNans = isnan(cueDataMean{gg,jj});
                            curCueData = cueDataMean{gg,jj}(curNans==0);
                            curAntiData = antiData{gg,jj}(curNans==0,:);
                            tileData{gg,jj} = zeros(size(curCueData,1),1);

                            for ii = 1:size(curCueData,1)
                                sortAnti = sort(curAntiData(ii,:));
                                curLess = sum(sortAnti<curCueData(ii));
                                curEqual = sum(sortAnti==curCueData(ii));

                                curTile = (curLess+curEqual/2)/sum(~isnan(sortAnti))*100;
                                if curTile>100
                                    return
                                end
                                tileData{gg,jj}(ii) = curTile;
                            end
                        end
                    end

                    % calcualte p values
                    [pACur,pMC] = anovaO2W_BH(tileData,1);
                    pAnova = pACur([1 3])'.*ones(nDays,2);

                    % plot data
                    subplot(1,nFields,ff); hold on
                    plotErrorSigCell(1:nDays,tileData,legs,pAnova,pMC',colors)

                    % set plot labels
                    title(fldNames{ff})
                    xlabel('NE Day')
                    ylabel('Percentile')
                    set(gca,'FontSize',12)
                end
                sgtitle(['Cue Identity Percentile: ' curCellType ' cells, ' curSexID ' mice'])
                savefig([svFolder '\cueDiscrim_Tile_' curCellType '_' curSexID])

            end
        end

        close all
    end
end


%% Plot all cue only

close all

fldNames = {'ampDiffMean','corrDiffMean'};
nFields = length(fldNames);

% define legend
legs = groupIDs;

% define plot colors
colors = {[0 0 1],[1 0 0]};

% define bar plot criteria
dayCats = {1,[7 8 9 10]};
Xlab = {'Pre-Learning','Post-Learning'};
yLabs = {'Activity Variation (cell): high is discrimination',...
    'Rank Correlation (population): low is discrimination'};

clear pt tt stt ff

nSubplot = ceil(nCellType^0.5);
spSz = [nSubplot,nSubplot];
spSz = [2,3];

for ss = 1:length(useSexes)
    curSex = useSexes(ss);
    curSexID = sexIDs{curSex};

    % loop through fields
    for ff = 1:nFields
        % initialize figures
        h1 = figure; hold on
        set(gca,'FontSize',12)
        h2 = figure; hold on
        set(gca,'FontSize',12)

        % initialize plot info
        
        yLimits = nan(2,2);

        % loop through cell types
        for cType = 1:nCellType
            % skip invalid combinations
            if ff==2 && cType>=5; continue; end

            % load cell type data
            curCellType = cellTypes{cType};
            load(['data\cueIdent\cueIdentityDataCat_' curCellType '.mat'])

            % select current data
            dataCur = cueIdentityDataCat.(fldNames{ff}).(curSexID).(pairTypes{1});

            %% Plot across days

            dataSub = cellfun(@(x) mean(x,2,'omitnan'),dataCur,'UniformOutput',false);

            % calcualte p values
            [pACur,pMC] = anovaO2W_BH(dataSub,1);
            pAnova = pACur([1 3])'.*ones(nDays,2);

            % plot line graphs
            figure(h1)
            subplot(spSz(1),spSz(2),cType); hold on
            plotErrorSigCell(1:nDays,dataSub,legs,pAnova,pMC',colors)

            % set plot labels
            title([curCellType ' cells'])
            xlabel('NE Day')
            ylabel(yLabs{ff})
            curYLimits = ylim;
            yLimits(1,1) = min(yLimits(1,1),curYLimits(1));
            yLimits(1,2) = max(yLimits(1,2),curYLimits(2));
            set(gca,'FontSize',12)


            %% Plot pre- versus post- learning

            % calculate statistics for all days
            subMeans = cellfun(@(x) mean(x,'omitnan'),dataSub);
            subSTD = cellfun(@(x) std(x,'omitnan'),dataSub);
            subSEM = cellfun(@(x) nansem(x,1),dataSub);
            subNs = cellfun(@(x) sum(~isnan(x)),dataSub);

            % calculate statistics for category days
            catM = zeros(2,nGroups);
            catSTD = zeros(2,nGroups);
            catSEM = zeros(2,nGroups);
            catN = zeros(2,nGroups);
            for cc = 1:2
                catM(cc,:) = mean(subMeans(:,dayCats{cc}),2,'omitnan');
                catSTD(cc,:) = mean(subSTD(:,dayCats{cc}),2,'omitnan');
                catSEM(cc,:) = mean(subSEM(:,dayCats{cc}),2,'omitnan');
                catN(cc,:) = mean(subNs(:,dayCats{cc}),2,'omitnan');
            end

            % plot bar graph
            figure(h2)
            subplot(spSz(1),spSz(2),cType); hold on
            pShow = [1 2;1 3;2 4;3 4];
            barMan(Xlab,catM,catSEM,catN,catSTD,colors,pShow);

            % set figure labels
            title([curCellType ' cells'])
            ylabel(yLabs{ff})
            curYLimits = ylim;
            yLimits(2,1) = min(yLimits(2,1),curYLimits(1));
            yLimits(2,2) = max(yLimits(2,2),curYLimits(2));
            set(gca,'fontsize',12)
        end

        % set axis limits
        for cType = 1:nCellType
            % set figure 1 axes
            figure(h1)
            subplot(spSz(1),spSz(2),cType)
            ylim(yLimits(1,:))

            % set figure 2 axes
            figure(h2)
            subplot(spSz(1),spSz(2),cType)
            ylim([yLimits(1,1),yLimits(2,2)])
            if ff==1
                ylim([yLimits(1,1) 0.8])
            end
        end

        % set title and save
        figure(h1)
        sgtitle([fldNames{ff} ' (all cell types, ' curSexID ' mice)'])
        savefig([svFolder '\cueDiscrim_allTypes-' fldNames{ff} '_' curSexID '.fig'])

        figure(h2)
        sgtitle([fldNames{ff} ': Pre- versus Post- Learning (all cell types, ' curSexID ' mice)'])
        savefig([svFolder '\cueDiscrim_allTypes-' fldNames{ff} '-bar_' curSexID '.fig'])
    end

    close all
end
