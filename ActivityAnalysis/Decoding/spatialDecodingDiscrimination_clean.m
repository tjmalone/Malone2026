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
binErrorThresh = 2; % correct decoding bin threshold
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

    load(['data\decoding\decodingData_' curCellType '.mat'],'decodingData','decodingDataGroup')


    %% Plot overall decoding

    close all

    % define pair details
    colors = {[0 0 1],[1 0 0]};
    legs = {'WT','AD'};
    useDayCats = {9,10,7:10,1:10};
    nCats = length(useDayCats);
    Xlab = {'Day 9','Day 10','Post-Leaning','All NE'};

    fldPlot = {'meanPerCorrect','meanPerCorrectYSub'};
    fldTitles = {'5cm - 595cm','150cm - 450cm'};
    xLimits = [0.5 nDays+0.5];
    yLimits = [0.2 0.65];
    nFields = length(fldPlot);
    yLab = '% correct bins';

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
            ylabel(yLab)
            xlim(xLimits)
            ylim(yLimits)
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
            ylabel(yLab)
            set(gca,'FontSize',12)

        end

        figure(h1)
        sgtitle(['Mean Decoding Correct Rate: ' curCellType ' cells, ' curSexID ' mice (nSub=30)'])
        savefig([svFolder '\decoding_All_' curCellType '_' curSexID])

        figure(h2)
        sgtitle(['Mean Decoding Correct Rate: ' curCellType ' cells, ' curSexID ' mice (nSub=30)'])
        savefig([svFolder '\decoding_All-bar_' curCellType '_' curSexID])


        %% Plot spatial decoding

        ttlSupps = {'Raw','Norm'};
        for nn = 1:2
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
                    tempData = cat(3,decodingDataGroup.meanPerCorrectSpatial.(curSexID){gg,curDays});
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

                if nn==2
                    for gg = 1:nGroups
                        curMean = mean(curData{gg}(:,nanCol),1);
                        valMin = min(curMean);
                        valMax = max(curMean);

                        curData{gg} = (curData{gg}-valMin)/(valMax-valMin);
                    end
                end

                % plot data
                plotErrorSig(cueX(nanCol),curData{1}(:,nanCol),curData{2}(:,nanCol),legs,pAnova,pMC,colors)

                % set labels
                xlabel('Track Location (cm)')
                ylabel(['% correct bins (' ttlSupps{nn} ')'])
                title(['Day ' num2str(curDays)])
                set(gca,'fontsize',12)
            end

            sgtitle([ttlSupps{nn} 'Spatial Decoding Correct Rate: ' curCellType ' cells (' curSexID ' mice)' ])
            savefig([svFolder '\decoding_Spatial_' curCellType '_' curSexID '_' ttlSupps{nn}])

        end
    end
end