%% trajectoryPlot_mouseCorr.m
% Calculate neuronal activity distribution as a function of track location
% per mouse for a given set of days. Caluclate and analyze intermouse map
% correlation split by sex and genotype
%
% Distribution Types:
%   Learning types - learn
%   Run types  - all runs
%   Distribution types - to others RBR
%
% Data Separations:
%   Sex - allSex, female, male
%   Cell type - all, common
%   Morphology type - allMorph, ste, pyr
%   Genotype - WT, AD
%


%% Load data

clear; clc; close all

p1 = 'D:\AD_Project\imagingData';
cd(p1)

% load distribution data
load('data\distActivity.mat','distActivity')

% load cell selections
load('data\cellSelect.mat','cellSelect')

% set save folder name
svFile = [p1 '\Figures\Figures' char(datetime('today','format','yyyyMMdd')) '\trajectory'];
if ~isfolder(svFile)
    mkdir(svFile);
    mkdir([svFile '\data']);
end


%% Set input parameters

% distribution key: learning type (lt), run type (rp), distribution type,
% (dt), success/fail (ds)
distTypes = {'dfof'};

% cell selection key: sex (ss), cell type (ct), morphology type (mt),
% genotype (gg)
sexIDs = {'female','male'};
cellTypes = {'common'};
% cellTypes = {'grid','cue','other'};
morphTypes = {'allMorph'};
% morphTypes = {'allMorph','ste','pyr'};
genotypes = {'WT','AD'};
nGeno = length(genotypes);
nSexes = length(sexIDs);
nMorphs = length(morphTypes);

% define averaging parameters
fovPerMouse = 2;

% define day category and labels
useDays = 1:11;
nDays = length(useDays);
dayNames = 1:nDays;
anovaIdx = 2:11;

% define colors
colors = {[0.5 0.5 0.5],[175 0 0]/255};
colors4 = {[102 105 255]/255,[255 102 102]/255;[0 0 175]/255,[175 0 0]/255};


% close all figures
clAll = 0;


%% Plot distributions

% get current distribution
distBase = distActivity.learn.allRuns.dfof.all;

% initialize stats and settings
outStats = {};
nUnits = 'mouse pairs';

% loop through cell types
for ct = cellTypes

    % loop through morphology types
    for mt = 1:length(morphTypes)
        %% Calculate mouse curves by group

        % skip invalid combination
        if strcmp(ct,'allType') && ~strcmp(morphTypes{mt},'allMorph')
            continue
        end

        % initialize quantification struct
        dataQuant = cell(nSexes,nGeno);

        % loop through sexes
        for ss = 1:nSexes

            % loop through genotypes
            for gg = 1:nGeno
                %% Process current data

                % get current cell selections
                curCellSel = cellSelect.learn.(sexIDs{ss}).(ct{:}).(morphTypes{mt}).(genotypes{gg});

                % correct empty cells
                maxIdx = cellfun(@max,curCellSel,'UniformOutput',false);
                emptyCellSel = cellfun(@isempty,maxIdx);
                emptyData = cellfun(@isempty,distBase);
                fixIdx = find(~emptyCellSel & emptyData);
                refIdx = find(~emptyData,1);
                szXCur = size(distBase{refIdx},2);
                for ii = 1:length(fixIdx)
                    distBase{fixIdx(ii)} = nan(maxIdx{fixIdx(ii)},szXCur);
                end

                % define current data
                curData = cellfun(@(x,y) x(y,:),distBase,...
                    curCellSel,'UniformOutput',false);
                nFOV = size(curData,1);

                % fill in empty cells
                emptyAll = cellfun(@isempty,curData);
                for ff = 1:nFOV
                    % get all empty and first full index
                    curEmpty = find(emptyAll(ff,:));
                    curFull = find(~emptyAll(ff,:),1);

                    % skip full or empty rows
                    if isempty(curEmpty) || isempty(curFull)
                        continue
                    end

                    % fill empty cells
                    for ii = curEmpty
                        curData{ff,ii} = nan(size(curData{ff,curFull}));
                    end
                end

                % fill in remaining empty cells
                curEmpty = find(cellfun(@(x) size(x,2),curData)==0);
                for ii = curEmpty'
                    curData{ii} = zeros(0,szXCur);
                end

                % take mouse mean
                curDataCat = zeros(nFOV/fovPerMouse,nDays,120);
                for ff = 1:nFOV/fovPerMouse

                    curFOV = ff*fovPerMouse-1:ff*fovPerMouse;

                    for dd = 1:nDays
                        curDay = useDays(dd);

                        % get current mouse/day data
                        curDataMouse = cat(1,curData{curFOV,curDay});

                        %take mouse/day mean
                        curDataMean = mean(curDataMouse,1,'omitnan');
                        curDataCat(ff,dd,:) = curDataMean;
                    end
                end

                % remove nan mice
                nanMice = all(isnan(curDataCat),[2 3]);
                curDataCatUse = curDataCat(~nanMice,:,:);

                % take mean across days
                dataQuant{ss,gg} = curDataCatUse;
            end
        end


        %% Calculate inter-mouse correlation

        % reshape data
        dataByGroup = reshape(dataQuant,[],1);

        % initialize combinations
        combs = combvec(1:4, 1:4)';
        combs(combs(:,1)>combs(:,2),:) = [];
        nCombs = size(combs,1);

        % define set info: Female WT (group 1), Female PS19 (group 2),
        % Male WT (group 3), Male PS19 (group 3)

        sets = {{[5 8 9];[2 4 7];1;3;6;10},{[4 7 9];[2 5 8];1;3;6;10},...
            {[2 7 8];[4 5 9];1;3;6;10},{[2 4 5];[7 8 9];1;3;6;10}};
        fullLabels = {'WT-Female','WT-Male','PS19-Female','PS19-Male'};

        for setNum = 1:length(sets)

            % define current set
            curSet = sets{setNum};
            curLabel = fullLabels{setNum};
            dataPooled = cell(nDays,size(curSet,1));

            for dd = 1:nDays
                % combine data
                curDayData = cellfun(@(x) squeeze(x(:,dd,:)),dataByGroup,'UniformOutput',false);
                dataComb = cat(1,curDayData{:})';
                nMouseGroups = cellfun(@(x) size(x,1),dataByGroup);
                groupIdxs =  [0; cumsum(nMouseGroups)];

                % calculate correlations
                corrAll = corr(dataComb,dataComb);

                % pool data by combination set
                for ii = 1:nCombs
                    setA = groupIdxs(combs(ii,1))+1:groupIdxs(combs(ii,1)+1);
                    setB = groupIdxs(combs(ii,2))+1:groupIdxs(combs(ii,2)+1);
                    corrs = corrAll(setA,setB);
                    corrsUnq =  corrs(tril(ones(size(corrs)),-1)==1);

                    % add to current set
                    for st = 1:length(curSet)
                        if ismember(ii,curSet{st})
                            dataPooled{dd,st} = [dataPooled{dd,st};corrs(:)];
                        end
                    end
                end
            end


            %% Plot time coursees

            F = figure;

            setGroups = {1:2,3:6};
            setGroupLabs = {{'Others-Others',['Others-' curLabel]}, fullLabels};
            setColors = {colors,colors4};
            setPShow = {[1 2],[1 2;1 3;2 4; 3 4]};
            setTtls = {'Inter-Group','Intra-Group'};
            tiledlayout(1,2)

            svLabelCur = [ct{:} '-' morphTypes{mt} '_' curLabel];
            xLims = [dayNames(1)-0.5 dayNames(end)+0.5];

            for sg = 1:2
                curSets = setGroups{sg};
                curLabs = setGroupLabs{sg};
                curColors = setColors{sg};
                curPShow = setPShow{sg};

                plotData = dataPooled(:,curSets);
                plotDataCat = cellfun(@(x) cat(2,x{:})', num2cell(plotData,1), 'UniformOutput', false);

                % plot curves
                nexttile(sg); hold on
                [~,hh,~,pCheck] = plotErrorMulti(dayNames,plotDataCat,curLabs,...
                    curColors,[],curPShow,curLabs,1,anovaIdx);

                % label pairwise subplot
                figure(hh)
                sgtitle([setTtls{sg} ': ' svLabelCur])
                xlim(xLims)

                % save pairwise subplot figure
                savefig([svFile '\interMouseCorr_timecourse_sub' setTtls{sg} '_' svLabelCur])

                figure(F)
                xlim(xLims)

                if sg==1
                    % calculate full statistics
                    outStats(end+1,:) = [svLabelCur anovaEffectSizeMinimal(plotDataCat{1}',plotDataCat{2}',nUnits)];
                end
            end

            % set figure title
            sgtitle(['Inter-mouse dfof correlation: ' svLabelCur],'FontSize',16)

            % save figure
            savefig([svFile '\interMouseCorr_timecourse_' svLabelCur])

            % save data
            save([svFile '\data\interMouseCorr_timecourse_' svLabelCur '.mat'],'dayNames','dataPooled')
        end
    end
end
