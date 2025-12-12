%% PS_plotDistributions.m
% Calculate and plot neuronal activity distributions as a function of track
% location for stellate/pyramidal cells.


%% Set input parameters

clear; clc
close all;

% set directory
foldMan = 'D:\AD_Project\imagingData\analysis_ManualSelection';

cd(foldMan);

% load alignment data
load('D:\AD_Project\imagingData\foldersLearning.mat')

nFOV = length(alignsLearning);
useFOV = 1:nFOV;

% set day groups
useDays = 1:11;
nDays = length(useDays);

% set save directories
dataPath = 'figures\activity\data\';
figPath = 'figures\activity\';

mkdir(dataPath)
mkdir(figPath)

% set sex splitting
useSexes = [1 2 3];

% whether to use all cells (1) or common cells (2)
cellType = 2;
cellTypeLabel = {'all','common'};


%% Load cell morphology info

% base code
load('D:\AD_Project\imagingData\analysisPyrSte\match_data\area_dataMatch_C158_R20_d4_o0-4_RBR_sigMean_common.mat','dataIdxs');

if cellType==1
    useCells{1} = dataIdxs.Act_Ste(:,:);
    useCells{2} = dataIdxs.Act_Pyr(:,:);
elseif cellType==2
    useCells{1} = dataIdxs.Com_Ste(:,1);
    useCells{2} = dataIdxs.Com_Pyr(:,1);
end


%% Find activity distributions for all mice

if ~isfile('activity_data\distActivity.mat')
    distActivity = struct();
    distActivity.dfof = cell(nFOV,nDays);
    distActivity.field = cell(nFOV,nDays);
    distActivity.RBR_mean = cell(nFOV,nDays);
    distActivity.RBR_others = cell(nFOV,nDays);
    distActivity.RBR_next = cell(nFOV,nDays);

    % cycle through FOV
    for ii = 1:nFOV
        disp(ii)

        % cycle through days
        for jj = 1:nDays
            cd(foldersLearning{ii}{jj})

            %% Store spatial dfof

            d = dir('dfofaveragesmooth_sig*cells.mat');
            load(d(1).name)
            distActivity.dfof{ii,jj} = dfofaveragesmooth_sig';


            %% Store spatial field distribution (percent of cells with field)

            % idenitfy field locations of all cells
            if isfile('gridAnalysis_sig\allCellsFields.mat')
                % load field locations
                load('gridAnalysis_sig\allCellsFields.mat')
            else
                % find field locations
                load('gridAnalysis_sig\allCellsCorrected.mat')

                % save field locations
                allCellsFields = allCellsCorrected.dfofaveragesmoothFields';
                save('gridAnalysis_sig\allCellsFields.mat','allCellsFields')
            end

            % threshold fields
            allCellsFields(allCellsFields~=0) = 1;

            % store field distribution
            distActivity.field{ii,jj} = allCellsFields;


            %% Store RBR consistency distribtuion

            load('RunByRun_sig\corrInfoLocation.mat')
            distActivity.RBR_mean{ii,jj} = corrInfoLocation.toMeanMean;
            distActivity.RBR_others{ii,jj} = corrInfoLocation.toOthersMean;
            distActivity.RBR_next{ii,jj} = corrInfoLocation.toNextMean;

        end
    end

    cd(foldMan)
    save('activity_data\distActivity.mat','distActivity')
else
    load('activity_data\distActivity.mat','distActivity')
end


%% Initialize plotting

% load genotypes {[WT],[AD]}
load('D:\AD_Project\imagingData\groupIDs.mat')
nGroups = length(groups);

fields = fieldnames(distActivity);

% set save folder name
svFile = 'match_data';


%% Load cues

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


%% Plot distributions

close all

% define fields and labels per use field
useFields = [1 2 4];
% useFields = [1];

ylabs = {'\DeltaF/F','Percent of cells with field','RBR correlation',...
    'RBR correlation','RBR correlation'};
svLabelField = fields(useFields);
X = {linspace(2.5,597.5,120),linspace(2.5,597.5,120),linspace(12.5,587.5,116),...
    linspace(12.5,587.5,116),linspace(12.5,587.5,116)};
ttlsFields = {'\DeltaF/F Distribution','Field Distribution','RBR Consistiency Distribution (to Mean)'...
    'RBR Consistiency Distribution (to Others)','RBR Consistiency Distribution (to Next)'};
yLimits = {[0 0.2], [0 0.2], [-0.1 0.1];...
    [0 0.5], [0 0.5], [-0.3 0.2];...
    [0.3 0.9], [0.3 0.9], [-0.2 0.5];...
    [0 0.7], [0 0.7], [-0.2 0.4];...
    [0.1 0.8], [0.1 0.8], [-0.4 0.6]};

% define day categories and labels
dayCats = {2:3,8:11};
ttlsDayCats = {'Pre-Learning','Post-Learning'};
svDayCats = {'NE-Pre','NE-Post'};
ttlsMorph = {'Stellate','Pyramidal'};

% define averaging labels
useSFX = 2;
sfx = {'cell','FOV'};
if cellType==1
    useSFX = 2;
    disp('All cells can only be analyzed by FOV')
end

% define general labels
xlabs = 'Track Location';
dayNames = ['FE' string(1:10),'rFam','rNov'];

% set plotting paramters
colors = {[0.5 0.5 1],[1 0.5 0.5];[0 0 1],[1 0 0]};
cueTypes = {'Left','Right'};

% close all figures
clAll = 0;

% loop through fields
for ff = 1:length(useFields)
    for ss = 1:length(sexes)
        % current sexes
        curSexes = sexes{useSexes(ss)};

        %% Define data distributions

        % initialize data array (morphology x day category x group)
        data = cell(2,length(dayCats),nGroups);

        % loop through morphology types
        for mm = 1:2
            % current distribution
            curField = cell(nFOV,nDays);
            for ii = 1:nFOV
                if cellType==2
                    % calculate all day common cell indices based on day 1
                    FOVCells = ismember(alignsLearning{ii}(:,1),useCells{mm}{ii});
                    if all(FOVCells==0)
                        continue
                    end

                    for jj = 1:nDays
                        % correct for "true" day
                        trueInv = find(trueDays(ii,:)==jj);

                        if isempty(trueInv)
                            curField{ii,jj} = nan(sum(FOVCells),length(X{ff}));
                        else
                            curCells = alignsLearning{ii}(FOVCells,trueInv);
                            curField{ii,jj} = distActivity.(fields{useFields(ff)}){ii,trueInv}(curCells,:);
                        end
                    end
                elseif cellType==1
                    for jj = 1:nDays
                        % use all cell indices for each day
                        FOVCells = useCells{mm}{ii,jj};

                        % correct for "true" day
                        trueInv = find(trueDays(ii,:)==jj);

                        if isempty(trueInv)
                            curField{ii,jj} = nan(length(FOVCells),length(X{ff}));
                        else
                            curField{ii,jj} = distActivity.(fields{useFields(ff)}){ii,trueInv}(FOVCells,:);
                        end
                    end

                end

                % cycle through dayCats
                for cc = 1:length(dayCats)

                    % set day category labels
                    curDayCats = dayCats{cc};

                    dayLabels = strcat({'Days '},dayNames(curDayCats(1)),'-',dayNames(curDayCats(end)));

                    % split data by genotype and day category
                    for gg = 1:nGroups
                        curMice = intersect(groups{gg},curSexes);
                        curData1 =  cell(length(curMice),length(curDayCats));
                        for kk = 1:length(curDayCats)
                            curData1(:,kk) = curField(curMice,curDayCats(kk));
                        end

                        % average by FOV
                        if useSFX==1
                            curData2 = curData1;
                        elseif useSFX==2
                            curData2 = cellfun(@(x) mean(x,1,'omitnan'),curData1,'UniformOutput',false);
                        end

                        curData3 = [];
                        % concatnate FOV
                        for  kk = 1:length(curDayCats)
                            curEmpty = cellfun(@isempty,curData2(:,kk));
                            curData3(:,:,kk) = cat(1,curData2{~curEmpty,kk});
                        end

                        data{mm,cc,gg} = mean(curData3,3,'omitnan');
                    end
                end
            end
        end


        %% Plot data

        % loop through plot types
        for ii = 1:3
            figure; hold on
            sgtitle([ttlsFields{useFields(ff)} ' (by ' sfx{useSFX} ', ' sexIDs{ss} ')'])

            for gg = 1:nGroups
                for jj = 1:2
                    %% Define data

                    plotData = cell(1,2);
                    for kk = 1:2
                        if ii==1
                            % stellate vs pyramidal
                            plotData{kk} = data{kk,jj,gg};
                            ttl = [groupIDs{gg} ' ' ttlsDayCats{jj} ': Ste vs. Pyr'];
                        elseif ii==2
                            plotData{kk} = data{jj,kk,gg};
                            ttl = [groupIDs{gg} ' ' ttlsMorph{jj} ': Pre vs. Post'];
                        elseif ii==3
                            plotData{kk} = data{jj,2,gg}-data{jj,1,gg};
                            ttl = [groupIDs{gg} ' ' ttlsMorph{jj} ': Learning Difference'];
                        end
                    end

                    %% Plot cues and rewards

                    % plot data
                    subplot(2,2,(gg-1)*2+jj); hold on
                    title(num2str((gg-1)*2+jj))

                    % calculate cue level
                    cueLvl = yLimits{useFields(ff),ii}(2);

                    % plot rewards
                    pos = [rewLoc(1) 0 diff(rewLoc) cueLvl];
                    rectangle('Position',pos,'FaceColor',colorsRew,'EdgeColor','none')

                    legs = {};
                    h = [];
                    idx = 1;
                    % plot cues
                    for jj = 1:size(cueTemp,2)
                        h(idx) = plotCues(cueX,cueTemp(:,jj),cueLvl,colorsCue(jj,:));
                        legs{idx} = cueTypes{jj};
                        idx = idx+1;
                    end


                    %% Plot data

                    % plot distributions as shaded lines
                    legs = {};
                    h = [];
                    idx = 1;
                    for kk = 1:2
                        h(idx) = semshade(plotData{kk},0.3,colors{kk,gg},X{useFields(ff)});

                        % set legend
                        if ii==1
                            legs{idx} = ttlsMorph{kk};
                        elseif ii==2
                            legs{idx} = ttlsDayCats{kk};
                        elseif ii==3
                            legs{idx} = 'Learning Difference';
                            continue
                        end
                        idx = idx+1;
                    end

                    % bring axes to top
                    set(gca,'Layer','top')

                    % set labels and axes limits
                    legend(h,legs,'Location','southwest')
                    xlabel(xlabs)
                    ylabel(ylabs{useFields(ff)})
                    title(ttl)
                    ylim(yLimits{useFields(ff),ii})

                end
            end

            % save figure
            savefig([figPath 'actDist-' svLabelField{ff} '_' sfx{useSFX}...
                '_' cellTypeLabel{cellType} '_' sexIDs{ss} '_' num2str(ii) '.fig'])

            % save data
            save([dataPath 'actDist-' svLabelField{ff} '_' sfx{useSFX}...
                '_' cellTypeLabel{cellType} '_' sexIDs{ss} '_' num2str(ii) '.mat'],'data')
        end


        %%

        % close figure
        if clAll
            close
        end
    end
end
