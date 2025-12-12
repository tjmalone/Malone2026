%% timecoursePlot_master.m
% Plot neuronal activity variables as a function of day and generate heat
% maps of the relative differences by subgroup
%
% Activity Types:
%   Learning types - learn, recall (not implemented)
%   Activity variables (implemented) - intra-day consistency, inter-day
%       consistency
%   Activity variables (not implemented) - activity correlation to cue
%       template, spatial decoding, cue decoding, activity variation
%   Run types  - all runs, true succes/fail, conditional success/fail
%   Distribution types - dfof, field, to mean RBR, to next RBR, to others
%       RBR
%
% Data Separations:
%   Sex - allSex, female, male
%   Cell type - allType (not implemented), common, grid (not implemented),
%       cue (not implemented)
%   Morphology type - allMorph, ste, pyr
%   Genotype - WT, AD
%
% Heatmap Types: pre-learning, post-learning, learning difference, all days
%


%% Load data

clear; clc; close all

p1 = 'D:\AD_Project\imagingData';
cd(p1)

globalData = struct();
globalFields = {'RBR','inter','dfofRdgBkg','FWidthSum','spatialSelectivity',...
    'decodeIORaw','decodeI','decodeO','decodeAll','dfofSig','speedScore',...
    'dfofInField','dfofNonField','dfofRdg','dfofBkg'};
globalTitles = {'Intra-day consistency (to others)','Inter-day consistency',...
    'Ridge/Background Ratio of \DeltaF/F','Total Field Coverage','Spatial Selectivity',...
    'Decoding IO ratio (raw)','In-Cue Decoding','Out-Cue Decoding','Overall Decoding',...
    'Mean \DeltaF/F_{sig} (while moving)','Speed Score','In-Field \DeltaF/F',...
    'Non-Field \DeltaF/F','\DeltaF/F Ridge','\DeltaF/F Background'};
globalYLabs = {'RBR correlation','inter-day correlation','rdg/bkg ratio',...
    'Field coverage (cm)','Spatial selectivity','IO ratio','% correct bins',...
    '% correct bins', '% correct bins','Mean \DeltaF/F','speed score',...
    '\DeltaF/F (%)','\DeltaF/F (%)','\DeltaF/F (%)','\DeltaF/F (%)'};
nGlobalFields = length(globalFields);

% set load and save folder name
loadFile = [p1 '\Figures\Figures_current\timecourse'];

svFile = [p1 '\Figures\Figures' char(datetime('today','format','yyyyMMdd')) '\timecourseStats'];
if ~isfolder(svFile)
    mkdir(svFile);
    mkdir([svFile '\data']);
end


%% Set input parameters

% cell selection key: sex (ss), cell type (ct), morphology type (mt),
% genotype (gg)
sexIDs = {'allSex','female','male'};
nSex = length(sexIDs);
cellTypes = {'common','grid','nongrid'};
morphTypes = {'allMorph','ste','pyr'};
nMorph = length(morphTypes);
genotypes = {'WT','AD'};

% define averaging method (by cell/by FOV)
avgMethods = 'cell';

% define day categories
dayNames = ['FE' string(1:10)];
dayCats = {1,2:11};
datCatsAll = {1,2:10};

% define colors and patterns
colors6 = {[0 0 1],[1 0 0];[102 105 255]/255,[255 102 102]/255;[0 0 175]/255,[175 0 0]/255};

% close all figures
clAll = 1;


%% Plot activity variables
allOvrIdx = 6:9;

useFields = 9:nGlobalFields;

% loop through data types
for gf = useFields
    % define day categories
    if gf==2
        useDayCats = datCatsAll;
    else
        useDayCats = dayCats;
    end

    % loop through cell types
    for ct = cellTypes
        %% Plot individual figure

        % overwrite decoding cell type
        if ismember(gf,allOvrIdx) && strcmp(ct,'common')
            ct = {'allType'};
        end

        % load data
        svLabelCur = ['\' globalFields{gf} '_' ct{:} '-' avgMethods];
        load([loadFile '\data' svLabelCur '.mat'],'mapData')

        % initialize tiled figure
        figure
        t = tiledlayout(3, 3);
        ttlCur = [globalTitles{gf} ': ' ct{:} ' cells (by ' avgMethods ')'];
        sgtitle(ttlCur)

        % initialize stat data
        statData = struct();
        statData.pAnova = cell(nSex,nMorph);
        statData.pMC = cell(nSex,nMorph);

        % loop through sexes
        for ss = 1:nSex

            % loop through morphology types
            for mt = 1:nMorph
                % skip invalid comparisons
                if ismember(gf,allOvrIdx) && mt~=1
                    continue
                end

                % initialize tile
                nexttile(mt+nSex*(ss-1)); hold on

                % define current data
                curData = squeeze(mapData.timecourse(ss,mt,:));
                curPlotData = cat(1,curData{:});

                % define color
                useColor = colors6(ss,:);

                % plot current data
                [pAnova,pMC] = errorSig(dayNames,curPlotData,useColor,genotypes,useDayCats,0,'RM2W');

                % set labels
                title([sexIDs{ss} '-' morphTypes{mt}])
                xlabel('Session')
                ylabel(globalYLabs{gf})

                % store stat data
                statData.pAnova{ss,mt} = pAnova;
                statData.pMC{ss,mt} = pMC;
            end
        end

        % add correlation data
        statData.pCorr = mapData.corrP;
        statData.rCorr = mapData.corrR;

        % save figure and data
        savefig([svFile svLabelCur '_stats.fig'])
        save([svFile '\data' svLabelCur '_stats.mat'],'statData')

        if clAll
            close all
        end


        %% Plot normalized data

        % initialize tiled figure
        figure
        t = tiledlayout(3, 3);
        ttlCur = [globalTitles{gf} ': ' ct{:} ' cells (by ' avgMethods ')'];
        sgtitle(ttlCur)

        % loop through sexes
        for ss = 1:nSex

            % loop through morphology types
            for mt = 1:nMorph
                % skip invalid comparisons
                if ismember(gf,allOvrIdx) && mt~=1
                    continue
                end

                % initialize tile
                nexttile(mt+nSex*(ss-1)); hold on

                % define current data
                curData = squeeze(mapData.timecourse(ss,mt,:));
                curNormData = cellfun(@(x,y) y/mean(x,'omitnan'),...
                    curData{1},curData{2},'UniformOutput',false);

                % define color
                useColor = colors6(ss,2);

                % plot current data
                normM = cellfun(@(x) mean(x,'omitnan'),curNormData);
                normSEM = cellfun(@(x) nansem(x,1),curNormData);
                errorbar(normM,normSEM)

                % set labels
                title([sexIDs{ss} '-' morphTypes{mt}])
                xlabel('Session')
                ylabel(globalYLabs{gf})
                xticks(1:length(dayNames))
                xticklabels(dayNames)
            end
        end

        % save figure and data
        savefig([svFile svLabelCur '_norm.fig'])

        if clAll
            close all
        end
    end

end

