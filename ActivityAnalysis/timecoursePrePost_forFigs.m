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

% initialize data struct
globalFields = {'RBR','inter','dfofRdgBkg','FWidthSum','spatialSelectivity',...
    'decodeIORaw','dfofSig','speedScore'};
globalTitles = {'Intra-day consistency (to others)','Inter-day consistency',...
    'Ridge/Background Ratio of \DeltaF/F','Total Field Coverage',...
    'Spatial Selectivity','Decoding IO ratio (raw)',...
    'Mean \DeltaF/F_{sig} (while moving)','Speed Score'};
globalYLabs = {'RBR correlation','inter-day correlation','rdg/bkg ratio',...
    'Field coverage (cm)','Spatial selectivity','IO ratio','Mean \DeltaF/F',...
    'speed score'};
nGlobalFields = length(globalFields);

% set load and save folder name
loadFile = [p1 '\Figures\Figures_current\timecourse'];

svFile = [p1 '\Figures\Figures' char(datetime('today','format','yyyyMMdd')) '\timecoursePrePost'];
if ~isfolder(svFile)
    mkdir(svFile);
    mkdir([svFile '\data']);
end


%% Set input parameters

% cell selection key: sex (ss), cell type (ct), morphology type (mt),
% genotype (gg)
sexIDs = {'allSex','female','male'};
nSex = length(sexIDs);
cellTypes = {'common'};
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
colors2 = {[0 0 1],[1 0 0]};

% close all figures
clAll = 0;


%% Plot activity variables

% loop through data types
for gf = 1%:nGlobalFields
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
        if gf==6 && strcmp(ct,'common')
            ct = {'allType'};
        end

        % load data
        svLabelCur = ['\' globalFields{gf} '_' ct{:} '-' avgMethods];
        load([loadFile '\data' svLabelCur '.mat'],'mapData')
        barData = mapData.Diff;


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
                if gf==6 && mt~=1
                    continue
                end

                % initialize tile
                nexttile(mt+nSex*(ss-1))
                title([sexIDs{ss} '-' morphTypes{mt}])

                % select current data
                curBarData = squeeze(barData(ss,mt,:))';

                % plot data
                if ~all(cellfun(@isempty,curBarData))
                    barGroup(genotypes,curBarData,colors2,[1 2])
                end

                % set labels
                ylabel(globalYLabs{gf})
            end
        end

        % save figure and data
        savefig([svFile svLabelCur '_bar.fig'])

        if clAll
            close all
        end
    end

end

