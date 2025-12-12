%% plotInterDatMatrices
% Plot inter-day activity matrics across days for all cell types.
%


%% Load data

clear; clc; close all

p1 = 'D:\AD_Project\imagingData';
cd(p1)

% load inter-day consistency data
load('data/interDayData_matrices.mat','interDay')

% load cell selections
load('data/cellSelect.mat','cellSelect')
cellSelect = cellSelect.learn;

% load alignments
load('foldersLearning.mat','alignsLearning','trueDays')
load('groupIDs.mat')

% set save folder name
svFile = [p1 '/Figures/Figures' char(datetime('today','format','yyyyMMdd')) '/activityExamples'];
if ~isfolder(svFile)
    mkdir(svFile);
    mkdir([svFile '/data']);
end


%% Set input parameters

% cell selection key: sex (ss), cell type (ct), morphology type (mt),
% genotype (gg)
sexIDs = {'allSex','female','male'};
nSex = length(sexIDs);
morphTypes = {'allMorph','ste','pyr'};
nMorph = length(morphTypes);
genotypes = {'WT','AD'};
nGeno = length(genotypes);

% define reference and plot days
plotDays = [1:11];
refDay = 2;         % index within plot days
dayNames = ['FE' string(1:10)];
nDays = length(plotDays);

% define use FOV
useFOV = 1:42;

% close all figures
clAll = 0;


%% Plot sorted activity matrices for all categories

% loop through sexes
for ss = 1:nSex

    % loop through morphology types
    for mt = 1:nMorph

        % initialize tiled figure
        figure
        t = tiledlayout(nGeno,nDays);
        sgtitle([sexIDs{ss} '-' morphTypes{mt}])

        % get current cell selections
        curCellSub = cellSelect.(sexIDs{ss}).common.(morphTypes{mt});

        maxVal = NaN;
        % loop through genotypes
        for gg = 1:nGeno
            %% Process current data

            % get current cell selection and make data copy
            curCellSel = curCellSub.(genotypes{gg});
            curCellSel = idx2common(curCellSel,alignsLearning);

            % define current data
            curData = cellfun(@(x,y) x(y,:),interDay(useFOV,plotDays),...
                curCellSel(useFOV,plotDays),'UniformOutput',false);

            % concatenate data
            curDataCat = cell(1,nDays);
            for dd = 1:nDays
                curDayData = cat(1,curData{:,dd});
                curDayMax = max(curDayData,[],2);
                curDataCat{dd} = curDayData./curDayMax;
                % curDataCat{dd} = curDayData;
            end


            %% Sort activity matrix

            % define reference data
            refData = curDataCat{refDay};

            % calculate reference bins
            sortVar = zeros(size(refData,1),1);
            for kk = 1:size(refData,1)
                [~,sortVar(kk)] = max(refData(kk,:));
            end

            % sort reference bins
            [~,orderMax] = sort(sortVar);

            % sort all day activity
            activityMatCatSort = cell(size(curDataCat));
            for kk = 1:nDays
                activityMatCatSort{kk} = curDataCat{kk}(orderMax,:);
            end


            %% Plot actiivty matrices

            curMax = cellfun(@(x) max(x,[],'all'),activityMatCatSort);
            maxVal = max(maxVal,max(curMax));

            for kk = 1:nDays
                % plot current matrix
                nexttile()
                imagesc(activityMatCatSort{kk})

                % define labels
                ttl = [genotypes{gg} ': Day ' char(dayNames(plotDays(kk)))];

                % add cell number to label
                if kk==1
                    ttl = [ttl ', ' num2str(size(activityMatCatSort{kk},1)) ' cells'];
                end

                % set label
                title(ttl)

                axis('off')
                clim([0 1])
            end

            % collect example traces
            if ss==1 && mt==1 && gg==1
                exampleAct = activityMatCatSort(2:3);
            end

        end

        % for ii = 1:nDays*nGeno
        %     nexttile(ii)
        %     clim([0 maxVal])
        % end
        % save figure
        savefig([svFile '/interday_' sexIDs{ss} '-' morphTypes{mt} '.fig'])

    end
end

