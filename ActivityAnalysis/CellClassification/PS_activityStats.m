%% PS_activityStats.m
% Calculate and plot the cue scores, cell type percentages, and other
% activity statistics for stellate/pyramidal cells.


%% Set input parameters

clear; clc
close all;

% set directory
foldMan = 'D:\AD_Project\imagingData\analysis_ManualSelection';
cd(foldMan);

% load alignment and cell type data
load('D:\AD_Project\imagingData\foldersLearning.mat')
load('D:\AD_Project\imagingData\data\globalCellTypes.mat')

nFOV = length(alignsLearning);
useFOV = 1:nFOV;

% set day groups
useDays = 1:11;
nDays = length(useDays);

% set save directories
dataPath = 'figures\activity\data\';
figPath = 'figures\activity\';

if ~isfolder(dataPath)
    mkdir(dataPath)
    mkdir(figPath)
end

% set sex splitting
useSexes = [1 2 3];

% set cue template types
cueTypes = {'Left','Right','All'};
% cueTypes = {'Left','Right'};
Xtypes = {'Left Cue','Right Cue','All Cue','Grid','Other'};

% whether to use all cells (1) or common cells (2)
cellSetType = 2;
cellSetLabel = {'all','common'};

% load morphology data
load(['D:\AD_Project\imagingData\analysis_ManualSelection\figures\match\data\'...
    'area_dataMatch_C158_R20_d4_o0-4_RBR_sigMean_common.mat','dataIdxs']);
morphTypes = {'ste','pyr'};

% load morphology cells
if cellSetType==1
    useCells{1} = dataIdxs.Act_Ste(:,:);
    useCells{2} = dataIdxs.Act_Pyr(:,:);
elseif cellSetType==2
    useCells{1} = dataIdxs.Com_Ste(:,1);
    useCells{2} = dataIdxs.Com_Pyr(:,1);
end

binWidth = 5;

% rerun all calculations
newRun = 0;


%% Get activity information for all mice

if newRun || ~isfile([dataPath 'cellTypes_' cellSetLabel{cellSetType} '.mat'])
    cellTypes = struct();

    % loop through morphology types
    for mm = 1:length(morphTypes)
        cellTypes.cueScore.(morphTypes{mm}) = cell(nFOV,nDays);
        cellTypes.SI.(morphTypes{mm}) = cell(nFOV,nDays);
        cellTypes.corrInfo.(morphTypes{mm}) = cell(nFOV,nDays);
        cellTypes.corrDrop.(morphTypes{mm}) = cell(nFOV,nDays);
        cellTypes.Local.(morphTypes{mm}) = cell(nFOV,nDays);
        cellTypes.Global.(morphTypes{mm}) = cell(nFOV,nDays);

        % loop through FOV
        for ii = 1:nFOV
            disp(ii)

            % loop through days
            for jj = 1:nDays
                % correct for "true" day
                trueInv = find(trueDays(ii,:)==jj);

                if isempty(trueInv)
                    continue
                end

                cd(foldersLearning{ii}{trueInv})


                %% Identify current cells

                % define use cells
                if cellSetType==2
                    % calculate all day common cell indices based on day 1
                    FOVCells = ismember(alignsLearning{ii}(:,1),useCells{mm}{ii});
                    curCells = alignsLearning{ii}(FOVCells,trueInv);

                elseif cellSetType==1
                    % use all cell indices for each day
                    FOVCells = useCells{mm}{ii,jj};
                    curCells = FOVCells;
                end

                if isempty(curCells)
                    continue
                end

                % number of matched cells
                nCurCells = length(curCells);


                %% Get cue scores

                % get true cell number
                roiSize = matfile('allRois.mat');
                nCells = size(roiSize,'roi',3);

                % extract cue score of each type
                allScores = nan(nCells,1);
                for kk = 1:length(cueTypes)
                    load(['cueAnalysis_sig\newScoreShuffleTemplate\' cueTypes{kk} '\cueCells.mat'])
                    allScores(cueCells.useIdx,kk) = cueCells.realScores;
                end

                % store cue score subset
                curScores = allScores(curCells,:);
                cellTypes.cueScore.(morphTypes{mm}){ii,jj} = curScores;
                cellTypes.cueScore.X = cueTypes;


                %% Get spatial information

                % extract SI information using old and new method and
                % spatial selectivity
                SIcomb = zeros(nCells,3);

                load('spatialInfo_sig\SI.mat')
                SIcomb(:,1) = SI;

                load('spatialInfo_sig\SInew.mat')
                SIcomb(:,2) = SI;

                load('spatialInfo_sig\spatialSelectivity.mat')
                SIcomb(:,3) = SS;

                % store cue score subset
                curSI = SIcomb(curCells,:);
                cellTypes.SI.(morphTypes{mm}){ii,jj} = curSI;
                cellTypes.SI.X = {'SI-old','SI-new','Spatial Selectivity'};


                %% Get cue correlation information

                % load binned activity and cue template
                d = dir('dfofaveragesmooth_sig*');
                load(d(1).name)
                load('cueAnalysis_sig\tempRL.mat')

                % exapand cues by x pixels
                x_pixels = 2;
                se = strel('line', 2 * x_pixels + 1, 0);
                curTemp  = imdilate(tempRL',se)';

                % get current activity data
                % curDfof = dfofaveragesmooth_sig(:,curCells);
                curDfof = mean(dfofaveragesmooth_sig(:,curCells),2);
                nCorrCells = size(curDfof,2);

                % calculate no lag cue correlation
                corrNoLags = corr(curDfof,tempRL);

                % calculate lag cue correlation
                maxLag = 20;
                corrLags = nan(nCorrCells,1);
                lagsRaw = nan(nCorrCells,1);
                lagsAbs = nan(nCorrCells,1);
                for kk = 1:nCorrCells
                    % perform cross correlation
                    [thisCorr,thisLags] = xcorr(curDfof(:,kk),tempRL,'normalized');

                    % get peak correlation
                    [mx,mxIdx] = max(thisCorr);
                    corrLags(kk) = mx;
                    lagsRaw(kk) = thisLags(mxIdx)*binWidth;
                    lagsAbs(kk) = abs(thisLags(mxIdx)*binWidth);
                end

                % store cue score subset
                curCorrInfo = [corrNoLags corrLags lagsRaw lagsAbs];
                cellTypes.corrInfo.(morphTypes{mm}){ii,jj} = curCorrInfo;
                cellTypes.corrInfo.X = {'No Lag Corr','Peak Corr','Peak Lag (raw)','Peak Lag (abs)'};


                %% Calculate dfof dropoff by distance

                % Find the start and end indices of each cluster
                clusterStart = strfind([0, tempRL'], [0 1]);
                clusterEnd = strfind([tempRL', 0], [1 0]);

                % Calculate the center index of each cluster
                clusterCenters = round((clusterStart + clusterEnd+1) / 2);

                % Pre-allocate a distance array
                distancesToCenters = zeros(1, length(tempRL));

                % Calculate the distance of each index to the nearest center
                for i = 1:length(tempRL)
                    distancesToCenters(i) = min(abs(clusterCenters - i));
                end

                actToCenters = {};
                maxDistance = 16; % true max for NE
                for i = 1:maxDistance
                    actToCenters{i} = mean(curDfof(distancesToCenters==i-1),'omitnan');
                end
                actMeans = cellfun(@(x) mean(x,'omitnan'),actToCenters);

                cellTypes.corrDrop.(morphTypes{mm}){ii,jj} = actMeans;
                cellTypes.corrDrop.X = {'Correlation dropoff'};


                %% Get local cell types

                % get local cue cell indices for left and right
                load('cueAnalysis_sig\newScoreShuffleTemplate\cueCellsAll.mat')
                subCueIdx = cell(1,2);
                for kk = 1:length(cueTypes)
                    subCueIdx{kk} = cueCellsAll.(cueTypes{kk}).cueCellRealIdx;
                end
                allCueIdx = cat(1,subCueIdx{:});

                % get local grid cell indices
                load('gridAnalysis_sig\gridIdx.mat')
                allGridIdx = setdiff(gridIdx,allCueIdx);

                % store local cells types
                typesLocal = false(length(curCells),length(Xtypes));
                for kk = 1:length(cueTypes)
                    typesLocal(ismember(curCells,subCueIdx{kk}),kk) = true;
                end
                typesLocal(ismember(curCells,allGridIdx),length(cueTypes)+1) = true;
                typesLocal(~any(typesLocal,2),length(Xtypes)) = true;

                % store values
                cellTypes.Local.(morphTypes{mm}){ii,jj} = typesLocal;
                cellTypes.Local.X = Xtypes;


                %% Get global cell types (common cells only)

                if cellSetType==1
                    typesGlobal = NaN;
                elseif cellSetType==2
                    % load global cell types
                    % load('commonCellTypes_sig.mat')

                    % store global cells types
                    typesGlobal = false(length(curCells),length(Xtypes));
                    for kk = 1:length(cueTypes)+1
                        curGlobalType = [commonTypes{kk} '_sig'];
                        typesGlobal(ismember(curCells,...
                            globalCellTypeIdx.(curGlobalType){ii}(:,trueInv)),kk) = true;
                    end
                    typesGlobal(~any(typesGlobal,2),length(Xtypes)) = true;
                end

                % store values
                cellTypes.Global.(morphTypes{mm}){ii,jj} = typesGlobal;
                cellTypes.Global.X = Xtypes;
            end
        end
    end

    cd(foldMan)
    save([dataPath 'cellTypes_' cellSetLabel{cellSetType} '.mat'],'cellTypes')
else
    load([dataPath 'cellTypes_' cellSetLabel{cellSetType} '.mat'],'cellTypes')
end


%% Process and concatenate data

useSFX = 1;
sfx = {'cell','FOV',};

% load genotypes {[WT],[AD]}
load('D:\AD_Project\imagingData\groupIDs.mat')
nGroups = length(groups);

fNames = fieldnames(cellTypes);
nCat = 6;
firstFOVAvg = 3;

close all

for ss = 1:length(sexes)
    % current sexes
    curSexes = sexes{useSexes(ss)};

    for ff = 1:nCat
        for mm = 1:length(morphTypes)
            curData1 = cellTypes.(fNames{ff}).(morphTypes{mm});

            emptyCell = cellfun(@isempty,curData1);
            firstData = find(emptyCell==0,1);
            nColumns = size(curData1{firstData},2);
            % fix empty cells
            for ii = 1:nFOV
                for jj = 1:nDays
                    if emptyCell(ii,jj)==0; continue; end

                    if emptyCell(ii,1)==0
                        curData1{ii,jj} = nan(size(curData1{ii,jj-1}));
                    else
                        curData1{ii,jj} = nan(1,nColumns);
                    end
                end
            end

            if useSFX==1 && ff<firstFOVAvg
                % average by cue score by cell if set
                curData2 = curData1;
            else
                % always average percents by FOV
                curData2 = cellfun(@(x) mean(x,1,'omitnan'),curData1,'UniformOutput',false);
            end

            curData3 = cell(nGroups,nDays);
            % concatnate FOV
            for gg = 1:nGroups
                for  jj = 1:nDays
                    curData3{gg,jj} = cat(1,curData2{intersect(groups{gg},curSexes),jj});
                end
            end

            % store concatenated values
            cellTypes.([fNames{ff} '_cat']).(morphTypes{mm}) = curData3;
            cellTypes.([fNames{ff} '_cat']).X = cellTypes.(fNames{ff}).X;
        end
    end


    %% Initialize plotting

    % define fields and labels per use field
    fNames = fieldnames(cellTypes);
    useFields = nCat+1:2*nCat;
    svLabelField = fNames(useFields);

    % set line fields
    lineFields = 4;
    dayFields = {2,{[3]}};

    ylabs = {'cue score','SI','correlation','dfof dropoff','Fraction of cells','Fraction of cells'};
    ttlsFields = {'Cue Score','Spatial Information','Correlations','dfof Dropoff'...
        'Cell type percentages (local)','Cell type percentages (global)'};

    % legemnd info
    legs = {'WT-Stellate','AD-Stellate','WT-Pyramidal','AD-Pyramidal'};
    legsShort = {'WS','AS','WP','AP'};

    % define day categories and labels
    % dayCats = 8:11;
    dayCats = 2:11;
    dayNames = ['FE' string(1:10),'rFam','rNov'];

    % set plotting paramters
    colors4 = {[0 0 1],[1 0 0],[0.3 0.3 1],[1 0.3 0.3]};
    styles4 = {'-','-','-.','-.'};
    Xend = 15;
    X = ((1:Xend)-1)*5;

    % close all figures
    clAll = 0;

    % cycle through fields
    for ff = 1:length(useFields)
        %% Split and format data

        curData = cellTypes.([fNames{useFields(ff)}]);

        % cycle through dayCats
        for cc = 1%:length(dayCats)

            % set day category labels
            curDayCats = dayCats;
            nG = size(curData.(morphTypes{1}){1,curDayCats(1)},2);

            % initialize data array ( x group)
            data = cell(nG,nGroups*length(morphTypes));
            curDataStore = cell(1,nGroups*length(morphTypes));

            dayLabels = strcat({'Days '},dayNames(curDayCats(1)),'-',dayNames(curDayCats(end)));

            % split data by genotype morphology and day category
            % loop through morphology types
            for mm = 1:2
                for gg = 1:nGroups

                    curData1 =  [];
                    for kk = 1:nDays
                        curData1(:,:,kk) = curData.(morphTypes{mm}){gg,kk};
                    end
                    curDataStore{(mm-1)*2+gg} = curData1;
                    curData2 = mean(curData1(:,:,curDayCats),3,'omitnan');

                    % store data
                    if ~ismember(ff,lineFields)
                        for kk = 1:nG
                            data{kk,(mm-1)*2+gg} = curData2(:,kk);
                        end
                    else
                        data{1,(mm-1)*2+gg} = curData2';
                    end
                end
            end


            %% Plot data

            figure; hold on

            if ~ismember(ff,lineFields)
                pShow = generatePShow(nG,[1 2;1 3;2 4;3 4]);
                [~,p] = barGroup(curData.X,data,colors4,pShow);

                legend(legs)
            else
                pShow = [1 3;2 4;1 2;3 4];

                % intialize combined figure
                endCut = 2;
                data = data(1,:);
                data = cellfun(@(x) x(1:Xend-endCut,:),data,'UniformOutput',false);

                % dataLine = data;
                % dataLine = cellfun(@(x) x./x(1,:),data,'UniformOutput',false);
                dataLine = cellfun(@(x) x/mean(x(1:2,:),'all','omitnan'),data,'UniformOutput',false);
                plotErrorMulti(X(1:Xend-endCut),dataLine,legs,colors4,styles4,pShow,legsShort,0)

                % set labels
                xlabel('distance to nearest cue (cm)')
            end

            % set labals
            ylabel(ylabs{ff})
            if ff<firstFOVAvg
                title([ttlsFields{ff} ': ' char(dayLabels) ' (by ' sfx{useSFX} ', ' sexIDs{ss} ')'])
            else
                title([ttlsFields{ff} ': ' char(dayLabels) ' (by ' sfx{2} ', ' sexIDs{ss} ')'])
            end
            set(gca,'FontSize',12)

            % save figure
            savefig([figPath 'actInfo-' svLabelField{ff} '_' cellSetLabel{cellSetType} '_' sfx{useSFX} '_' sexIDs{ss} '.fig'])

            % save data
            save([dataPath 'actInfo-' svLabelField{ff} '_' cellSetLabel{cellSetType} '_' sfx{useSFX} '_' sexIDs{ss} '.mat'],'data')

            %% Plot statistics by day
            
            % check if field has day plot variables
            [logi,indi] = ismember(ff,dayFields{1});
            if logi
                % get internal day plot variables
                plotIdx = dayFields{2}{indi};

                % loop through day plot variables
                for ii = plotIdx
                    curDayData = cellfun(@(x) squeeze(x(:,ii,:))',curDataStore,'UniformOutput',false);

                    figure; hold on
                    pShow = [1 3;2 4;1 2;3 4];
                    plotErrorMulti(1:nDays,curDayData,legs,colors4,styles4,pShow,legsShort,0)
                    title('Spatial Selectivity')
                    xlabel('Day')
                    ylabel('Spatial Selectivity')

                    % save figure
                    savefig([figPath 'actInfo-spSelectivityDay_' cellSetLabel{cellSetType} '_' sfx{useSFX} '_' sexIDs{ss} '.fig'])

                    % save data
                    save([dataPath 'actInfo-spSelectivityDay_' cellSetLabel{cellSetType} '_' sfx{useSFX} '_' sexIDs{ss} '.mat'],'data')

                end
            end

            %% Plot ridge background ratio

            if ismember(ff,lineFields)
                rdgRange = 1:3;
                bkgRange = 7:9;

                % calculate ridge background ratio
                ratioRB = cell(1,length(legs));
                for ii = 1:length(legs)
                    dataRB = dataLine{ii};
                    rdgM = mean(dataRB(rdgRange,:),'omitnan');
                    bkgM = mean(dataRB(bkgRange,:),'omitnan');
                    ratioRB{ii} = rdgM./bkgM;
                end

                % plot ridge background ratio
                figure; hold on
                pShow = generatePShow(1,[1 2;1 3;2 4;3 4]);
                [~,p] = barGroup(legs,ratioRB,colors4,pShow);

                ylabel('Ridge/Background Ratio')
                title(['Ridge/Background Ratio: ' char(dayLabels)...
                    ' (average days -> RB calc) (by ' sfx{useSFX} ', ' sexIDs{ss} ')'])

                % save figure
                savefig([figPath 'actInfo-RB-m1_' cellSetLabel{cellSetType} '_' sfx{useSFX} '_' sexIDs{ss} '.fig'])

                % save data
                save([dataPath 'actInfo-RB-m1_' cellSetLabel{cellSetType} '_' sfx{useSFX} '_' sexIDs{ss} '.mat'],'data')


                %% method 2

                rdgRange = 1:3;
                bkgRange = 7:9;

                % calculate ridge background ratio
                ratioRB = cell(1,length(legs));
                for ii = 1:length(legs)
                    dataRB = curDataStore{ii}(:,:,curDayCats);
                    rdgM = squeeze(mean(dataRB(:,rdgRange,:),2,'omitnan'));
                    bkgM = squeeze(mean(dataRB(:,bkgRange,:),2,'omitnan'));
                    ratioRB{ii} = mean(rdgM./bkgM,2);
                end

                % plot ridge background ratio
                figure; hold on
                pShow = generatePShow(1,[1 2;1 3;2 4;3 4]);
                [~,p] = barGroup(legs,ratioRB,colors4,pShow);

                ylabel('Ridge/Background Ratio')
                title(['Ridge/Background Ratio: ' char(dayLabels)...
                    ' (RB calc -> avgerage days) (by ' sfx{useSFX} ', ' sexIDs{ss} ')'])

                % save figure
                savefig([figPath 'actInfo-RB-m2_' cellSetLabel{cellSetType} '_' sfx{useSFX} '_' sexIDs{ss} '.fig'])

                % save data
                save([dataPath 'actInfo-RB-m2_' cellSetLabel{cellSetType} '_' sfx{useSFX} '_' sexIDs{ss} '.mat'],'data')
            end

        end

        % close figure
        if clAll
            close all
        end
    end
end