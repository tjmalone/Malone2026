%% identifyPyrStePlot_area
% Generate indices for pyramidal and stellate cells based on area

clear; clc
close all;

% set directory
foldMan = 'D:\AD_Project\Revision\PyrSteAnalysis';

cd(foldMan);

% load morphology information
load('data_PyrSte_rev.mat','dataPS');
load('morphIdxs.mat','morphIdxs');

% load mutaul index information
load('infoMatch_rev.mat','nCells','infoMatch')
useDays = 1:11;
nDays = length(useDays);
nFOV = length(infoMatch);
useFOV = 1:nFOV;

load('D:\AD_Project\imagingData/foldersLearning.mat')

% load groups
load('D:\AD_Project\imagingData\groupIDs.mat')
for nn = 1:length(groups)
    groups{nn} = intersect(useFOV,groups{nn});
end
nGroups = length(groups);

% set save directories
dataPath = 'Results/';
figPath = 'Results/';

% set activity averaging method
useSFX = 1;
sfx = {'cell','fov','mouse'};

% set sex splitting
useSexes = [1 2 3];
% useSexes = [1];

threshRange = 10;            % threshold for morphology range
threshDist = 4;             % threshold for cell matching. Selected based on visualization
threshOvr = 0.4;            % threshold for matching area overlap

% FOV scale factor (converts to um^2)
sizeScale = 1.47^2/1000^2; % m=m^2/pixel^2

% Load and scale FOV areas
load('D:\AD_Project\imagingData/data_FOV/FOVarea.mat')
areasNorm = FOVareas.Manual*sizeScale;

% whether to average by mouse
byMouse = 0;

% whether to close all figures
closeAll = 0;


%% Calculate overlap percentages
% calculate man, act, com, ste, pyr, man/act, ste/act, pyr/act,
% ste/com, pyr/com

% initialize stats parameters and array
outStats = {};
nUnits = 'FOV';
testName = 'two-tailed unpaired Students t-test';
testPair = 0;
testMC = 0;
testLimitP = 0;

% initialize data structures for current parameter set
dataIdxs = struct();
dataN = struct();
dataRBR = struct();
dataDist = struct();

% loop through FOV
for ff = 1:nFOV
    % current morphology edges
    thLow = dataPS.allThresh(ff)-threshRange;
    thHigh = dataPS.allThresh(ff)+threshRange;

    % current FOV morphology classes
    diamCur = zeros(size(dataPS.allArea{ff}));
    diamCur(dataPS.allArea{ff}<=thLow) = 2;
    diamCur(dataPS.allArea{ff}>=thHigh) = 1;

    % loop through days
    for gg = 1:nDays
        % store cell indices for base categories
        dataIdxs.manual{ff,gg} = 1:nCells.manual(ff);

        dataIdxs.stellate{ff,gg} = find(diamCur==1);
        dataIdxs.pyramidal{ff,gg} = find(diamCur==2);

        % current match information
        if isempty(diamCur)
            curMatch = [];
        else
            curMatch = infoMatch{ff}{gg};
        end

        if ~isempty(curMatch)
            % apply distance threshold
            useIdx = curMatch(:,4)<=threshDist & curMatch(:,5)>=threshOvr;

            % match stellate and pyramidal indices
            isStellate = ismember(curMatch(:,1),dataIdxs.stellate{ff,gg});
            isPyramidal = ismember(curMatch(:,1),dataIdxs.pyramidal{ff,gg});

            % store match indices
            dataIdxs.active{ff,gg} = curMatch(useIdx & (isStellate | isPyramidal),2);
            dataIdxs.common{ff,gg} = curMatch(useIdx & curMatch(:,3)==1,2);
            dataIdxs.Act_Man{ff,gg} = curMatch(useIdx,2);
            dataIdxs.Act_Ste{ff,gg} = curMatch(useIdx & isStellate,2);
            dataIdxs.Act_Pyr{ff,gg} = curMatch(useIdx & isPyramidal,2);
            dataIdxs.Com_Ste{ff,gg} = curMatch(useIdx & isStellate & curMatch(:,3)==1,2);
            dataIdxs.Com_Pyr{ff,gg} = curMatch(useIdx & isPyramidal & curMatch(:,3)==1,2);

        else
            dataIdxs.active{ff,gg} = [];
            dataIdxs.common{ff,gg} = [];
            dataIdxs.Act_Man{ff,gg} = [];
            dataIdxs.Act_Ste{ff,gg} = [];
            dataIdxs.Act_Pyr{ff,gg} = [];
            dataIdxs.Com_Ste{ff,gg} = [];
            dataIdxs.Com_Pyr{ff,gg} = [];
        end
    end
end

% calculate numbers of cells
fNames = fieldnames(dataIdxs);
for mm = 1:length(fNames)
    curNs = cellfun(@length,dataIdxs.(fNames{mm}));
    curNs(curNs==0) = NaN;
    dataN.(fNames{mm}) = curNs;
end

% calculate percentages
perData = struct();
perData.Man = dataN.manual;
perData.ManSte = dataN.stellate;
perData.ManPyr = dataN.pyramidal;
perData.ManNorm = dataN.manual./areasNorm;
perData.ManSteNorm = dataN.stellate./areasNorm;
perData.ManPyrNorm = dataN.pyramidal./areasNorm;
perData.ActMan_Man = dataN.Act_Man./dataN.manual;
perData.ActSte_Ste = dataN.Act_Ste./dataN.stellate;
perData.ActPyr_Pyr = dataN.Act_Pyr./dataN.pyramidal;
perData.ActSte_Act = dataN.Act_Ste./dataN.active;
perData.ActPyr_Act = dataN.Act_Pyr./dataN.active;
S = perData.ActSte_Act;
P = perData.ActPyr_Act;
% perData.ActDif = S./P;
perData.ActDif = P./S;
% perData.ActDif = (S-P)./(P+S);
% perData.ActDif = log(P./S);
% perData.ActDif(perData.ActDif==inf) = NaN;

perData.ComSte_Com = dataN.Com_Ste./dataN.common;
perData.ComPyr_Com = dataN.Com_Pyr./dataN.common;
perData.ComDif = perData.ComPyr_Com./perData.ComSte_Com;
perData.ComDif(abs(perData.ComDif)==inf) = NaN;

% save data
save([dataPath 'area_dataMatch_forFigs.mat'],'dataIdxs','dataN','perData')


%% Plot percentage overlaps

% loop through sexes
for nn = 1:length(useSexes)
    % current sexes
    curSexes = sexes{useSexes(nn)};

    % get field names
    perNames = fieldnames(perData);
    useFields = {1:3,4:6,7:9,10:12,13:15};
    nFieldGroup = length(useFields);

    % define plot information
    plotDay = 1;
    colors = {[0 0 1],[1 0 0]};
    X = {{'ALL','STE','PYR'},{'ALL','STE','PYR'},{'ALL','STE','PYR'},...
        {'STE','PYR','DIFF'},{'STE','PYR','DIFF'}};
    yLabs = {'# cells','# cells per mm^2','Proportion','Proportion','Proportion'};
    yLims = [0 600; 0 1000; 0 1; 0 1.2; 0 1.2];
    ttlBase = {'# Manually Identified Cells','# Manually Identified Cells',...
        'Percent of Manual Population that is Active',...
        'Active Cell Classification','Common Cell Classification'};

    % define p-vaules to plot
    pShow = {[1 4; 2 5; 3 6], [1 4; 2 5; 3 6], [1 4; 2 5; 3 6],...
        [1 2; 1 4; 2 5; 3 6; 4 5], [1 2; 1 4; 2 5; 3 6; 4 5]};

    % initialize figure
    figure('Position',[100,150,1000,750])

    for ff = 1:length(useFields)
        % initialize data array
        nCurField = length(useFields{ff});
        curPerData = cell(nCurField,nGroups);

        % separate data by group and data category
        for gg = 1:nGroups
            for mm = 1:nCurField
                shapeData = perData.(perNames{useFields{ff}(mm)})(intersect(groups{gg},curSexes),plotDay);

                % average by mouse or FOV
                if byMouse==1
                    shapeData = reshape(shapeData,2,[]);
                    shapeData = mean(shapeData,1)';
                    avgType = 'mouse';
                else
                    avgType = 'FOV';
                end

                curPerData{mm,gg} = shapeData;
            end
        end

        % plot in subplot
        subplot(2,3,ff); hold on
        [h,p] = barGroup(X{ff},curPerData,'violin',colors,pShow{ff},'pair');

        % set plot info
        ylabel(yLabs{ff})
        ylim(yLims(ff,:))
        legend(h,groupIDs)
        set(gca,'FontSize',12)

        if ismember(ff,[3 4])
            title([ttlBase{ff} ': p=' num2str(p(1,4),2) ', ' num2str(p(3,6),2)])
        else
            title(ttlBase{ff})
        end

        % calculate full statistics
        if ff==5
            statData = curPerData(1:2,:);
        else
            continue
        end
        testCat = [ttlBase{ff} ': ' sexIDs{nn}];
        outStats(end+1,:) = [testCat ttestEffectSize(...
            statData(:,1),statData(:,2),testName,nUnits,testPair,testMC,testLimitP)];

    end

    % set global figure title
    sgtitle(['Manual vs. Active Cells: (' sexIDs{nn} ')'])
    set(gca,'FontSize',12)

    % save figure
    savefig(gcf,[figPath 'area_dataMatch_forFigs_' sexIDs{nn} '_' avgType '.fig'])

    if closeAll
        close all
    end
end

