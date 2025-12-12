%% plotExampleSelfCorr
% plots an example joint run distribution demonstrating the difference in
% correlation between success/fail runs and the adjacent next run
%


%% Load data

clear; clc; close all

warning('off')

p1 = 'D:\AD_Project\imagingData';
cd(p1)

% load alignment data
load('foldersLearning.mat')

% load self correlation data
load('data\dataByCell.mat','dataAll')

% load cell selections
load('data\cellSelect.mat','cellSelect')
useCells = cellSelect.learn.allSex.common.allMorph;

% define genotypes
genotypes = {'WT','AD'};
nGeno = length(genotypes);

% set save folder name
svFile = [p1 '\Figures\Figures' char(datetime('today','format','yyyyMMdd')) '\distributionsJoint'];
if ~isfolder(svFile)
    mkdir(svFile);
    mkdir([svFile '\data']);
end


%% Define cue templates

% define cue/reward info
colorsCue = [0 0 0;0.5 0.5 0.5];
colorsRew = [0.75 1 1];
cueX = 2.5:5:1197.5;
runOffset = 600;

% define cue template
cueFolder = 'D:\AnalysisCode\PostAnalysis\Cues\6mEnv2\';
rewLoc = [510 560];

% load cue templates
load([cueFolder 'tempL.mat'])
load([cueFolder 'tempR.mat'])
cueTemp = [tempL,tempR;tempL,tempR];


%% Identify distribution to plot

% calculate success fail differences
allDiff = nan(size(useCells));
for gg = 1:nGeno
    % get current cells
    curCells = useCells.(genotypes{gg});

    % calculate success/fail difference
    curSuccess = cellfun(@(x,y) mean(x(y),'omitnan'),dataAll.selfRun_S,curCells);
    curFail = cellfun(@(x,y) mean(x(y),'omitnan'),dataAll.selfRun_F,curCells);
    curDiff = curSuccess-curFail;

    % combine with global variable
    allDiff = max(allDiff,curDiff,'omitnan');
end

% identify peak difference
[~,idxMax] = max(curDiff,[],'all','omitnan');

% identify FOV
[row,col] = ind2sub(size(curDiff),idxMax);

% row = 42;
% col = 6;
% % 
% row = 39;
% col = 10;

% row = 35;
% col = 5;

% row = 38;
% col = 10;

row = 5;
col = 9;

%% Get example distributions

% move to peak directory
cd(foldersLearning{row}{col})

selfBins = [1:102;121:222];

% initialize plot distribution struct
plotDist = struct();
runTypes = {'Success','Fail'};
nRT = length(runTypes);

for rt = 1:nRT
    % load joint activity
    load(['successFailAnalysis\' runTypes{rt} '\jointActivity.mat'],'jointActivity');
    curDfofM = jointActivity.dfofMInterpM(alignsLearning{row}(:,col));
    mxSz = max(cellfun(@(x) size(x,1),curDfofM),[],'omitnan');

    % fill missing cell data
    curDfofM(cellfun(@isempty,curDfofM)) = {nan(mxSz,240)};

    % perform run-by-run correlation for all runs
    curRunCorrAll = cellfun(@(x) corr(x(:,selfBins(1,:))',x(:,selfBins(2,:))',...
        'rows','pairwise'),curDfofM,'UniformOutput',false);

    % take mean of self correlations
    curRunCorrDiag = cellfun(@(x) diag(x),curRunCorrAll,'UniformOutput',false);
    curRunCorrCat = cat(2,curRunCorrDiag{:});
    curRunCorrMean = mean(curRunCorrCat,2,'omitnan');

    % select median lap
    if rt==1
        curMed = max(curRunCorrMean);
    else
        curMed = median(curRunCorrMean);
    end
    [~,curMedIdx] = min(abs(curRunCorrMean-curMed));

    % store distribution for selected lap
    curDfofMCat = cat(3,curDfofM{:});
    plotDist.(runTypes{rt}).dist = squeeze(curDfofMCat(curMedIdx,:,:))';
    plotDist.(runTypes{rt}).lapCorr = curMed;
    
    load(['successFailAnalysis\' runTypes{rt} '\useRuns.mat'],'useRunIdxAll');
    plotDist.(runTypes{rt}).lapIdx = useRunIdxAll(curMedIdx);
end

cd(p1)


%% Plot example distribution

% set base plot labels
ylabBase = {'\DeltaF/F'};
xBase = linspace(2.5,1197.5,240);
nBins = 120;

% define colors
colors2 = {[0 0 1],[1 0 0]};
colors4 = {[0 0 1],[1 0 0];[0.5 0.5 1],[1 0.5 0.5]};

% initialize figure
figure;
tiledlayout(2,2)

% plot joint runs
nexttile(1,[1 2]); hold on
for rt = 1:nRT
    curPlotData = plotDist.(runTypes{rt}).dist;
    semshade(curPlotData,0.3,colors2{rt},xBase);
end

% plot individual run types overlayed
for rt = 1:nRT
    nexttile(rt+2); hold on
    curPlotData = plotDist.(runTypes{rt}).dist;

    % plot individual lap average
    for ii = 1:2
        subRange = (nBins*(ii-1)+1):(nBins*(ii-1)+nBins);
        semshade(curPlotData(:,subRange),0.3,colors4{ii,rt},xBase(1:nBins));
    end

    curLap = plotDist.(runTypes{rt}).lapIdx;
    title([runTypes{rt} ' run example: laps ' num2str(curLap) '-' num2str(curLap+1)...
        ', r=' num2str(plotDist.(runTypes{rt}).lapCorr,2)])
end

axIdx = [1 3 4];
axLim = [240, 120, 120];
for ax = 1:3
    nexttile(axIdx(ax)); hold on
    ylim([0 0.25])

    % plot cues
    for ii = 1:size(cueTemp,2)
        plotCues(cueX(1:axLim(ax)),cueTemp(1:axLim(ax),ii),max(ylim),colorsCue(ii,:),min(ylim));
    end

    % plot rewards
    pos = [rewLoc(1) min(ylim) diff(rewLoc) diff(ylim)];
    rectangle('Position',pos,'FaceColor',colorsRew,...
        'EdgeColor','none')
    if ax==1
        pos = [rewLoc(1)+runOffset min(ylim) diff(rewLoc) diff(ylim)];
        rectangle('Position',pos,'FaceColor',colorsRew,...
            'EdgeColor','none')
    end

    % flip children order to plot reward behind distributions
    axCur = gca;
    axCur.Children = flipud(axCur.Children);
end

savefig([svFile '\selfCorr_example.fig'])
