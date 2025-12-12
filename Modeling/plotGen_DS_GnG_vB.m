%% plotGen_DS
% plot downsampling distributions


%% Initialize plotting

clear; close all; clc

% set folder
p1 = 'D:\AD_Project\Behavior';
cd(p1)

load('data/model_vB/modelData_DS_vB.mat')

useVar = 'R2';
barLab = 'R^2';
cellTypes = {'common','grid','nongrid'};
colors = {[0 0 1],[1 0 0.5],[1 0.5 0]};
typeNames = {'Train','Test'};
ttls = {'Training Data Predictions','Test Data Predictions'};
refCType = 1;
xLims = {[0.45 0.7],[0.2 0.6]};

ksRange = 0:0.001:1;

% set save folder name
svFile = [p1 '/Figures/Figures' char(datetime('today','format','yyyyMMdd')) '/DSmodel_vB'];
if ~isfolder(svFile)
    mkdir(svFile);
end


%% Plot all distributions

figure;
tiledlayout(1,2)

% initialize p comparisons
pShow = [1 2;1 3;2 3];
nP = length(pShow);

% initialize minumum p settings
nRange = 1:1000;
pThresh = 0.05;
sigPerThresh = 0.8;
rng(42)

% deinfe plot setting
plotN = 100;
plotItr = randsample(200,plotN);
sigPerAll = cell(2,nP);

% initialize stats parameters and array
outStats = {};
nUnits = 'subsamples';
testName = 'two-tailed unpaired Students t-test';
testPair = 0;
testMC = 0;
testLimitP = 0;

for ii = 1:2
    % get current data
    curType = typeNames{ii};
    typeData = modelData.(useVar).(curType);

    nexttile(); hold on

    % plot distributions
    for jj = 1:3
        refData = typeData{jj}(plotItr);
        [Y,X] = ksdensity(refData,ksRange);
        plot(X,Y,'Color',colors{jj})
    end

    % set plot labels
    xlim(xLims{ii})
    xlabel(barLab)
    ylabel('pdf')
    legend(cellTypes)

    curP = zeros(nP,1);
    curD = zeros(nP,1);
    curN = zeros(nP,1);
    statsSumP = 'p = ';
    statsSumN = 'N = ';
    statsSumD = 'd = ';
    for kk = 1:nP
        curData1 = typeData{pShow(kk,1)};
        curData2 = typeData{pShow(kk,2)};

        % subsample data for plotting
        curDataSub1 = typeData{pShow(kk,1)}(plotItr);
        curDataSub2 =  typeData{pShow(kk,2)}(plotItr);

        % calculate p-value
        [~,curP(kk)] = ttest2(curDataSub1,curDataSub2);

        % calculat means and stanard deviations
        m1    = mean(curData1,'omitnan');
        m2    = mean(curData2,'omitnan');
        s1    = nansem(curData1,1);
        s2    = var(curData2);

        % calculate pooled variance
        sPool = sqrt((s1+s2)/2);

        % perform power analysis
        curSigPer = sampsizepwr('t2',[m1 sPool],m2,[],nRange,'Alpha',pThresh);
        sigPerAll{ii,kk} = curSigPer;

        % find 80% power
        idxPow = find(curSigPer>=sigPerThresh,1,'first');
        if ~isempty(idxPow)
            curN(kk) = nRange(idxPow);
        else
            curN(kk) = NaN;
        end

        % calculate effect size
        effect = meanEffectSize(curData1,curData2,'Effect','cohen');
        curD(kk) = effect{"CohensD","Effect"};

        % update title
        statsSumP = [statsSumP ', ' num2str(curP(kk),2)];
        statsSumN = [statsSumN ', ' num2str(curN(kk),2)];
        statsSumD = [statsSumD ', ' num2str(curD(kk),2)];

        % calculate full statistics
        testCat = [curType ': ' cellTypes{pShow(kk,1)} ' vs. ' cellTypes{pShow(kk,2)}];
        outStats(end+1,:) = [testCat ttestEffectSize(curDataSub1,curDataSub2,testName,nUnits,testPair,testMC,testLimitP)];
    end

    title([ttls{ii} ': ' statsSumP '; ' statsSumN '; ' statsSumD])
end

% save figure
savefig([svFile '/cmpDistribution_GnG.fig'])


%% Plot power analysis

figure
tiledlayout(1,2)
cellTypesCmp = {'common-grid','common-nongrid','grid-nongrid'};


for ii = 1:2
    nexttile(); hold on
   
    h = zeros(1,nP);
    for kk = 1:nP
        % plot current power analysis
        curPower = sigPerAll{ii,kk};
        h(kk) = plot(nRange,curPower,'Color',colors{kk});
    end

    % set labels
    xlabel('N')
    ylabel('power')
    yline(sigPerThresh)
    xscale('log')
    legend(h,cellTypesCmp)
    xlim([1 1000])

end

% save figure
savefig([svFile '/powerAnalysis_GnG.fig'])

