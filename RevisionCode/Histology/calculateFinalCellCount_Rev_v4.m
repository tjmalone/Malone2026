%% calculateFinalCellCount_Manual
% for 1.5 mm
% for reelin/calbindin/tau analysis


%% Load previous data

clear; close all; clc

% load global mouse info
load('Z:\labMembers\KC\AD_Project\Histology\miceInfo.mat','mice','miceKey')

% load sub mouse info
load('Z:\labMembers\KC\AD_Project\Histology\tau_RC\ADMouseOrderShort.mat','mouseOrder')
nMice = size(mouseOrder,1);

% initialize data matrix
dataOld = struct();
% dataOld.geno = zeros(2*nMice,1);
% dataOld.sex = cell(2*nMice,1);
% dataOld.age = zeros(2*nMice,1);

% find mouse info
for mm = 1:nMice
    % find current mouse match
    mouseID = mouseOrder{mm,1};
    mouseIdx = find(strcmp(mouseID,mice(:,1)));

    % store mouse info
    idx = mm*2-1:mm*2;
    [dataOld(idx).ID] = deal(mouseID);
    [dataOld(idx).geno] = deal(mice{mouseIdx,3});
    [dataOld(idx).sex] = deal(mice{mouseIdx,4});
    [dataOld(idx).age] = deal(mice{mouseIdx,6});
end

% load cell counts
baseFold1 = 'Z:\labMembers\KC\AD_Project\Histology\tau_RC\FinalVersion\';
load('Z:\labMembers\KC\AD_Project\Histology\tau_RC\FinalVersion\calculateFinalCellCount.mat','cellCountFinal');

% store cell counts
tmp = num2cell(cellCountFinal(:,1)); [dataOld.reel]     = tmp{:};
tmp = num2cell(cellCountFinal(:,2)); [dataOld.cal]      = tmp{:};
tmp = num2cell(cellCountFinal(:,3)); [dataOld.tau]      = tmp{:};
tmp = num2cell(cellCountFinal(:,4)); [dataOld.reel_tau] = tmp{:};
tmp = num2cell(cellCountFinal(:,5)); [dataOld.cal_tau]  = tmp{:};


%% Calculate new overlaps

p1 = 'Z:\labMembers\TM\AD_Project\Revision\Histology\revisionHistology';
cd(p1)

% find all folders
foldNew = findSubF('reel_cal_AT8',1,[],0);
nMice = length(foldNew);

% load mouse info
load('mouseInfo.mat','mouseInfo')
if length(mouseInfo.ID)~=nMice
    error('Mouse number mismatch')
end

% initialize data matrix
dataNew = struct();

sfxCodes = {'A','B'};
overlayCodes = {'overlay_reel-cal','overlay_reel-tau','overlay_cal-tau','overlay_full'};
nOC = length(overlayCodes);
tmpDir = tempname;

for mm = 1:nMice
    cd(foldNew{mm})

    % find image files
    mec = dir('MEC_*.tif');
    if length(mec)~=2
        error('Incorrect image number')
    end

    for ff = 1:2
        curIdx = 2*(mm-1)+ff;

        % store mouse info
        dataNew(curIdx).ID = mouseInfo.ID{mm};
        dataNew(curIdx).geno = 2;
        dataNew(curIdx).sex = mouseInfo.sex{mm};
        dataNew(curIdx).age = mouseInfo.age{mm};

        % find suffix numbers
        mecID = str2double(mec(ff).name(end-4));

        if isfile(['selectSkip_' sfxCodes{ff} '.txt'])
            % set overlaps for invalid FOV
            dataNew(curIdx).reel = NaN;
            dataNew(curIdx).cal = NaN;
            dataNew(curIdx).tau = NaN;
            dataNew(curIdx).reel_tau = NaN;
            dataNew(curIdx).cal_tau = NaN;
        elseif isfile(['selectNA_'  sfxCodes{ff} '.txt'])
            % set overlaps when no tau was present
            dataNew(curIdx).reel = 1;
            dataNew(curIdx).cal = 1;
            dataNew(curIdx).tau = 0;
            dataNew(curIdx).reel_tau = 0;
            dataNew(curIdx).cal_tau = 0;
        else
            % move to overlay folder
            cd(['OverlayMEC' num2str(mecID)])

            % load total cell numbers
            load('numRois.mat','numRois')

            % process fiels
            nCells = zeros(1,nOC);
            for cc = 1:nOC
                if isfile([overlayCodes{cc} '.roi'])
                    % all .roi files are 1 roi
                    nCells(cc) = 1;
                elseif isfile([overlayCodes{cc} '.zip'])
                    % count number of rois in zip folder
                    nCells(cc) = numel(unzip([overlayCodes{cc} '.zip'],tmpDir));
                    rmdir(tmpDir, 's');
                else
                    % absent files is assumed 0 cells
                    nCells(cc) = 0;
                end
            end

            % set overlaps (be careful as these indices are hardcoded)
            dataNew(curIdx).reel = numRois(1)-nCells(4);
            dataNew(curIdx).cal = numRois(2)-nCells(4);
            dataNew(curIdx).tau = numRois(3);
            dataNew(curIdx).reel_tau = nCells(2)-nCells(4);
            dataNew(curIdx).cal_tau = nCells(3)-nCells(4);

            cd(foldNew{mm})
        end
    end
end

cd(p1)


%% Combine data

% combine all data
dataAll = cat(2,dataOld,dataNew)';
nMice = length(dataAll);

% calculate percentages
for mm = 1:nMice
    dataAll(mm).pToR = dataAll(mm).reel_tau./dataAll(mm).reel*100;
    dataAll(mm).pToC = dataAll(mm).cal_tau./dataAll(mm).cal*100;
    dataAll(mm).pRoT = dataAll(mm).reel_tau./dataAll(mm).tau*100;
    dataAll(mm).pCoT = dataAll(mm).cal_tau./dataAll(mm).tau*100;
end

% calcualte combined mouse info
sexAll = cellfun(@(x) x(1)=='F',{dataAll.sex})'+1;
ageAll = [dataAll.age]';

% combine all data
flds = {'pToR','pToC'};
combData = [];
for ii = 1:2
    combData = [combData; [dataAll.(flds{ii})]'];
end

% find fit of combined data
fitData = combData(~isnan(combData));
tmpAge = [ageAll;ageAll];
fitAge = tmpAge(~isnan(combData));
[mdl,gof] = fitGompertz(fitAge,fitData);

% define fit input parameters
Y0 = mdl.Y0;
YMax = mdl.Ymax;
StartPoint = [Y0 YMax 1];
Lower = [Y0 YMax 0.01];
Upper = [Y0 YMax 10];


%% Plot overall

figure; hold on
ttls = 'cells that are tau+';
flds = {'pToR','pToC'};
colors = {[1 0 0],[0 0 1]};
labs = {'reelin','calbindin'};
ylims = [0 30];
xMax = 10;
xFit = linspace(0,xMax,200);

% intialize final statistic
outStats = {};

% plot scatter
h = zeros(1,2);
x = [];
y = [];
g = [];
for ii = 1:2
    % current data
    curData = [dataAll.(flds{ii})]';
    plotAge = ageAll(~isnan(curData));
    plotData = curData(~isnan(curData));
    h(ii) = scatter(plotAge,plotData,30,colors{ii},'filled');

    % store data for fit
    x = [x;plotAge];
    y = [y;plotData];
    g = [g;ii*ones(size(plotAge))];
end

% [Y0 Ymax K], 1 = shared, 0 = group-specific
sharedParam = [1 1 0];

[p,yFits,b1,b2,outCur] = fitGompertzByGroup(x,y,g,xFit,StartPoint,Lower,Upper,sharedParam,' FOVs');
% [p,yFits,b1,b2] = fitGompertzByGroup(x,y,g,xFit,[],[],[],sharedParam);
for ii = 1:2
    plot(xFit,yFits{ii},'Color',colors{ii})
end

% plot mutual plateau
yline(b1(2),'Color',[0.5 0.5 0.5])

% set plot limits and labels
legend(h,labs)
ylim(ylims)
xlabel('Age (months)');
xlim([0 xMax])
title([ttls ', p = ' num2str(p,2)])

% add final stats
testCat = 'Global fit';
outStats(end+1,:) = [{testCat},outCur(:)'];


%% Plot by cell type

% initialize figure
figure; tiledlayout(2,2)
ttls = {'reelin+ that are tau+','calbindin+ that are tau+'};
flds = {'pToR','pToC'};
statsLab = {'reelin','calbindin'};
colors = {[175 0 0]/255;[255 102 102]/255};
labs = {'Male','Female'};
ylims = [0 30];
xMax = 10;
xFit = linspace(0,xMax,200);

% plot scatter
sseSex = zeros(1,2);
dfSex = zeros(1,2);
for ii = 1:2
    nexttile(ii); hold on

    % current data
    curData = [dataAll.(flds{ii})]';

    % plot scatter
    h = zeros(1,2);
    x = [];
    y = [];
    g = [];
    for ss = 1:2
        plotData = curData(sexAll==ss & ~isnan(curData));
        plotAge = ageAll(sexAll==ss & ~isnan(curData));

        h(ss) = scatter(plotAge,plotData,30,colors{ss},'filled');

        % store data for fit
        x = [x;plotAge];
        y = [y;plotData];
        g = [g;ss*ones(size(plotAge))];
    end

    % [Y0 Ymax K], 1 = shared, 0 = group-specific
    sharedParam = [1 1 0];

    [p,yFits,b1,b2,outCur] = fitGompertzByGroup(x,y,g,xFit,StartPoint,Lower,Upper,sharedParam,' FOVs');
    % [p,yFits,b1,b2] = fitGompertzByGroup(x,y,g,xFit,[],[],[],sharedParam);
    for ss = 1:2
        plot(xFit,yFits{ss},'Color',colors{ss})
    end

    % plot mutual plateau
    yline(b1(2),'Color',[0.5 0.5 0.5])

    % set plot limits and labels
    legend(h,labs)
    ylim(ylims)
    xlabel('Age (months)');
    xlim([0 xMax])
    title([ttls{ii} ', p = ' num2str(p,2)])

    % add final stats
    testCat = [statsLab{ii} ' fit'];
    outStats(end+1,:) = [{testCat},outCur(:)'];
end


%% Plot by sex

ttls = {'cells that are tau+','cells that are tau+'};
flds = {'pToR','pToC'};
colors4 = {[0 0 175]/255,[175 0 0]/255;[102 105 255]/255,[255 102 102]/255};
labs = {'reelin','calbindin'};
sexTypes = {'Male','Female'};

% plot scatter
sseType = zeros(1,2);
dfType = zeros(1,2);
for ss = 1:2
    nexttile(ss+2); hold on

    % plot scatter
    h = zeros(1,2);
    x = [];
    y = [];
    g = [];
    for ii = 1:2
        % current data
        curData = [dataAll.(flds{ii})]';
        plotAge = ageAll(sexAll==ss & ~isnan(curData));
        plotData = curData(sexAll==ss & ~isnan(curData));
        h(ii) = scatter(plotAge,plotData,30,colors4{ss,ii},'filled');

        % store data for fit
        x = [x;plotAge];
        y = [y;plotData];
        g = [g;ii*ones(size(plotAge))];
    end

    % [Y0 Ymax K], 1 = shared, 0 = group-specific
    sharedParam = [1 1 0];

    [p,yFits,b1,b2,outCur] = fitGompertzByGroup(x,y,g,xFit,StartPoint,Lower,Upper,sharedParam,' FOVs');
    % [p,yFits,b1,b2] = fitGompertzByGroup(x,y,g,xFit,[],[],[],sharedParam);
    for ii = 1:2
        plot(xFit,yFits{ii},'Color',colors4{ss,ii})
    end

    % plot mutual plateau
    yline(b1(2),'Color',[0.5 0.5 0.5])

    % set plot limits and labels
    legend(h,labs)
    ylim(ylims)
    xlabel('Age (months)');
    xlim([0 xMax])
    title([ttls{ii} ', p = ' num2str(p,2)])

    % add final stats
    testCat = [sexTypes{ss} ' fit'];
    outStats(end+1,:) = [{testCat},outCur(:)'];
end


%% Plot combined bar graph

ageLim = 6.5;

ttls = {'cells that are tau+','tau+ that are reelin/calbindin+'};
fld1 = {'pToR','pRoT'};
fld2 = {'pToC','pCoT'};
labs = {'reelin','calbindin'};
colors = {[1 0 0],[0 0 1]};

figure; tiledlayout(1,2)
mapAll = cell(2,2);
for ii = 1:2
    nexttile(ii); hold on
    stat1 = [dataAll(ageAll>=ageLim).(fld1{ii})];
    stat2 = [dataAll(ageAll>=ageLim).(fld2{ii})];

    [~,p] = ttest(stat1,stat2);
    barGroup(labs,{stat1,stat2},'violin',colors,[1 2],'pair',[],0)

    title(ttls{ii})

    % store data for heatmap
    mapAll{ii,1} = stat1;
    mapAll{ii,2} = stat2;

    % add final statistics
    nUnits = ' FOVs';
    testName = 'two-tailed paired Students t-test';
    testPair = 1;
    testMC = 0;
    testLimitP = 0;
    testCat = [fld1{ii} ' vs. ' fld2{ii}];
    outStats(end+1,:) = [testCat ttestEffectSize(stat1',stat2',testName,nUnits,testPair,testMC,testLimitP)];
end


%% Plot heat map

for ii = 1:2
    morphData = mapAll(ii,:)';

    % define data
    curMap = cell(3,3,2);
    for kk = 1:numel(curMap)
        curMap{kk} = nan;
    end
    curMap(1,2:3,:) = cat(2,{0;0},morphData);

    % make heat map
    h = heatmapGridSplit(curMap,{'Sex','All','Female','Male'},...
        {'Morphology','All','ste','pyr'},-1,{1},{2:3},{'eastoutside'});
    sgtitle([ttls{ii}],'FontSize',18)
end


%% Plot bar graph by sex

ageLim = 6.5;

ttls = {'cells that are tau+','tau+ that are reelin/calbindin+'};
fld1 = {'pToR','pRoT'};
fld2 = {'pToC','pCoT'};
% labs = {'reelin','calbindin'};
xLabs = {'Male','Female'};

colors = {[1 0 0],[0 0 1]};

figure; tiledlayout(1,2)
for ii = 1:2
    nexttile(ii); hold on

    % collate plot data
    plotData = cell(2,2);
    for ss = 1:2
        plotData{ss,1} = [dataAll(sexAll==ss & ageAll>=ageLim).(fld1{ii})]';
        plotData{ss,2} = [dataAll(sexAll==ss & ageAll>=ageLim).(fld2{ii})]';
    end

    sigPair = [1,2;1,3;2,4;3,4];
    [~,p] = barGroup(xLabs,plotData,'violin',colors,sigPair,'pair',[],0);

    title(ttls{ii})

    % add final statistics
    nUnits = ' FOVs';
    testName = 'two-tailed paired Students t-test';
    testPair = 1;
    testMC = 0;
    testLimitP = 0;
    testCat = [fld1{ii} ' vs. ' fld2{ii} ': by sex'];
    statData = plotData;
    outStats(end+1,:) = [testCat ttestEffectSize(statData(:,1),statData(:,2),testName,nUnits,testPair,testMC,testLimitP)];

    % add final statistics
    nUnits = ' FOVs';
    testName = 'two-tailed unpaired Students t-test';
    testPair = 0;
    testMC = 0;
    testLimitP = 0;
    testCat = [fld1{ii} ' vs. ' fld2{ii} ': by sex (cell type)'];
    statData = plotData';
    outStats(end+1,:) = [testCat ttestEffectSize(statData(:,1),statData(:,2),testName,nUnits,testPair,testMC,testLimitP)];
end


%% Plot normalized bar graph by sex

ageLim = 6.5;

ttls = {'cells that are tau+','tau+ that are reelin/calbindin+'};
fld1 = {'pToR','pRoT'};
fld2 = {'pToC','pCoT'};
% labs = {'reelin','calbindin'};
xLabs = {'Male','Female'};

colors = {[1 0 0],[0 0 1]};

figure; tiledlayout(1,2)
for ii = 1:2
    nexttile(ii); hold on

    % collate plot data
    plotData = cell(2,2);
    for ss = 1:2
        plotData{ss,1} = [dataAll(sexAll==ss & ageAll>=ageLim).(fld1{ii})]';
        plotData{ss,2} = [dataAll(sexAll==ss & ageAll>=ageLim).(fld2{ii})]';

        % normalize data
        normData = plotData{ss,1};
        plotData{ss,1} = plotData{ss,1}./normData;
        plotData{ss,2} = plotData{ss,2}./normData;

        % remove nans
        plotData{ss,1} = plotData{ss,1}(normData~=0);
        plotData{ss,2} = plotData{ss,2}(normData~=0);
    end

    sigPair = [1,3;2,4;3,4];
    barGroup(xLabs,plotData,'violin',colors,sigPair,[],[],0)

    title(ttls{ii})
end
