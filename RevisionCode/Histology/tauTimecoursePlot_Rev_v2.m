%% Load previous data

clear; close all; clc

% set base folder
p1 = 'Z:\labMembers\TM\AD_Project\Revision\Histology\Leica';
cd(p1)

% load new data
load('tauTimecourse.mat')

% load old data
load('Z:\labMembers\KC\AD_Project\Histology\AD_timecourse_Leica\FinalVersion\mouseInfoComb.mat',...
    'combMiceName','combGeno','combSex','combTrueAge','combTimecourseIdx');
nMice = length(combMiceName);

% initialize data matrix
dataOld = struct();
dataOld(nMice,1) = struct();

% store cell counts
[dataOld.ID] = combMiceName{:};
tmp = num2cell(combGeno); [dataOld.geno] = tmp{:};
[dataOld.sex] = combSex{:};
tmp = num2cell(combTrueAge); [dataOld.age] = tmp{:};
tmp = num2cell(combTimecourseIdx); [dataOld.batch] = tmp{:};
tmp = num2cell(nan(nMice,1)); [dataOld.pixelInt] = tmp{:};

% load raw intenisty
load('Z:\labMembers\KC\AD_Project\Histology\AD_timecourse_Leica\FinalVersion\raw\tauDataGrouped.mat', 'tauInt')
raw = tauInt(:,2);
tmp = num2cell(raw); [dataOld.raw] = tmp{:};

% load positive area
load('Z:\labMembers\KC\AD_Project\Histology\AD_timecourse_Leica\FinalVersion\perArea\tauDataGrouped.mat', 'tauInt')
area = tauInt(:,2);
tmp = num2cell(area); [dataOld.area] = tmp{:};

% load positive area intensity
load('Z:\labMembers\KC\AD_Project\Histology\AD_timecourse_Leica\FinalVersion\backgroundSub\tauDataGrouped.mat', 'tauInt')
posInt = tauInt(:,2);
tmp = num2cell(posInt); [dataOld.posInt] = tmp{:};


%% Combine and plot data

% combine all data
dataAll = cat(1,dataOld,dataNew');
nMice = length(dataAll);

% calcualte combined mouse info
sexAll = cellfun(@(x) x(1)=='F',{dataAll.sex})'+1;
ageAll = [dataAll.age]';
genoAll = [dataAll.geno]';

% define labels
typesGeno = {'WT','PS19'};
colors2 = {[0 0 1],[1 0 0]};
nGeno = length(typesGeno);
typesSex = {'Male','Female'};
nSex = length(typesSex);
colors4 = {[0 0 175]/255,[175 0 0]/255;[102 105 255]/255,[255 102 102]/255};

% define age groups
ageThresh = 6.51;
ageM = [1,-1];
typesAge = {'young','old'};

% initialize figure
nPlot = 4;
figure; tiledlayout(2,nPlot)
flds = {'area','posInt'};
yLims = {[0 100],[0 0.25]};
yLabs = {'% Area Tau+','Tau intensity (A.U.)'};
xMax = 10;
xFit = linspace(0,xMax,200);

% initialize full stats
outStats = {};

for ff = 1:length(flds)
    %% Plot combined sex scatter

    % current data type
    curPlot =[dataAll.(flds{ff})]';

    % initialize panel
    nexttile(nPlot*(ff-1)+1); hold on

    % plot scatter
    for gg = 1:nGeno
        scatter(ageAll(genoAll==gg),curPlot(genoAll==gg),30,colors2{gg},'filled')

        if gg==2
            % perform fit of PS19 data
            x = ageAll(genoAll==gg & ~isnan(curPlot));
            y = curPlot(genoAll==gg & ~isnan(curPlot));
            mdl = fitGompertz(x,y);

            % plot fit
            yFit = mdl(xFit);
            plot(xFit,yFit,'Color',colors2{gg})
            yline(mdl.Ymax)

            % define final statistics
            testCat = 'Group fit';
            testType = 'Fit of Gompertz exponential';
            names = string(coeffnames(mdl));
            vals = string(coeffvalues(mdl))';
            terms = compose("%s=%.2g", names, vals);
            outCoeff = strjoin(terms, ", ");
            Ns = ['n=' num2str(length(y)) ' mice'];
            outStats(end+1,1:8) = {testCat,'',testType,outCoeff,Ns,'','',''};
        end
    end

    % set plot limits and labels
    ylim(yLims{ff})
    ylabel(yLabs{ff});
    xlabel('Age (months)');
    xlim([0 xMax])


    %% Plot age based bar graph

    % define bar data
    barData = cell(2,nGeno);
    for gg = 1:nGeno
        for aa = 1:2
            barData{aa,gg} = curPlot(genoAll==gg & ageM(aa)*ageAll<ageM(aa)*ageThresh);

            if all(isnan(barData{aa,gg}))
                barData{aa,gg} = 0;
            end
        end
    end

    % plot bar graph
    nexttile(nPlot*(ff-1)+2); hold on
    sigPair = [1,3;2,4;3,4];
    [~,p] = barGroup(typesAge,barData,'violin',colors2,sigPair,[],[],0);

    % set labels
    ylim(yLims{ff})
    ylabel(yLabs{ff});
    xlabel ('Group')

    % perform full statistics
    nUnits = 'mice';
    testName = 'two-tailed unpaired Students t-test';
    testPair = 0;
    testMC = 0;
    testLimitP = 0;
    statData = [barData;[barData(1,2) barData(2,2)]];
    testCat = ['Group bar: ' flds{ff}];
    outStats(end+1,:) = [testCat ttestEffectSize(statData(:,1),statData(:,2),testName,nUnits,testPair,testMC,testLimitP)];


    %% Plot sex based bar graph

    % plot scatter
    dataSex = cell(nSex,2);
    for ss = 1:nSex
        for aa = 1:2
            dataSex{ss,aa} = curPlot(genoAll==2 & sexAll==ss & ageM(aa)*ageAll<ageM(aa)*ageThresh);
        end
    end

    % plot bar graph
    nexttile(nPlot*(ff-1)+3); hold on
    sigPair = [1,3;2,4;3,4];
    [~,p] = barGroup(typesSex,dataSex,'violin',colors4,sigPair,[],[],0);

    % set labels
    ylim(yLims{ff})
    ylabel(yLabs{ff});
    xlabel ('Group')

    % perform full statistics
    nUnits = 'mice';
    testName = 'two-tailed unpaired Students t-test';
    testPair = 0;
    testMC = 0;
    testLimitP = 0;
    testCat = ['Sex bar:' flds{ff}];
    statData = [dataSex;[dataSex(1,2) dataSex(2,2)]];
    outStats(end+1,:) = [testCat ttestEffectSize(statData(:,1),statData(:,2),testName,nUnits,testPair,testMC,testLimitP)];


    %% Plot timecourse by sex

    % plot scatter
    nexttile(nPlot*(ff-1)+4); hold on
    sseSex = zeros(1,nSex);
    dfSex = zeros(1,nSex);
    x = [];
    y = [];
    g = [];
    for ss = 1:nSex
        useIdx = genoAll==2 & sexAll==ss & ~isnan(curPlot);
        useAge = ageAll(useIdx);
        useData = curPlot(useIdx);
        scatter(useAge,useData,30,colors4{ss,2},'filled')

        % store data for fit
        x = [x;useAge];
        y = [y;useData];
        g = [g;ss*ones(size(useAge))];
    end

    % [Y0 Ymax K], 1 = shared, 0 = group-specific
    sharedParam = [1 1 0];

    % [p,yFits] = fitGompertzByGroup(x,y,g,xFit,StartPoint,Lower,Upper);
    [p,yFits,b1,b2,outCur] = fitGompertzByGroup(x,y,g,xFit,[],[],[],sharedParam,' mice');
    for ss = 1:2
        plot(xFit,yFits{ss},'Color',colors4{ss,2})
    end

    % plot mutual plateau
    yline(b1(2),'Color',[0.5 0.5 0.5])

    % set plot limits and labels
    ylim(yLims{ff})
    ylabel(yLabs{ff});
    xlabel('Age (months)');
    xlim([0 xMax])
    title(['p = ' num2str(p,2)])

    % define final statistics
    testCat = 'Sex comparison fit';
    outStats(end+1,:) = [{testCat},outCur(:)'];
end

