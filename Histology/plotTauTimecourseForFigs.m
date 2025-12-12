%% quantTauLayersKC_AllTimecourses
% For Timecourses A B C
%%%%%
% Quantifies tau for the timecourse
% Updates from Version 2
% - new threshold calculated by finding the xth percentile
% - outputs results for 3 options
        % - Raw
        % - Background subtraction tau intensity
        % - Percent area
%%%%%


%% Define inputs

clear; close all; clc

% change to analysis directory
baseFold = 'Z:\labMembers\KC\AD_Project\Histology\AD_timecourse_Leica\FinalVersion';
cd(baseFold)

% analysis types
analysisType = {'raw','backgroundSub','perArea'};

% dfine axes labels
yAxis = {'Tau intensity (A.U.)','Tau intensity (A.U.)','% Area Tau+'};

% define colors
colors2 = {[0 0 1],[1 0 0]};

% load categories
load('cats.mat')
nAge = length(catAge);
nGeno = length(catGeno);
useLayer = 2;

% define age categories
ageCatBins = {1:4,5:7};
ageCatVals = {[2 3 5 6],[7 8 9]};
xGroup = {'Young', 'Old'};

% define sexes
sexes = {'Female','Male'};
nSex = length(sexes);

% load mouse data
load('mouseInfoComb.mat')

% define analyses to run
rUse = 2:3;


%% Plot time course

% initialize stats parameters and array
outStats = {};
nUnits = 'mice';
testName = 'two-tailed unpaired Students t-test';
testPair = 0;
testMC = 0;
testLimitP = 0;

% loop through analysisType (raw, backsub, perArea)
for r = rUse
    cd(analysisType{r})

    load('tauDataGrouped.mat', 'tauInt', 'tauData')

    % plot by layer
    for kk = useLayer
        %% Colect layer data
        dataLayer = cellfun(@(x) x(:,kk),tauData,'UniformOutput',false);

        % calculate significance
        [curAnova,curMC,pLabel] = anovaO2W_BH(dataLayer,0);
        pMC = curMC';
        pAnova = curAnova(1)*ones(size(pMC));


        %% Plot line

        figure; hold on

        % plot averages
        plotErrorSigCell(catAge,dataLayer,catGeno,pAnova,pMC,colors2)

        % plot scatter
        for gg = 1:nGeno
            scatter(combTrueAge(combGeno==gg),tauInt(combGeno==gg,kk),30,colors2{gg},'filled')
        end

        % set plot limits and labels
        ylim([0 max(tauInt(:,kk))])
        title(['Tau Intensity (Layer' num2str(kk) ')']);
        ylabel(yAxis(r));
        xlabel('Age (months)');

        % save figure
        savefig(['tauIntensityTimecourse_L' num2str(kk)]);


        %% Plot bar

        % define bar data
        barData = cell(length(ageCatBins),nGeno);
        for gg = 1:nGeno
            for aa = 1:length(ageCatBins)
                barData{aa,gg} = cat(1,dataLayer{gg,ageCatBins{aa}});
            end
        end

        % plot bar graph
        figure;
        sigPair = [1,3;2,4;3,4];
        barGroup(xGroup,barData,'violin',colors2,sigPair,[],[],0)

        % set labels
        ylabel(yAxis(r))
        xlabel ('Group')
        title([analysisType{r} ': Layer ' num2str(kk) ' tau intensity'])

        % save figure
        savefig(['ageGroupBar-byGeno_L' num2str(kk)])

        % perform full statistics
        statData = [barData;[barData(1,2) barData(2,2)]]';
        testCat = ['ageGroupBar-byGeno_L' num2str(kk)];
        outStats(end+1,:) = [testCat ttestEffectSize(statData(1,:),statData(2,:),testName,nUnits,testPair,testMC,testLimitP)];


        %% Split AD mice by sex

        % plot scatter
        dataSex = cell(nSex,nGeno);
        for ss = 1:nSex
            for aa = 1:length(ageCatVals)
                dataSex{ss,aa} = tauInt(combGeno==2 & strcmp(combSex,sexes{ss})...
                    & ismember(combAge,ageCatVals{aa}),kk);
            end
        end

        % plot bar graph
        figure;
        sigPair = [1,3;2,4;3,4];
        barGroup(sexes,dataSex,'violin',colors2,sigPair,[],[],0)

        % set labels
        ylabel(yAxis(r))
        xlabel ('Group')
        title([analysisType{r} ': Layer ' num2str(kk) ' tau intensity'])

        % save figure
        savefig(['ageGroupBar-bySex_L' num2str(kk)])

        % perform full statistics
        if r==2
            statData = [dataSex(2,2); dataSex(1,2)];
        elseif r==3
            statData = [flipud(dataSex);[dataSex(2,2) dataSex(1,2)]]';
        end
        testCat = ['ageGroupBar-bySex_L' num2str(kk)];
        outStats(end+1,:) = [testCat ttestEffectSize(statData(1,:),statData(2,:),testName,nUnits,testPair,testMC,testLimitP)];

    end

    
    cd(baseFold)
end

