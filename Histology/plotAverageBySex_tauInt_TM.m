%%  plotAverageBySex_tauInt2
% plots tau intensity for comparing Males and Females

% Corrects "youngVsOldBarGraph" and "layer2YvsO"
% Makes Bar graph called "male vs female"

% For Timecourses ABC
%%%%%
% plots by sex
% - outputs results for 3 options
        % - Raw
        % - Background subtraction tau intensity
        % - Percent area
%%%%%

%% Part 1: Set values

clear; close all; clc

%%%%%%%%%%%%%%%%%%%Source of Data
foldVers = '99PercentileV1';
saveFold = 'bySex';
%%%%%%%%%%%%%%%%%%%%
nLy = 3;

foldTimecourse = 'Z:\labMembers\KC\AD_Project\Histology\AD_timecourse_Leica\FinalVersion';
cd(foldTimecourse)

% identify all timecourses
f = 'Timecourse*';
d = dir(f);
numTimecourse = size(d,1);
timecourseNames = {};
baselineInt = {};
analysisType = {'raw','backgroundSub','perArea'};
yAxis = {'Tau intensity per area (A.U.)', '% Area Tau+'};

% Load and set details about mice
%Make folder for each analysis type
foldVersion = pwd;
load('mouseInfoComb.mat')
load('mouseInfo.mat');
load('cats.mat');

% Ages
catAge = [2, 3, 5, 6, 7, 8, 9];
nAge = length(catAge);

% genotypes
% WT=1, AD=2
catGeno = {'WT','AD'};
nGeno = length(catGeno);

% Define 4 groups 
xGroupNames = [{'wtMale'} {'adMale'}  {'wtFemale'} {'adFemale'}];
% groupNum -> (1, 2, 3, 4)
yGroupNames = [{'young','old'}];
% groupNum -> (1, 2)
nXGroups = length(xGroupNames);
nYGroups = length(yGroupNames);
xGroups = [1:length(xGroupNames)];

%% Part 2. sort into 8 groups
% order of splitting:
    % Sex
    % Geno
    % age
groupNum = []; 
for i = 1:length(combSex)
     % split by sex and geno
    if combSex{i,1}(1,1) == 'M'
        if combGeno(i,1) == 1
            groupNum(i,1) = 1;
        else
            groupNum(i,1) = 2;
        end
    else
        if combGeno(i,1) == 1
            groupNum(i,1) = 3;
        else
            groupNum(i,1) = 4;
        end
    end

    % split by age
    if combAge(i,1) < 7
        groupNum(i,2) = 1;
    else
        groupNum(i,2) = 2;
    end
end

% %2. sort into 4 groups (young vs old)
% layer2YvsO

%% Part 3: Loop through each analysis type and make graphs

% loop through analysisType (raw, backsub, perArea)
for r = 3:length(analysisType)

% testing
% r = 1;
     if r == 3
        yAxisType = 2;
     else
        yAxisType = 1;
     end
    % Load files
    cd(foldVersion)
    cd(analysisType{r})
    load('tauDataGrouped.mat')

    mkdir(saveFold)
    cd(saveFold)
    
    % Organize data for based on sex, geno, and age (layer 2 only)
    tauDataBySexAge = cell(nYGroups,nXGroups);
    ageBySex = cell(nYGroups,nXGroups);
    tauMean = zeros(nXGroups,nAge);
    tauSEM = zeros(nXGroups,nAge);

    for ii = 1:size(tauInt,1)
        tauDataBySexAge{groupNum(ii, 2),groupNum(ii,1)}(end+1,1) = tauInt(ii,2);
        ageBySex{groupNum(ii, 2),groupNum(ii,1)}(end+1,1) = ageIndex(ii);
    end
    save('tauDataBySex_imaging.mat','tauDataBySexAge','ageBySex')
    
    % Graph 1: scatter plot by sex
    line = 1;
    ll = 1;
    figure; hold on
    colors = {[0 0 1],[1 0 0],[.7 .7 1],[1 .7 .7]};
    if line == 2
        for gg = 1:nXGroups
            errorbar(catAge,tauMean(gg,:,ll),tauSEM(gg,:,ll),...
                'Color',colors{1,gg},...
                'DisplayName',[xGroupNames{gg} ': Layer ' catLayer{ll}])
        end
        nameEnding = 'line';
    else
        nameEnding = 'noLine';
    end

    for i = 1:nXGroups
        for j = 1:nYGroups
            scatter(catAge(ageBySex{j,i}), tauDataBySexAge{j,i}, 30, colors{1,i});
        end
    end

    titleString = ['Layer 2 ' analysisType{r}];
    title(titleString);
    ylabel(yAxis(yAxisType));
    legend(xGroupNames,'Location','Northwest');
    xlabel 'Age (months)';
    hold off
    savefig([titleString '_' nameEnding '_imaging']);

    %% Graph 2: bar by sex for old mice

    colors = {[1 0 0],[1 .7 .7]};

    oldMice = {};
    oldMice = tauDataBySexAge(:,[2 4]);
    oldMiceX = xGroupNames([2 4]);
    figure;
    [~,pMC] = barGroup({'young','old'},oldMice,colors,[1 3; 2 4; 1 2; 3 4]);
    ylabel(yAxis(yAxisType))
    xlabel ('Group')
    title([analysisType{r} ': Layer 2 tau intensity (p = ' num2str(pMC(2,4),2) ')'])
    legend(oldMiceX,'Location','northwest')
    % savefig(['maleVsFemaleBarGraph_imaging'])

    % % Graph 3: young vs old mice bar
    % oldMice = {};
    % oldMice = tauDataBySex(2,:);
    % % sigPair = [1,4;2,4;3,4];
    % sigPair = nchoosek(xGroups,2);
    % figure;
    % barGroup(xGroupNames,oldMice,[],sigPair)
    % ylabel(yAxis(yAxisType))
    % xlabel ('Group')
    % title([analysisType{r} ': Layer 2 tau intensity'])
    % savefig(['maleVsFemaleBarGraph'])
end

% close all;


            
