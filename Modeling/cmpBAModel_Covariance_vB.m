%% cmpBAModel_Covariance
%  Plot the covariance between activity variables

clear; close all; clc

% set folder
p1 = 'D:\AD_Project\Behavior';
cd(p1)

% load genotypes and sexes
load('groupIDs.mat')
groupNames = [groupIDs, 'All'];
% load data
load('data/model_vB/corrData_vB.mat','corrData','corrLabels')

% define use days
useDays = 2:11;

% set regression parameters
paramsUse = 1:25;
paramsUse = [1 3 15 9 17 14 6 4 5 19 24];

nParams = length(paramsUse);
meanStartIdx = 2;

% define factor names
varNames = corrLabels(paramsUse);
% varNames = {'Success Rate','Absolute Speed Score','Spatial Selectivity',...
%     'Activity Level','Raw Speed Score','Within-Day Consistency',...
%     'Cross-Day Consistency','Inter-Mouse Consistency',...
%     'Decoding Cue Anchoring','Spatial Activity Cue Anchoring'};

% define colors
colors2 = {[0 0 1],[1 0 0]};


%% Calculate and plot correlation

dataVectorTL = cell(1,3);
dataCV = zeros(nParams,3);
dataCorrB = zeros(nParams-1,3);
dataCorrA = zeros(nParams-1,3);

for gg = 1:3
    % load current data
    if gg<3
        useGroup = groups{gg};
    else
        useGroup = 1:size(corrData,1);
    end

    dataUse = reshape(corrData(useGroup,useDays,paramsUse),[],nParams);

    %% Plot covariance matrix

    corrAbs = (corr(dataUse,'rows','pairwise'));
    corrAbs(eye(nParams)==1) = NaN;

    % get lower triangle
    tl = tril(true(nParams-1));
    dataAct = abs(corrAbs(meanStartIdx:end,meanStartIdx:end));
    dataVectorTL{gg} = dataAct(tl);

    % calculate mean covariance
    meanCV = mean(dataVectorTL{gg},'omitnan');

    % corrAbs(corrAbs<0.55) = 0;

    % make heat map with transparency
    figure; hold on
    imAlpha = ones(nParams);
    imAlpha(isnan(corrAbs))=0;
    imagesc(corrAbs,'AlphaData',imAlpha);

    % set background color
    set(gca,'color',0*[1 1 1]);
    xlim([0.5 nParams+0.5])
    ylim([0.5 nParams+0.5])

    % set color map
    curCMap = customcolormap([0 0.25 0.5 0.75 1],[0 0 1; 0.5 0.5 1; 1 1 1; 1 0.5 0.5; 1 0 0]);
    colormap(curCMap)
    colorbar
    clim([-1 1])

    xticks(1:nParams)
    yticks(1:nParams)
    xticklabels(varNames)
    yticklabels(varNames)
    axis('square')
    title(['Variable Covariance: ' groupNames{gg} ', mean = ' num2str(meanCV,2)])

    % store standard deviation
    dataCV(:,gg) = std(dataUse,1,'omitnan')./abs(mean(dataUse,1,'omitnan'));

    % store correlation to behavior
    dataCorrB(:,gg) = abs(corrAbs(2:end,1));

    % store correlation to activity
    dataCorrA(:,gg) = mean(abs(corrAbs(2:end,2:end)),1,'omitnan');


end

[~,pCorr] = ttest(dataVectorTL{1},dataVectorTL{2});
fprintf('Absolute Correlation Difference: p = %.3f\n',pCorr)

% %% Plot standard deviation
% 
% figure; hold on
% barGroup(varNames,num2cell(dataCV(:,1:2)),'bar',colors2)
% title('Variable variation')
% 
% 
% %% Plot behavior correlation
% 
% figure; hold on
% barGroup(varNames(2:end),num2cell(dataCorrB(:,1:2)),'bar',colors2)
% title('Absolute correlation with behavior')
% 
% 
% %% Plot activity correlation
% 
% figure; hold on
% barGroup(varNames(2:end),num2cell(dataCorrA(:,1:2)),'bar',colors2)
% title('Absolute correlation with activity variables')
% 
