function [hh,hhSub,pAnova,pMC] = plotErrorMulti(X,data,legLabels,colors,styles,pShow,legLabelsShort,mcON,useIdx)
%% plotErrorMulti
% Plots a line graph with error bars for an arbitrary number of data sets.
% For two datasets, plotErrorSig should be used. Significance p-values for
% group differences are calculated and printed on the plot. Multiple
% comparisons are not used. Signifance symbols corresponding to particular
% p-values are hard coded. Currently, n.s. will only be shown for group
% statistics.
%
% Inputs:
%       X - X axis labels for plot. Must be the same length as the number
%           of data rows in each cell.
%       data - Data input. Should be a 1xn cell array. Each cell should be
%           a numberical matrix with replicates as columns. NaN values are
%           ignored. Required
%       legLabels - Legend labels.Can be any valid input to the standard
%           legend function.
%       colors - Sets the colors to plot data with. Must be cell array of
%           numerical triplets (0 to 1) of length equal to the number of
%           groups (or bars if there is only 1 group). Optional.
%       pShow - which p-value to plot. Must be 2D numerical array, where
%           each row is a desired pair of bars (based on the single-value
%           index in the matrix Y). Optional.
%       legLabelsShort- Short legend labels. Must be a cell array of
%           strings the same length as pShow
%       mcON - whether multiple corrections will be performed on p values.
%           All possible bar combinations are included in multiple
%           correction. Default is off. Optional.
%       useIdx - which indices to use for ANOVA
%


%% Process inputs

% temporarily turn hold on
tf = ishold; hold on
curFig = gcf;

% get sizes
nGroup = length(data);
nSample = length(X);

if nargin<8 || isempty(mcON)
    mcON = 0;
end

if nargin<9 || isempty(useIdx)
    useIdx = 1:nSample;
end

% calculate plotting values
dataMean = cellfun(@(x) mean(x,2,'omitnan'),data,'UniformOutput',false);
dataSEM = cellfun(@(x) nansem(x,2),data,'UniformOutput',false);


%% Plot data

% plot data
h = cell(1,nGroup);
for ii = 1:nGroup
    h{ii} = errorbar(X,dataMean{ii},dataSEM{ii},'LineWidth',1);

    if nargin>=4 && ~isempty(colors)
        h{ii}.Color = colors{ii};
    end

    if nargin>=5 && ~isempty(styles)
        h{ii}.LineStyle = styles{ii};
    end
end

% set legend
if nargin>=3 && ~isempty(legLabels)
    legend(legLabels,'Location','Southeast')
end


%% Perform pairwise ANOVA calculations

% print p-values for selected data set pairs
if nargin>=6 && ~isempty(pShow)
    %% Print ANOVA p-values

    nCmps = size(pShow,1);
    nCmpsSub = ceil(nCmps^0.5);
    pAnova = nan(nCmps,2);
    pMC = cell(nCmps,1);

    % get p-values
    for cc = 1:nCmps
        data1 = data{pShow(cc,1)}';
        data2 = data{pShow(cc,2)}';

        % calculate multiple comparisons with all data
        [~,pMC{cc}] = anovaRM2W_full_BH(data1,data2,mcON);

        % calculate ANOVA with partial data
        [pACur,~] = anovaRM2W_full_BH(data1(:,useIdx),data2(:,useIdx),mcON);
        pAnova(cc,:) = pACur([1 3])';
    end

    % print p values
    printMultiP(pAnova,legLabelsShort,pShow)

    % get y limits
    yLimsUse = ylim.*[1 1.2];


    %% Generate pairwise subplot

    hhSub = figure; hold on
    for cc = 1:nCmps
        subplot(nCmpsSub,nCmpsSub,cc)

        % set inputs
        curCmp = pShow(cc,:);
        curData1 = data{curCmp(1)}';
        curData2 = data{curCmp(2)}';
        if nargin>=3 && ~isempty(legLabels)
            curLegs = legLabels(curCmp);
        else
            curLegs = [];
        end
        curAnova = nan(nSample,2);
        curAnova(useIdx,:) = pAnova(cc,:).*ones(length(useIdx),2);
        curMC = pMC{cc};
        curColors = colors(curCmp);

        % plot line graph
        plotErrorSig(X,curData1,curData2,curLegs,curAnova,curMC,curColors,0)

        % set limits
        ylim(yLimsUse)
    end
else
    pAnova = NaN;
    pMC = NaN;
    hhSub = NaN;
end

%% Set outputs

% return to previous hold state
figure(curFig);
if tf~=ishold
    hold;
end

% return figure handle only if output is requested
if nargout>0
    hh = h;
end

end

