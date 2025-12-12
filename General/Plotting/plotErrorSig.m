function hh = plotErrorSig(X,data1,data2,legLabels,pGroup,pMC,colors,indON)
%% plotErrorSig
% Plot a line graph with error bars with significance labels for two groups
% of data. Signifance symbols corresponding to particular p-values are hard
% coded. Currently, n.s. will only be shown for group statistics.
%
% Inputs:
%       X - X axis labels for plot. Must be the same length as the number
%           of data rows.
%       data1 - Data input 1. Must be numberical matrix with replicates as
%           columns. NaN values are ignored.
%       data2 - Data input 2. Must be numberical matrix with replicates as
%           columns. NaN values are ignored.
%       legLabels - Legend labels.Can be any valid input to the standard
%           legend function. Can be skipped with empty array
%       pGroup - Group difference p-values (i.e. ANOVA). Must be 2D
%           numerical array. Musy have the same number of rows as data
%           inputs. Each element should be the p-value for the
%           corresponding significance group. Significance groups are
%           determined by matching p-values, so groups may be combined if
%           p-values are identical across groups. Columns correspond to the
%           type of signifiance group. Currently hard coded so that column
%           1 (required) is the group difference between data inputs and
%           column 2 (optional) is the interaction between group difference
%           and row variable (i.e. time). Groups must be the same across
%           columns
%       pMC - Multiple comparisons p values. Must be a 1D numerical array
%           the same length as the numbner of data rows.
%       colors - Sets the colors to plot data groups with. Must be cell
%           array of numerical triplets (0 to 1). Optional.
%       indON - Whether to plot individual subject lines. Individual lines
%           will be plotted with transparency, but the transparency will
%           not be kept if a .fig is repopened. Saving as a .tif will
%           preserve transparency.
% 


%% Process inputs

% temporarily turn hold on
tf = ishold; hold on

% get significance value sizes
nT = length(pMC);
nPG = size(pGroup);

% set default colors
if nargin<7 || isempty(colors)
    colors = {[0 0 1],[1 0 0]};
end

% set default individual plot state
if nargin<8 || isempty(indON)
    indON = 0;
end

% significance sybols and ranges
sigSymG = {'n.s.','*','**','***'};
sigSymMC = {'','*','**','***'};
sigThresh = [1,0.05,0.01,0.001];

% find group significance symbols
symGroup = cell(nPG);
for m = 1:nPG(1)
    for n = 1:nPG(2)
        idx = find(pGroup(m,n)<=sigThresh,1,'last');
        symGroup{m,n} = sigSymG(idx);
    end
end

% find unique significance groups
[~,~,sets] = unique(pGroup(:,n),'stable');

% find multiple comparisons symbols
symMC = cell(nT,1);
for m = 1:nT
    idx = find(pMC(m,:)<=sigThresh,1,'last');
    symMC{m} = sigSymMC(idx);
end

% threshold for shaded errorbars
shadeThresh = 25;


%% Process data

% concatenate data
data = {data1,data2};
dataMean = zeros(2,nT);
dataSEM = zeros(2,nT);

% get data means and error
for g = 1:2
    dataMean(g,:) = nanmean(data{g},1);
    dataSEM(g,:) = nansem(data{g},1);
end


%% Plot means

h = zeros(1,2);
for g = 1:2
    if length(X)<shadeThresh
        % plot groups data with errorbars
        h(g) = errorbar(X,dataMean(g,:),dataSEM(g,:),'color',colors{g},'LineWidth',1);
    else
        % plot groups data with shading
        h(g) = semshade(data{g},0.3,colors{g},X);
    end

    % plot individual subject lines
    if indON
        plot(X,data{g}','color',[colors{g} 0.25],'LineWidth',0.25)
    end
end


%% Plot p values

% set new y limits to make space for significance labels
yLims = ylim;
yMax = yLims(2)+diff(yLims)*0.2;    % new max y
yMC = yLims(2)+diff(yLims)*0.05;    % y for MC p value
yG = yLims(2)+diff(yLims)*0.12;    % y for group p value and line

ylim([yLims(1) yMax])

% cycle through multiple significance groups
for ii = 1:max(sets)
    % get current group set
    curX = find(sets==ii);
    if length(curX)==1; continue; end
    
    % get symbol for each group
    curP = cell(1,nPG(2));
    for n = 1:nPG(2)
        curP(n) = symGroup{curX(1),n};
    end
    
    % set significance text
    if nPG(2)>1
        gText = ['G: ' curP{1} '   GxT: ' curP{2}];
    else
        gText = ['G: ' curP{1}];
    end
    
    % plot group significance bar and text
    plot([X(curX(1))+0.1 X(curX(end))-0.1],[yG yG],'-k','LineWidth',0.5)
    text(median(X(curX)),yG,gText,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','FontSize',12)
end

% plot multiple comparisons significance text
text(X,yMC*ones(1,nT),symMC,'HorizontalAlignment','center','FontSize',12)

% set legend
if ~isempty(legLabels)
    legend(h,legLabels,'Location','Southeast')
end

% return to previous hold state
if tf~=ishold
    hold;
end

% return figure handle only if output is requested
if nargout>0
    hh = h;
end

end

