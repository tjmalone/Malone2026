function h = heatmapGrid(dataAll,labelY,labelX)
%%  h = heatmapGrid(dataAll,labelY,labelX)
% Generates a heatmap to display the difference between a set of data by
% two variables independently and combined. Intended for use with sex,
% morphology, and sex/morphlogy with AD project data. Can be used for
% alternative comparisons if input data format is followed.
%
% Inputs:
%   dataAll - data sorted by variable X and Y. Must be 3 x 3 x 2 cell
%       array. Dimension 1 corresponds to variable Y choices. Dimension 2
%       corresponds to variable X choices. Variable 3 corresponds to
%       primary group difference (i.e. genotype).
%   labelX and labelY - labels for variable groups and elements. Must be 1
%       x 4 cell array. First element is variable name. Other elements are
%       elements names.
%

if nargin<2 || isempty(labelY)
    labelY = {'Variable Y','Y1','Y2','Y3'};
end
if nargin<3 || isempty(labelX)
    labelX = {'Variable X','X1','X2','X3'};
end

%%
% 
% clear; close all; clc
% 
% data1 = num2cell(zeros(3));
% data1{1,1} = [0 0];
% data2 = {-1:-1 -2 -3;-4 1 2;3 4 0:0};
% dataAll = cat(3,data1,data2);
% 
% labelY = {'Variable Y','Y1','Y2','Y3'};
% labelX = {'Variable X','X1','X2','X3'};


%% Calculate group differences and p values

% calculate current difference
mapMean = cellfun(@(x) mean(x,'omitnan'),dataAll);
mapDiff = -diff(mapMean,[],3);

% calculate current p-values
[~,mapP] = cellfun(@(x,y) ttest2(x,y),dataAll(:,:,1),dataAll(:,:,2));

% find difference signs
mapSign = mapDiff>0;


%% Plot difference heatmap

% initialize figure axes
h = figure;
tiledlayout(1, 2);
ax = cell(1,2);

for aa = 1:2
    % set axis
    ax{aa} = nexttile(); hold on
    if aa==1
        title('Difference Heatmap')
        plotData = mapDiff;
    elseif aa==2
        title('p-value Heatmap')
        plotData = -log10(mapP).*(mapSign*2-1);
    end

    % make heat map with transparency
    imAlpha = ones(size(plotData));
    imAlpha(isnan(plotData))=0;
    imagesc(plotData,'AlphaData',imAlpha,'YData',[size(plotData,1) 1]);

    % set background color
    set(gca,'color',0*[1 1 1]);

    col = colorbar;

    if aa==1
        % set color limits
        diffMax = max(abs(plotData),[],'all','omitnan');
        curCLim = [-diffMax diffMax];

        % set color map
        curCMap = customcolormap([0 0.25 0.5 0.75 1],[0 0 1; 0.5 0.5 1; 1 1 1; 1 0.5 0.5; 1 0 0]);
        curCMap = customcolormap([0 0.25 0.5 0.75 1],[0 0.75 0; 0.5 1 0.5; 1 1 1; 1 0.5 1; 0.75 0 0.75]);

        % set colorbar label
        ylabel(col,'Difference (WT-AD)')
    else
        % set color limits
        curCLim = [-4,4];

        % set color map
        cBar = [0.0501 0.05];
        cBarMap = abs(log10(cBar))/range(curCLim);
        curCMap = customcolormap([0 0.5-cBarMap(1) 0.5-cBarMap(2) 0.5+cBarMap(2) 0.5+cBarMap(1) 1],...
            [0 0 1; 1 0.95 0.95; 1 1 1; 1,1,1; 1 0.95 0.95; 1 0 0],1000);
        curCMap = customcolormap([0 0.5-cBarMap(1) 0.5-cBarMap(2) 0.5+cBarMap(2) 0.5+cBarMap(1) 1],...
            [0 0.75 0; 0.95 1 0.95; 1 1 1; 1,1,1; 1 0.95 1; 0.75 0 0.75]);

        % set colorbar label
        ylabel(col,'p-value (log10)')
    end

    % set colors
    clim(curCLim)
    ax{aa}.Colormap = colormap(ax{aa},curCMap);

    % set plot limits
    axis('equal')
    xlim([0.5 3.5])
    ylim([0.5 3.5])
    set(gca,'FontSize',16)
    set(col,'FontSize',16)

    % set x axis labels
    ax{aa}.XAxisLocation = 'top';
    ax{aa}.XTick = 1:3;
    ax{aa}.XTickLabel = labelX(2:4);
    xlabel(labelX{1})

    % set y axis labels
    ax{aa}.YTick = 1:3;
    ax{aa}.YTickLabel = labelY(4:-1:2);
    ax{aa}.YTickLabelRotation = 90;
    ylabel(labelY{1});
end

end

