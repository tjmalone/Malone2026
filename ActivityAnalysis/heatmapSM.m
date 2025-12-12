function h = heatmapSM(dataAll,dataY,dataX,dataYX,labelY,labelX)
%%  h = heatmapSM(dataSex,dataMorph,dataMix)
% generates a heatmap to display the difference between a set of data by
% two variables independently and combined. Intended for use with sex,
% morphology, and sex/morphlogy with AD project data. Can be used for
% alternative comparisons if input data format is followed.
%
% Inputs:
%   dataX - data sorted by variable X. Must be 2 x 2 cell array.
%       Dimension 1 corresponds to variable X choices. Dimension 2
%       corresponds to primary group difference (i.e. genotype)
%   dataY - data sorted by variable Y. Must be 2 x 2 cell array.
%       Dimension 1 corresponds to variable X choices. Dimension 2
%       corresponds to primary group difference (i.e. genotype)
%   dataXY - data sorted by variable X. Must be 2 x 2 x 2 cell array.
%       Dimension 1 corresponds to variable X choices. Dimension 1
%       corresponds to variable Y choices. Dimension 3 corresponds to
%       primary group difference (i.e. genotype)
%   labelX and labelY - labels for variable groups and elements. Must be 1
%       x 3 cell array. First element is variable name. Other elements are
%       elements names
%

if nargin<4 || isempty(labelY)
    labelY = {'Variable Y','Y1','Y2'};
end
if nargin<5 || isempty(labelX)
    labelX = {'Variable X','X1','X2'};
end

%%

% clear; close all; clc
% 
% dataY = {[0.9 .9 .9 1],[2 2 2 2];[1 1 1 1],[1 1.1 1 1]};
% dataX = {[1 1 1 1],[3 3 3 3];[1 1 1 .11],[-1 -1 -1 -1]};
% dataYX = cat(3,dataY,dataX);
% labelY = {'Variable Y','Y1','Y2'};
% labelX = {'Variable X','X1','X2'};


%% Calculate group differences and p values

% initialize output arrays
mapDiff = nan(3,3);
mapP = nan(3,3);

for dd = 1:4
    % define current data
    if dd==1
        curData = dataAll;
        dim = 2;
    elseif dd==2
        curData = dataY;
        dim = 2;
    elseif dd==3
        curData = dataX;
        dim = 2;
    elseif dd==4
        curData = dataYX;
        dim = 3;
    end

    % calculate current difference
    curMean = cellfun(@(x) mean(x,'omitnan'),curData);
    curDiff = -diff(curMean,[],dim);

    % calculate current p-values
    if ~any(cellfun(@isempty,curData),'all')
        if dim==2
            [~,curP] = cellfun(@(x,y) ttest2(x,y),curData(:,1),curData(:,2));
        elseif dim==3
            [~,curP] = cellfun(@(x,y) ttest2(x,y),curData(:,:,1),curData(:,:,2));
        end
    elseif dim==2
        curP = nan(size(curData,1),1);
    elseif dim==3
        curP = nan(2,2);
    end

    % store current data
    if dd==1
        iR = 1;
        iC = 1;
    elseif dd==2
        iR = 2:3;
        iC = 1;
    elseif dd==3
        iR = 1;
        iC = 2:3;
    elseif dd==4
        iR = 2:3;
        iC = 2:3;
    end
    mapDiff(iR,iC) = curDiff;
    mapP(iR,iC) = curP;
end

% correct map positions

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
        diffMax = max(abs(plotData),[],'all');
        curCLim = [-diffMax diffMax];

        % set color map
        % curCMap = customcolormap([0 0.25 0.5 0.75 1],[0 0 1; 0.5 0.5 1; 1 1 1; 1 0.5 0.5; 1 0 0]);
        % curCMap = customcolormap([0 0.25 0.5 0.75 1],[0 1 0; 0.5 1 0.5; 1 1 1; 1 0.5 1; 1 0 1]);
        % curCMap = customcolormap([0 0.25 0.5 0.75 1],[0 0.7 0; 0.5 0.75 0.5; 1 1 1; 0.75 0.5 0.75; 0.75 0 0.75]);
        curCMap = customcolormap([0 0.25 0.5 0.75 1],[0 0.75 0; 0.4 0.75 0.4; 1 1 1; 0.75 0.4 0.75; 0.75 0 0.75]);
        curCMap = customcolormap([0 0.25 0.5 0.75 1],[0 0.75 0; 0.5 1 0.5; 1 1 1; 1 0.5 1; 0.75 0 0.75]);


        % set colorbar label
        ylabel(col,'Difference (WT-AD)')
    else
        % set color limits
        curCLim = [-4,4];

        % set color map
        cBar = [0.051 0.05];
        cBarMap = abs(log10(cBar))/range(curCLim);
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
    ax{aa}.XTick = [2 3];
    ax{aa}.XTickLabel = labelX(2:3);
    xlabel(labelX{1})

    % set y axis labels
    ax{aa}.YTick = [1 2];
    ax{aa}.YTickLabel = labelY(3:-1:2);
    ax{aa}.YTickLabelRotation = 90;
    ylabel(labelY{1});
end

