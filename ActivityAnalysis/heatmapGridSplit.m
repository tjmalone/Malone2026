function [h,mapP] = heatmapGridSplit(dataAll,labelX,labelY,valueSign,idxX,idxY,barLocs)
%%  h = heatmapGrid(dataAll,labelX,labelY)
% Generates a heatmap to display the difference between a set of data by
% two variables independently and combined. Intended for use with sex,
% morphology, and sex/morphlogy with AD project data. Can be used for
% alternative comparisons if input data format is followed.
%
% Inputs:
%   dataAll - data sorted by variable X and Y. Must be 3 x 3 x 2 cell
%       array. Dimension 1 corresponds to variable X choices. Dimension 2
%       corresponds to variable Y choices. Variable 3 corresponds to
%       primary group difference (i.e. genotype).
%   labelX and labelY - labels for variable groups and elements. Must be 1
%       x 4 cell array. First element is variable name. Other elements are
%       elements names.
%   valueSign - if a higher WT than PS19 value corresponse to AD deficit,
%       equals 1. Otherwise, equal -1. Default is 1.
%   idxX - defines which heatmap X indices to show in each of two plots.
%       Must be 1xN cell array. Default is {1:3,2:3}
%   idxY - defines which heatmap Y indices to show in each of two plots.
%       Must be 1xN cell array. Default is {1,2:3}
%   barLocs - Where colorbars should be located. Must be 1xN cell array
%       containing valid positions. Default is
%       {'southoutside', 'eastoutside'}
%

if nargin<2 || isempty(labelX)
    labelX = {'Variable X','X1','X2','X3'};
end

if nargin<3 || isempty(labelY)
    labelY = {'Variable Y','Y1','Y2','Y3'};
end

if nargin<4 || isempty(valueSign)
    valueSign = 1;
elseif abs(valueSign)~=1
    error('Invalid value sign')
end

if nargin<5 || isempty(idxX)
    idxX = {1:3,2:3};
end

if nargin<6 || isempty(idxY)
    idxY = {1,2:3};
end

if nargin<7 || isempty(barLocs)
    barLocs = {'southoutside','eastoutside'};
end

nTypes = length(idxX);
if length(idxY)~=nTypes || length(barLocs)~=nTypes
    error('Invalid input sizes')
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

%% Calculate group differences and p values

% calculate current difference
mapMean = cellfun(@(x) mean(x,'omitnan'),dataAll);
mapDiff = -valueSign*diff(mapMean,[],3);

% calculate current p-values
[~,mapP] = cellfun(@(x,y) ttest2(x,y),dataAll(:,:,1),dataAll(:,:,2));

% find difference signs
mapSign = mapDiff>0;

% define significance thresholds
sigSym = {'n.s.','*','**','***'};
sigThresh = [1,0.05,0.01,0.001];

syms = cell(3);
for m = 1:3
    for n = 1:3
        idx = find(mapP(m,n)<=sigThresh,1,'last');
        if ~isempty(idx)
            syms{m,n} = sigSym{idx};
        else
            syms{m,n} = '';
        end
    end
end


%% Plot difference heatmap

% initialize figure axes
h = figure;
tiledlayout(nTypes,2);
ax = cell(nTypes,2);

% loop through groupings
for hg = 1:nTypes
    nX = length(idxX{hg});
    nY = length(idxY{hg});
    curXRange = [0.5 nX+0.5];
    curYRange = [0.5 nY+0.5];

    for aa = 1:2
        % set axis
        ax{aa,hg} = nexttile(); hold on
        if aa==1
            title('Difference Heatmap')
            plotData = mapDiff(idxX{hg},idxY{hg});
        elseif aa==2
            title('p-value Heatmap')
            plotData = -log10(mapP(idxX{hg},idxY{hg})).*(mapSign(idxX{hg},idxY{hg})*2-1);
        end

        % transpose matrix
        plotData = plotData';

        % make heat map with transparency
        imAlpha = ones(size(plotData));
        imAlpha(isnan(plotData))=0;
        imagesc(plotData,'AlphaData',imAlpha,'YData',[size(plotData,1) 1]);

        % set background color
        set(gca,'color',0*[1 1 1]);

        col = colorbar('Location',barLocs{hg});

        if aa==1
            % set color limits
            diffMax = max(abs(plotData),[],'all','omitnan');
            if isnan(diffMax) || diffMax==0
                diffMax = 1;
            end

            curCLim = [-diffMax diffMax];
            
            % set color map
            % curCMap = customcolormap([0 0.25 0.5 0.75 1],[0 0 1; 0.5 0.5 1; 1 1 1; 1 0.5 0.5; 1 0 0]);
            % curCMap = customcolormap([0 0.25 0.5 0.75 1],[0 0.75 0; 0.5 1 0.5; 1 1 1; 1 0.5 1; 0.75 0 0.75]);
            curCMap = customcolormap([0 0.25 0.5 0.75 1],[0.75 0 0.75; 1 0.5 1; 1 1 1; 0.5 1 0.5; 0 0.75 0]);

            % set colorbar label
            ylabel(col,'AD deficit')
        else
            % set color limits
            curCLim = [-4,4];

            % set color map
            cBar = [0.0501 0.05];
            cBarMap = abs(log10(cBar))/range(curCLim);
            % curCMap = customcolormap([0 0.5-cBarMap(1) 0.5-cBarMap(2) 0.5+cBarMap(2) 0.5+cBarMap(1) 1],...
            %     [0 0 1; 1 0.95 0.95; 1 1 1; 1,1,1; 1 0.95 0.95; 1 0 0],1000);
            % curCMap = customcolormap([0 0.5-cBarMap(1) 0.5-cBarMap(2) 0.5+cBarMap(2) 0.5+cBarMap(1) 1],...
            %     [0 0.75 0; 0.95 1 0.95; 1 1 1; 1,1,1; 1 0.95 1; 0.75 0 0.75],1000);
            curCMap = customcolormap([0 0.5-cBarMap(1) 0.5-cBarMap(2) 0.5+cBarMap(2) 0.5+cBarMap(1) 1],...
                [0.75 0 0.75; 1 0.95 1; 1 1 1; 1,1,1; 0.95 1 0.95; 0 0.75 0],1000);

            % set colorbar label
            ylabel(col,'p-value (log10)')

            % print significance symbols
            curX = idxX{hg}-idxX{hg}(1)+1;
            curY = fliplr(idxY{hg}-idxY{hg}(1)+1);
            for xx = 1:nX
                for yy = 1:nY
                    tX = idxX{hg}(xx);
                    tY = idxY{hg}(yy);
                    curText = [syms{tX,tY} ', p=' num2str(mapP(tX,tY),2)];
                    text(curX(xx),curY(yy),curText,'FontSize',16,...
                        'HorizontalAlignment','center','VerticalAlignment','middle')
                end
            end
        end

        % set colors
        clim(curCLim)
        ax{aa,hg}.Colormap = colormap(ax{aa,hg},curCMap);

        % set plot limits
        axis('equal')
        xlim(curXRange)
        ylim(curYRange)
        set(gca,'FontSize',16)
        set(col,'FontSize',16)

        % set x axis labels
        ax{aa,hg}.XAxisLocation = 'top';
        ax{aa,hg}.XTick = 1:nX;
        ax{aa,hg}.XTickLabel = labelX(idxX{hg}+1);
        xlabel(labelX{1})

        % set y axis labels
        ax{aa,hg}.YTick = 1:nY;
        ax{aa,hg}.YTickLabel = labelY(fliplr(idxY{hg})+1);
        ax{aa,hg}.YTickLabelRotation = 90;
        ylabel(labelY{1});
    end
end

end

