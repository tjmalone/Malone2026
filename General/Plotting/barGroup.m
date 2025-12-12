function [hh,p] = barGroup(X,Y,plotType,colors,pShow,tType,mcON,indON,ovrMean,ovrSEM,varargin)
%% barGroup
%
% Creates a bar graph using X and Y inputs. Valid with a two-dimensional
% cell array of one-dimensional numerical arrays or two-dimensional matrix.
% Columns of a matrix will be treated as groups. All NaN values will be
% ignored. A t-test is performed on every combination of bars. By default,
% the t-test will be unpaired without multiple comparisons.
%
% Inputs:
%       X - x-axis labels. Must be a valid input to xticklabels.
%       Y - input data. Must be a 2D cell array of 1D numerical arrays or a
%           2D matrix. In a 2D cell array each row is a group of clustered
%           bars.
%       type - plot type. Must be bar, boxplot, or violin. Default is bar.
%       colors - Sets the colors to plot data with. Must be cell array of
%           numerical triplets (0 to 1) of length equal to the number of
%           groups (or bars if there is only 1 group). Optional.
%       pShow - which p-value to plot. Must be 2D numerical array, where
%           each row is a desired pair of bars (based on the single-value
%           index in the matrix Y). Optional.
%       tType - whether a paired ('pair') or unpaired (default) t-test is performed.
%           Optional.
%       mcON - whether multiple corrections will be performed on p values.
%           All possible bar combinations are included in multiple
%           correction. Default is off. Optional.
%       indON - Whether to plot individual data points. (0=off, 1=on, 2=on
%           with jitter)
%       ovrMean/ovrSEM - overrides mean and sem calculated from Y. p-values
%           are still based on Y. Both must be specified together to
%           overide Y. This can allow plotting based on a transform of Y,
%           while preserving p-values.
%      varargin - Used to set optional inputs for violinPlots.
%
% Outputs:
%       hh - figure handle
%       p - p value matrix from t-test
%
% Taylor Malone, Yi Gu Lab, 2025
%

%% Process data

% set plot type
if nargin<3 || isempty(plotType)
    plotType = 'bar';
elseif ~ismember(plotType,{'bar','boxplot','violin'})
    error('Invalid plot type')
end

% set default individual plot state
if nargin<8 || isempty(indON)
    indON = 0;
end

% temporarily turn hold on
tf = ishold; hold on

% convert input to cell array
if ~iscell(Y)
    numInput = 1;
    if ndims(Y)==1
        error('Incorrect inputs')
    end
    Ytemp = cell(1,size(Y,2));
    for ii = 1:size(Y,2)
        Ytemp{ii} = Y(:,ii);
    end
    Y = Ytemp';
else
    numInput = 0;
end

% get data means and error
if nargin<10 || isempty(ovrMean) || isempty(ovrSEM)
    bMean = cellfun(@(y) mean(y,'omitnan'),Y);
    bSEM = cellfun(@(y) std(y,0,'omitnan')/sqrt(sum(~isnan(y))),Y);
else
    bMean = ovrMean;
    bSEM = ovrSEM;
end


%% Plot data

% Calculate the number of groups and number of bars in each group
if numInput==0
    [ngroups,nbars] = size(bMean);
elseif numInput==1
    [nbars,ngroups] = size(bMean);
end

% plot data
switch plotType
    case 'bar'
        h = bar(bMean,'grouped');

        % Get the x coordinate of the bars
        x = nan(ngroups,nbars);

        % process by group if there are multiple groups, otherwise process by bar
        for ii = 1:nbars
            if ngroups>1
                % get x values
                x(:,ii) = h(ii).XEndPoints;

                % set color
                if nargin>3 && ~isempty(colors)
                    h(ii).FaceColor = 'flat';
                    if min(size(colors))==1
                        h(ii).CData = colors{ii};
                    else
                        for gg = 1:ngroups
                            h(ii).CData(gg,:) = colors{gg,ii};
                        end
                    end
                end
            else
                % get x values
                x(:,ii) = h.XEndPoints(ii);

                % set color
                if nargin>3 && ~isempty(colors)
                    h.FaceColor = 'flat';
                    h.CData(ii,:) = colors{ii};
                end
            end
        end

        % plot errorbars
        errorbar(x,bMean,bSEM,'Color',[0 0 0],'LineStyle','none');

        % calculate bar widths
        if ngroups>1
            groupSpacing = mean(diff(x(1,:)));
            barWidth = h(1).BarWidth*groupSpacing;
        else
            barWidth = 0.8;
        end

        % plot idividual values
        if indON~=0
            for ii = 1:nbars
                for jj = 1:ngroups
                    curX = x(jj,ii);
                    curY = Y{jj,ii};
                    % define jitter
                    if indON==1
                        XJitter = 'None';
                    elseif indON==2
                        XJitter = 'density';
                    end

                    % plot points
                    swarmchart(curX*ones(length(curY),1),curY,30,'filled',...
                        'MarkerEdgeColor','None','MarkerFaceColor','k',...
                        'XJitter',XJitter,'XJitterWidth',barWidth*0.4)
                end
            end
        end

        % set peak values
        peakY = bMean+bSEM;

    case 'boxplot'
        error('Not implemented')

    case 'violin'
        if nargin>3
            vColors = colors;
        else
            vColors = [];
        end

        [hAll,x] = violinplotWrapper(Y,vColors,indON,varargin{:});
        h = [hAll(:).ViolinPlot];

        % set peak values
        peakY = cellfun(@(x) max(x,[],'omitnan'),Y);

        % set axis properties
        box off
end


%% Set labels

% add x labels
if ngroups>1
    xticks(1:ngroups)
else
    xticks(1:nbars)
end
xticklabels(X)


%% Find p values

% initialize p value matrix
nY = numel(Y);
p = zeros(nY);

% set t-test type
if nargin>5 && strcmp(tType,'pair')
    pair = 1;
else
    pair = 0;
end

% perform t-tests for all bar combinations
flag = 0;
for ii = 1:nY
    for jj = 1:nY
        if pair
            if isequal(size(Y{ii}),size(Y{jj}))
                [~,p(ii,jj)] = ttest(Y{ii},Y{jj});
            else
                % perform unpaired t-test if sizes are different (will
                % catch "mixed" t-tests unless multiple groups have the
                % same number of elements (should only be used very
                % carefully)
                [~,p(ii,jj)] = ttest2(Y{ii},Y{jj});
                if flag==0
                    disp('Unpaired test used for some comparisons')
                    flag = 1;
                end
            end
        else
            [~,p(ii,jj)] = ttest2(Y{ii},Y{jj});
        end
    end
end

% perform multiple corrections
if nargin>6 && ~isempty(mcON) && mcON==1
    % extract lower diagonal (unique p values)
    m  = tril(true(size(p)),-1);
    pVec  = p(m);

    % perform multiple corrections
    [~,pMCVec] = bonferroni_holm(pVec);

    % generate corrected matrix
    pMC = zeros(size(p));
    pMC(m) = pMCVec;
    pMC = pMC + pMC';
else
    % generate uncorrected matrix
    pMC = p;
end


%% Plot p values

if nargin>4 && ~isempty(pShow)
    % significance definition
    sigSym = {'n.s.','*','**','***'};
    sigThresh = [1,0.05,0.01,0.001];

    % set new y limits
    ymax = min(ylim) + 1.15*diff(ylim);
    ylim([min(ylim) ymax])

    % get bar peaks and set y buffer above bars
    ystep = 0.05*ymax;

    globYs = [];

    % cycle through each bar pair
    for ii = 1:size(pShow,1)

        % get current x, y, and p value
        xCur = x(pShow(ii,:));
        yCur = max([0 max(peakY(pShow(ii,:)))]);
        pCur = pMC(pShow(ii,1),pShow(ii,2));

        % set significance text
        pSym = sigSym(find(pCur<=sigThresh,1,'last'));

        % set y positions to not overlap
        yCurThis = yCur+ystep;
        while ismember(yCurThis,globYs)
            yCurThis = yCurThis+ystep;
        end
        globYs(end+1) = yCurThis;

        % plot significance line
        plot(xCur,yCurThis*ones(1,2),'k')

        % plot significance text
        text(mean(xCur),yCurThis,pSym,'HorizontalAlignment',...
            'center','VerticalAlignment','bottom','FontSize',12)
    end

end


%% Set outputs

% return to previous hold state
if tf~=ishold
    hold;
end

% return figure handle only if output is requested
if nargout>0
    hh = h;
end

end

