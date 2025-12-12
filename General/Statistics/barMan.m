function [hh,p] = barMan(X,M,SEM,N,STD,colors,pShow,mcON)
%% barGroup
%
% Creates a bar graph using mean, std, sem, and n directly.
% X and Y inputs. Valid with a two-dimensional matrix.
% Columns of a matrix will be treated as groups. All NaN values will be
% ignored. A t-test is performed on every combination of bars. By default,
% the t-test will be unpaired without multiple comparisons.
%
% Inputs:
%       X - x-axis labels. Must be a valid input to xticklabels.
%       M,SEM,N,STD - matrix containing means, sems, ns, and stds, of bars.
%           Current implementation requires these to be 2D numerical
%           arrays.
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
%       indON - Whether to plot individual data points.
%
% Outputs:
%       hh - figure handle
%       p - p value matrix from t-test
%

%% Process inputs

if ~isequal(size(M),size(SEM),size(N),size(STD))
    error('Error: data inputs must be the same size')
elseif ~isequal(size(M,2),length(X))
    error('Error: Size of X must match size of data')
end

% temporarily turn hold on
tf = ishold; hold on


%% Plot graph

% plot bar graph
h = bar(M,'grouped');

% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(M);

% Get the x coordinate of the bars
x = nan(ngroups,nbars);

% process by group if there are multiple groups, otherwise process by bar
for ii = 1:nbars
    % get x values
    x(:,ii) = h(ii).XEndPoints;

    % set color
    if nargin>5 && ~isempty(colors)
        % set color
        h(ii).FaceColor = 'flat';
        h(ii).CData = colors{ii};
    end
end

% plot errorbars
errorbar(x,M,SEM,'Color',[0 0 0],'LineStyle','none');

% add x labels
xticks(1:ngroups)
xticklabels(X)


%% Find p values

% initialize p value matrix
nY = numel(M);
p = zeros(nY);

% perform t-tests for all bar combinations
for ii = 1:nY
    for jj = 1:nY
        p(ii,jj) = ttestMan(M(ii),STD(ii),N(ii),M(jj),STD(jj),N(jj));
    end
end

% perform multiple corrections
if nargin>7 && ~isempty(mcON) && mcON==1
    
    % extract lower diagonal (unique p values)
    m  = tril(true(size(p)),-1);
    pVec  = p(m);

    % perform multiple corrections
    [~,pMCVec] = bonferroni_holm(pVec);

    % generate corrected matrix
    pMC = zeros(size(p));
    pMC(m) = pMCVec;
    pMC = pMC + pMC';
    % end

else
    % generate uncorrected matrix
    pMC = p;
end


%% Plot p values

if nargin>6 && ~isempty(pShow)
    % significance definition
    sigSym = {'n.s.','*','**','***'};
    sigThresh = [1,0.05,0.01,0.001];
    
    % set new y limits
    ymax = 1.15*diff(ylim);
    ylim([min(ylim) ymax])
    
    % get bar peaks and set y buffer above bars
    peakY = M+SEM;
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

