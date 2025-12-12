function [hh, p] = barGroupOneSample(X,Y,tVal,colors)
%% barGroupOneSample
%
% Creates a bar graph using X and Y inputs. Y must be a two-dimensional
% cell array of one-dimensional numerical arrays. All NaN values will be
% ignored.
%
% Will perform one sample t-test on all bars compared to tVal (default = 0)
%
% Outputs:
%   hh = figure handle
%   p = p values from t-test array
%

tf = ishold; hold on

bMean = cellfun(@(y) mean(y,'omitnan'),Y);
bSEM = cellfun(@(y) std(y,0,'omitnan')/sqrt(sum(~isnan(y))),Y);

h = bar(bMean,'grouped');

% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(bMean);

% Get the x coordinate of the bars
x = nan(nbars,ngroups);
for i = 1:nbars
    x(i,:) = h(i).XEndPoints;
    if nargin>3 && ~isempty(colors)
        h(i).FaceColor = 'flat';
        h(i).CData = colors{i};
    end
end
x = x';

% Plot the errorbars
errorbar(x,bMean,bSEM,'Color',[0 0 0],'LineStyle','none');

xticks(1:ngroups)
xticklabels(X)


%% Find p values

sigSym = {'n.s.','*','**','***'};
sigThresh = [1,0.05,0.01,0.001];

if nargin<3 || isempty(tVal)
    tVal = 0;
end

peakY = max(bMean+bSEM);
ymax = 1.15*diff(ylim);
ylim([min(ylim) ymax])
ystep = 0.05*ymax;

nY = numel(Y);
p = zeros(nY,1);

for ii = 1:nY
    [~,p(ii)] = ttest(Y{ii},tVal);
    
    xCur = x(ii);
    yCur = max([0 peakY]);
    pSym = sigSym(find(p(ii)<=sigThresh,1,'last'));
    
    text(xCur,yCur+2*ystep,pSym,'HorizontalAlignment','center',...
        'FontSize',12)
end


%%

if tf~=ishold
    hold;
end

if nargout>0
    hh = h;
end

end

