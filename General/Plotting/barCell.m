function [hh, p] = barCell(X,Y,varargin)
%% barCell
%
% Creates a bar graph using X and Y inputs. Y must be a one-dimensional
% cell array of one-dimensional numerical arrays or a two-dimensional
% numeric arrray. Numerical arrays will be averaged across the first
% dimension. All NaN values will be ignored.
%
% If varargin is "pair", a paried t-test will be performed on all sets
% pariwise. Otherwise, an unpaired t-test will be performed.
%
% Outputs:
%   hh = figure handle
%   p = p values from t-test array
%

tf = ishold; hold on

if iscell(Y)
    bMean = cellfun(@(y) mean(y,'omitnan'),Y);
    bSEM = cellfun(@(y) std(y,0,'omitnan')/sqrt(sum(~isnan(y))),Y);
else
    bMean = mean(Y,1,'omitnan');
    bSEM = std(Y,0,1,'omitnan')./sqrt(sum(~isnan(Y),1));
end

h = bar(bMean);
errorbar(X,bMean,bSEM,'Color',[0 0 0],'LineStyle','none')

if tf~=ishold
    hold;
end

nY = size(Y,2);
p = zeros(nY);

if ~isempty(varargin) && strcmp(varargin{1},'pair')
    pair = 1;
else
    pair = 0;
end
    
for ii = 1:nY
    for jj = 1:nY
        if iscell(Y)
            if pair
                [~,p(ii,jj)] = ttest(Y{ii},Y{jj});
            else
                [~,p(ii,jj)] = ttest2(Y{ii},Y{jj});
            end
        else
            if pair
                [~,p(ii,jj)] = ttest(Y(:,ii),Y(:,jj));
            else
                
                [~,p(ii,jj)] = ttest2(Y(:,ii),Y(:,jj));
            end
        end
    end
end

if nargout>0
    hh = h;
end

end


