function plotErrorSigCell(X,data,lgnd,pGroup,pMC,colors,indON)
%%

tf = ishold; hold on

nT = size(data,2);
nPG = size(pGroup);

if nargin<5 || isempty(colors)
    colors = {[0 0 1],[1 0 0]};
end

% set default individual plot state
if nargin<7 || isempty(indON)
    indON = 0;
end

sigSymG = {'n.s.','*','**','***'};
sigSymMC = {'','*','**','***'};
sigThresh = [1,0.05,0.01,0.001];

symGroup = cell(nPG);
for m = 1:nPG(1)
    for n = 1:nPG(2)
        idx = find(pGroup(m,n)<=sigThresh,1,'last');
        symGroup{m,n} = sigSymG(idx);
    end
end

sets = zeros(nPG);
for n = 1:nPG(2)
    [~,~,sets(:,n)] = unique(pGroup(:,n),'stable');
end

symMC = cell(nT,1);
for m = 1:nT
    idx = find(pMC(m,:)<=sigThresh,1,'last');
    symMC{m} = sigSymMC(idx);
end


%% Process data

dataMean = cellfun(@(x) nanmean(x),data);
dataSEM = cellfun(@(x) nansem(x,1),data);


%% Plot means

h = zeros(1,2);
for g = 1:2
    h(g) = errorbar(X,dataMean(g,:),dataSEM(g,:),'color',colors{g},'LineWidth',1);

    % plot individual subject lines
    if indON
        plot(X,cat(2,data{g,:}),'color',[colors{g} 0.25],'LineWidth',0.25)
    end
end


%% Plot p values

yLims = ylim;
yMax = yLims(2)+diff(yLims)*0.2;
yMC = yLims(2)+diff(yLims)*0.05;
yGL = yLims(2)+diff(yLims)*0.12;
yGP = yLims(2)+diff(yLims)*0.15;

ylim([yLims(1) yMax])

if nPG(1)==1
    % explicit for 1 or 2 groups, must be modified if more are added
    if nPG(2)>1
        gText = ['G: ' symGroup{1,1} '; GxT: ' symGroup{1,2}];
    else
        gText = ['G: ' symGroup{1,1}];
    end

    plot([X(1)+0.1 X(end)-0.1],[yGL yGL],'-k','LineWidth',0.5)
    text(median(X),yGP,gText,'HorizontalAlignment','center','FontSize',12)
else
    for ii = 1:max(sets(:,1))
        curX = X(sets(:,1)==ii);
        if length(curX)==1; continue; end
        
        curP = cell(1,nPG(2));
        for n = 1:nPG(2)
            curP(n) = symGroup{find(sets(:,n)==ii,1),n};
        end

        if nPG(2)>1
            gText = ['G: ' curP{1} '   GxT: ' curP{2}];
        else
            gText = ['G: ' curP{1}];
        end
        
        plot([curX(1)+0.1 curX(end)-0.1],[yGL yGL],'-k','LineWidth',0.5)
        text(median(curX),yGP,gText,'HorizontalAlignment','center','FontSize',12)
    end
end

text(X,yMC*ones(1,nT),symMC,'HorizontalAlignment','center','FontSize',12)

if nargin>=3 && ~isempty(lgnd)
    legend(h,lgnd,'Location','Southeast')
end

if tf~=ishold
    hold;
end


end

