function plotLineP(X,pGroup,pMC,sigSymG,sigSymMC,sigThresh,rawP,pgLabels)
%% plotLineP
% print group p-values (i.e. Anova) and multiple comparison p-values (i.e.
% t-test) on a previously generated line plot.


%% Process inputs

% save hold state
tf = ishold; hold on

nT = length(X);
nPG = size(pGroup);

% process X input
if sum(size(X)>1)>2
    error('X input must be 1-dimensional')
elseif ~isnumeric(X)
    X = 1:length(X);
end

% symbol key
if nargin<4 || isempty(sigSymG)
    sigSymG = {'n.s.','*','**','***'};
end

if nargin<5 || isempty(sigSymMC)
    sigSymMC = {'','*','**','***'};
end

if nargin<6 || isempty(sigThresh)
    sigThresh = [1,0.05,0.01,0.001];
end

if nargin<7 || isempty(rawP)
    rawP = 0;
end

if nargin<8 || isempty(pgLabels)
    pgLabels = {'G','GxT'};
elseif ischar(pgLabels)
    pgLabels = {pgLabels};
end


%% Set symbols

% set group symbols
symGroup = cell(nPG);
for m = 1:nPG(1)
    for n = 1:nPG(2)
        if rawP==0
            idx = find(pGroup(m,n)<=sigThresh,1,'last');
            symGroup{m,n} = sigSymG{idx};
        else
            symGroup{m,n} = num2str(pGroup(m,n),2);
        end
    end
end

% combine group symbols by unique p value (proxy for "true" group)
sets = zeros(nPG);
for n = 1:nPG(2)
    [~,~,sets(:,n)] = unique(pGroup(:,n),'stable');
end

symMC = cell(nT,1);
if nargin>=3 && ~isempty(pMC)
    for m = 1:nT
        if rawP==0
            idx = find(pMC(m,:)<=sigThresh,1,'last');
            symMC{m} = sigSymMC{idx};
        else
            symMC{m} = num2str(pMC(m,:),2);
        end
    end
end


%% Print p-values

yLims = ylim;
yMax = yLims(2)+diff(yLims)*0.2;
yMC = yLims(2)+diff(yLims)*0.05;
yGL = yLims(2)+diff(yLims)*0.12;
yGP = yLims(2)+diff(yLims)*0.15;

ylim([yLims(1) yMax])

if nPG(1)==1
    % explicit for 1 or 2 groups, must be modified if more are added
    if nPG(2)>1
        gText = [pgLabels{1} ': ' symGroup{1,1} '; ' pgLabels{2} ': ' symGroup{1,2}];
    else
        gText = [pgLabels{1} ': ' symGroup{1,1}];
    end

    plot([X(1)+0.1 X(end)-0.1],[yGL yGL],'-k','LineWidth',0.5)
    text(median(X),yGP,gText,'HorizontalAlignment','center',...
        'VerticalAlignment','baseline','FontSize',12)
else
    for ii = 1:max(sets(:,1))
        curX = X(sets(:,1)==ii);
        if length(curX)==1; continue; end
        
        curP = cell(1,nPG(2));
        for n = 1:nPG(2)
            curP(n) = symGroup{find(sets(:,n)==ii,1),n};
        end

        if nPG(2)>1
            gText = [pgLabels{1} ': ' curP{1} '   ' pgLabels{2} ': ' curP{2}];
        else
            gText = [pgLabels{1} ': ' curP{1}];
        end
        
        plot([curX(1)+0.1 curX(end)-0.1],[yGL yGL],'-k','LineWidth',0.5)
        text(median(curX),yGP,gText,'HorizontalAlignment','center',...
            'VerticalAlignment','baseline','FontSize',12)
    end
end

text(X,yMC*ones(1,nT),symMC,'HorizontalAlignment','center','FontSize',12)

% restore hold state
if tf~=ishold
    hold;
end