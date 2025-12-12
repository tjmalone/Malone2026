function printMultiP(pGroup,lgnd,cmps,sigSymG,sigSymMC,sigThresh)
%% plotLineP
% print group p-values (i.e. Anova) and multiple comparison p-values (i.e.
% t-test) on a previously generated line plot.


%% Process inputs

% save hold state
tf = ishold; hold on

nPG = size(pGroup);

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


%% Set symbols

% set group symbols
symGroup = cell(nPG);
for m = 1:nPG(1)
    for n = 1:nPG(2)
        idx = find(pGroup(m,n)<=sigThresh,1,'last');
        symGroup{m,n} = sigSymG{idx};
    end
end


%% Print p-values

pText = cell(nPG(1),1);
for ii = 1:nPG(1)
    pText{ii} = [lgnd{cmps(ii,1)} '-' lgnd{cmps(ii,2)} ': ' symGroup{ii,1} ', ' symGroup{ii,2}];
end

yLims = ylim;
XLims = xlim;

text(XLims(1),yLims(2),pText,'HorizontalAlignment','left',...
    'VerticalAlignment','top','FontSize',12)

% restore hold state
if tf~=ishold
    hold;
end