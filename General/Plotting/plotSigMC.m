function plotSigMC(pMC,xMC,yMC)
% plot significance for multiple comparisons

% get significance value sizes
if nargin<2
    xMC = 1:length(pMC);
end
if nargin<3
    yMC = 0.95*max(ylim);
end

nT = length(pMC);

% significance sybols and ranges
sigSymMC = {'','*','**','***'};
sigThresh = [1,0.05,0.01,0.001];

% find multiple comparisons symbols
symMC = cell(nT,1);
for m = 1:nT
    idx = find(pMC(m,:)<=sigThresh,1,'last');
    symMC{m} = sigSymMC(idx);
end


%% Plot p values

% plot multiple comparisons significance text
text(xMC,yMC*ones(1,nT),symMC,'HorizontalAlignment','center','FontSize',12)
