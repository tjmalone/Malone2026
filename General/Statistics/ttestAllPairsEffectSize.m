function [tStats] = ttestAllPairsEffectSize(data1,testName,nUnits,pair,MC,testIdxs)
% Performs all pairwise two-tailed Student's t-tests for a set of data.
% Outputs all statistics in a cell array. Can be used to manually generate
% correlation information for compliance with journal manuscript
% requirements. Can perform paired or unparied t-test with or without
% multiple comparisons.
%
% Inputs: (must have same number of points)
%       data1 - first data matrix (subject x category).
%       testName - name of statistical test
%       nUnits - unit of comparsion
%       pair - whether to perform paired or unpaired test (1 or 0). Default
%           is 0.
%       MC - whether to perform correction for multiple comparisons (1 or
%           0). Default is 0.
%       testIdxs - which pairwise column comparions to perform. Default is
%           all combinations
%
% Output:
%       tStats - cell array containing the following corr information
%           {p-value, test name, test statistic, sample size, effect
%           size, confidence interval, degrees of freedom}
%

if nargin<4 || isempty(pair)
    pair = 0;
end

if nargin<5 || isempty(MC)
    MC = 0;
end

if nargin<6 || isempty(testIdxs)
    % number of categories
    nCats = size(data1,2);

    testIdxs = table2array(combinations(1:nCats,1:nCats));
    testIdxs = testIdxs(testIdxs(:,1)<testIdxs(:,2),:);
end

% number of tests
nTests = size(testIdxs,1);


%% Perform t-test

% initialize output variables
pval = zeros(nTests,1);
tstat = zeros(nTests,1);
df = zeros(nTests,1);
cohenD = zeros(nTests,1);
CI = zeros(nTests,2);
szAll = zeros(nTests,2);

for ii = 1:nTests
    % get current data
    if iscell(data1)
        curData1 = cat(1,data1{:,testIdxs(ii,1)});
        curData2 = cat(1,data1{:,testIdxs(ii,2)});
    else
        curData1 = cat(1,data1(:,testIdxs(ii,1)));
        curData2 = cat(1,data1(:,testIdxs(ii,2)));
    end

    % perform t-test
    if pair
        [~,pval(ii),~,stats] = ttest(curData1,curData2);
    else
        [~,pval(ii),~,stats] = ttest2(curData1,curData2);
    end

    % find test statistics
    tstat(ii) = stats.tstat;
    df(ii) = stats.df;

    % find effect size
    effect = meanEffectSize(curData1,curData2,'Effect','cohen','Paired',pair);
    cohenD(ii) = effect{"CohensD","Effect"};
    CI(ii,:) = effect{"CohensD","ConfidenceIntervals"};

    szAll(ii,:) = [length(curData1),length(curData2)];
end

% perform multiple corrections
if MC
    [~,pval] = bonferroni_holm(pval');
end


%% Generate output array

tStats = cell(nTests,7);

for ii = 1:nTests
    % set text values
    pValue = ['p=' num2str(pval(ii),'%.2g')];
    testStatistic = ['t=' num2str(tstat(ii),2)];
    effectSize = ['d=' num2str(cohenD(ii),2)];
    confidenceInterval = ['CI=' num2str(CI(ii,1),2) ',' num2str(CI(ii,2),2)];
    degreesOfFreedom = ['df=' num2str(df(ii))];

    % set test name
    if nargin<2 || isempty(testName)
        testName = '';
    end

    % get unique sample numbers
    szAllUnq = unique(szAll(ii,:),'stable');
    if length(szAllUnq)==1
        szNum = num2str(szAllUnq);
    else
        szNum = [num2str(szAll(ii,1)) ',' num2str(szAll(ii,2))];
    end

    % set sample type
    if nargin<3 || isempty(nUnits)
        nUnits = ' biologically independent samples';
    end

    % set sample size
    sampleSize = ['n=' szNum ' ' nUnits];

    % generate output array
    tStats(ii,:) = {pValue,testName,testStatistic,sampleSize,effectSize,...
        confidenceInterval,degreesOfFreedom};

end
