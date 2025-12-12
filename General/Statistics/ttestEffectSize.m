function [tStats] = ttestEffectSize(data1,data2,testName,nUnits,pair,MC,limitP,testType)
% Performs two-tailed Student's t-test. Outputs all statistics in a cell
% array. Can be used to manually generate correlation information for
% compliance with journal manuscript requirements. Can perform paired or
% unparied t-test with or without multiple comparisons.
%
% Inputs: (must have same number of points)
%       data1 - first data matrix (subject x category)
%       data2 - second data matrix (subject x category)
%       testName - name of statistical test
%       nUnits - unit of comparsion
%       pair - whether to perform paired or unpaired test (1 or 0). Default
%           is 0.
%       MC - whether to perform correction for multiple comparisons (1 or
%           0). Default is 0.
%       limitP - whether to only print significant p-values (1 or 0).
%           Default is 0.
%       testType - what type of test to perform. Default is "t-test".
%           "ranksum" can be performed as an unpaired test.
%
% Output:
%       tStats - cell array containing the following corr information
%           {p-value, test name, test statistic, sample size, effect
%           size, confidence interval, degrees of freedom}
%

%% Note: Currently cannot perform multiple t-tests.

if nargin<5 || isempty(pair)
    pair = 0;
end

if nargin<6 || isempty(MC)
    MC = 0;
end

if nargin<7 || isempty(limitP)
    limitP = 0;
end

if nargin<8 || isempty(testType)
    testType = 'ttest';
end


% correct data input format
nTests = size(data1,2);
if nTests>1 && ~iscell(data1)
    data1 = num2cell(data1,1);
    data2 = num2cell(data2,1);
elseif nTests==1 && iscell(data1)
    if numel(data1)~=1 || numel(data2)~=1
        data1 = data1';
        data2 = data2';
        nTests = size(data1,2);
        % error('Single t-test must be performed using numerical array')
    else
        data1 = data1{1};
        data2 = data2{1};
        nTests = size(data1,2);

        if nTests>1
            data1 = num2cell(data1,1);
            data2 = num2cell(data2,1);
        end
    end
end

if size(data1,2) ~= size(data2,2)
    error('Mismatched data sizes')
end


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
        curData1 = cat(1,data1{:,ii});
        curData2 = cat(1,data2{:,ii});
    else
        curData1 = data1;
        curData2 = data2;
    end

    % remove NaN values
    curData1 = curData1(~isnan(curData1));
    curData2 = curData2(~isnan(curData2));

    % perform t-test
    if pair
        if strcmp(testType,'ttest')
            [~,pval(ii),~,stats] = ttest(curData1,curData2);
        else
            error('Paired test must be t-test')
        end
    else
        if strcmp(testType,'ttest')
            [~,pval(ii),~,stats] = ttest2(curData1,curData2);
            % [~,pval(ii),~,stats] = ttest2(curData1,curData2,'tail','right');
        elseif strcmp(testType,'ranksum')
            [pval(ii),~,stats] = ranksum(curData1,curData2);
        end
    end

    if strcmp(testType,'ttest')
        % find test statistics
        tstat(ii) = stats.tstat;
        df(ii) = stats.df;

        % find effect size
        effect = meanEffectSize(curData1,curData2,'Effect','cohen','Paired',pair);
        cohenD(ii) = effect{"CohensD","Effect"};
        CI(ii,:) = effect{"CohensD","ConfidenceIntervals"};
    elseif strcmp(testType,'ranksum')
        % find test statistics
        tstat(ii) = stats.ranksum;

        df(ii) = NaN;
        cohenD(ii) = NaN;
        CI(ii,:) = NaN;
    end

    szAll(ii,:) = [length(curData1),length(curData2)];
end

% perform multiple corrections
if MC
    [~,pval] = bonferroni_holm(pval');
end


%% Generate output array


% initialize text variables
if strcmp(testType,'ttest')
    testStatistic = 't=';
elseif strcmp(testType,'ranksum')
    testStatistic = 'ranksum=';
end
pValue = 'p=';
effectSize = 'd=';
confidenceInterval = 'CI=';
degreesOfFreedom = 'df=';

% get unique sample numbers
szAllUnq = unique(szAll,'stable');
if length(szAllUnq)==1
    szNum = num2str(szAllUnq);
else
    szNum = [];
end

% define significant p-values
if limitP
    useIdx = find(pval<=0.05);
else
    useIdx = 1:nTests;
end

for ii = 1:nTests
    % skip non-significant p-values
    if ~ismember(ii,useIdx)
        continue
    end

    % set text values
    pValue = [pValue num2str(pval(ii),'%.2g')];
    if strcmp(testType,'ttest')
        testStatistic = [testStatistic num2str(tstat(ii),2)];
    elseif strcmp(testType,'ranksum')
        testStatistic = [testStatistic num2str(tstat(ii))];
    end
    effectSize = [effectSize num2str(cohenD(ii),2)];
    confidenceInterval = [confidenceInterval...
        num2str(CI(ii,1),2) ',' num2str(CI(ii,2),2)];
    degreesOfFreedom = [degreesOfFreedom num2str(df(ii))];

    % add delimiters
    if ii~=nTests
        pValue = [pValue ', '];
        testStatistic = [testStatistic ', '];
        effectSize = [effectSize ', '];
        confidenceInterval = [confidenceInterval '; '];
        degreesOfFreedom = [degreesOfFreedom ', '];
    end

    % set sample numbers
    if length(szAllUnq)~=1
        szNum = [szNum num2str(szAll(ii,1)) ',' num2str(szAll(ii,2))];
        if ii~=nTests
            szNum = [szNum '; '];
        end
    end
end

% set test name
if nargin<3 || isempty(testName)
    testName = '';
end

% set sample type
if nargin<4 || isempty(nUnits)
    nUnits = ' biologically independent samples';
end

% set sample size
sampleSize = ['n=' szNum ' ' nUnits];

% generate output array
tStats = {pValue,testName,testStatistic,sampleSize,effectSize,...
    confidenceInterval,degreesOfFreedom};

end
