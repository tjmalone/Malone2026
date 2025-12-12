function [corrStats] = corrEffectSize(data1,data2,testName,nUnits)
% Performs two-tailed pearson linear correlation. Outputs all statistics in
% a cell array. Can be used to manually generate correlation information
% for compliance with journal manuscript requirements. If only one data
% matrix is input, a correlation of the mean of this matrix against time
% will be calculated. 
%
% Inputs: (must have same number of points)
%       data1 - first data matrix (subject x 1).
%           Alternate: (repeat x time). Mean and transpose are taken
%       data2 - second data matrix (subject x 1)
%       testName - name of statistical test
%       nUnits - unit of comparsion
%
% Output:
%       corrStats - cell array containing the following corr information
%           {p-value, test name, test statistic, sample size, effect
%           size, confidence interval, degrees of freedom}
%

if size(data1,2)>1
    data1 = mean(data1,1,'omitnan')';
end

if isempty(data2)
    data2 = (1:length(data1))';
end

% remove NaN values
nanidx = isnan(data1) | isnan(data2);
data1 = data1(~nanidx);
data2 = data2(~nanidx);


%% Perform Correlation

% perform correlation
useIdx = 2;
[rMat,pMat,rL,rU] = corrcoef(data1,data2,'Rows','complete');

% find correlation parameters
rho = rMat(useIdx);
pval = pMat(useIdx);

% find effect size
rsquare = rho^2;

% find confidence interval
CI = [rL(useIdx),rU(useIdx)].^2;


%% Generate output array

% set p-value
pValue = ['p=' num2str(pval,'%.2g')];

% set test name
if nargin<3 || isempty(testName)
    testName = '';
end

% set test statistic
testStatistic = ['r=' num2str(rho,2)];

% set sample number
szNum = num2str(length(data1));

% set sample type
if nargin<4 || isempty(nUnits)
    nUnits = ' biologically independent samples';
end

% set sample size
sampleSize = ['n=' szNum ' ' nUnits];

% set effect size
effectSize = ['r^2=' num2str(rsquare,2)];

% set confidence interval
confidenceInterval = ['CI=' num2str(CI(1),2) ',' num2str(CI(2),2)];

% set degrees of freedom
degreesOfFreedom = ['df=' num2str(length(data1)-2)];

% generate output array
corrStats = {pValue,testName,testStatistic,sampleSize,effectSize,...
    confidenceInterval,degreesOfFreedom};

end

