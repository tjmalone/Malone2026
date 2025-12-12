function [anovaStats] = anovaEffectSizeMinimal(data1,data2,nUnits)
% Performs repeated measure 2-Way ANOVA assuming sphericity or linear
% mixed-effects model. Outputs all statistics in a cell array. Can be used
% to manually generate ANOVA information for compliance with journal
% manuscript requirements. Currently only the group difference is output.
% Outputs minimal information to be compatible with linear mixed-effects
% models. Use anovaEffectSize for full capabilities.
%
% Inputs: (must have same number of time points)
%       data1 - first data matrix (as for anovaRM2W_full_BH)
%       data2 - second data matrix (as for anovaRM2W_full_BH)
%       nUnits - unit of comparsion
%
% Output:
%       anovaStats - cell array containing the following ANOVA information
%           {p-value, test name, test statistic, sample size, effect
%           size, confidence interval, degrees of freedom}
%
% Note: effect size and confidence interval are NaN for compatibility with
% other functions.
%

%% Perform ANOVA

% remove full NaN rows
data1 = data1(~all(isnan(data1),2),:);
data2 = data2(~all(isnan(data2),2),:);

[pAnova,~,~,stats] = anovaRM2W_full_BH(data1,data2,[],1);


%% Generate output array

% set p-value
pValue = ['p=' num2str(pAnova(1),'%.2g')];

% set test name
testName = stats.test;

% set test statistic
testStatistic = ['F=' num2str(stats.F(1),2)];

% set sample number
szAll = unique([size(data1,1),size(data2,1)],'stable');
if isscalar(szAll)
    szNum = num2str(szAll);
else
    szNum = [num2str(szAll(1)) ',' num2str(szAll(2))];
end

% set sample type
if nargin<3 || isempty(nUnits)
    nUnits = 'biologically independent samples';
end

% set sample size
sampleSize = ['n=' szNum ' ' nUnits];

% set effect size
if strcmp(testName,'two-way repeated measures ANOVA')
    % get effect size and confidence interval for ANOVA
    tempStats = anovaEffectSize(data1,data2,testName,nUnits);
    effectSize = tempStats{5};
    confidenceInterval = tempStats{6};
else
    % set lme effect size and confidence interval for ANOVA
    effectSize = ['eta=' num2str(stats.eta2,'%.2g')];
    confidenceInterval = ['CI=' num2str(stats.CI(1),'%.2g') ',' num2str(stats.CI(2),'%.2g')];
end

% set degrees of freedom
dfAll = stats.df(1,:);
if isscalar(dfAll)
    dfNum = num2str(dfAll);
else
    dfNum = [num2str(dfAll(1),'%.0f') ',' num2str(dfAll(2),'%.0f')];
end
degreesOfFreedom = ['df=' dfNum];

% generate output array
anovaStats = {pValue,testName,testStatistic,sampleSize,effectSize,...
    confidenceInterval,degreesOfFreedom};

end

