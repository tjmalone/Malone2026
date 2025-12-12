function [pAnova,pMC] = anovaRM2W_BH(data1,data2)
% Performs repeated measure 2-Way ANOVA assuming sphericity. Calculates
% multiple comparisons using a non-paired ttest with Bonferonni correction.
% Anova p-vaule matches Prism output for predictor variable (X/Learning
% Type). 
%
% Inputs: (must have same number of time points)
%       data1 - first data matrix (subject x time)
%       data2 - second data matrix (subject x time)
%

% concatenate data
data = cat(1,data1,data2);

% generate predictor categories (X)
X = cat(1,repmat({'A'},size(data1,1),1),repmat({'B'},size(data2,1),1));

% generate data table
Q = [table(X) array2table(data)];

% generate repeated measures table 
W = table((1:size(data,2))','VariableNames',{'Time'});  

% create repeated measures model
rm = fitrm(Q,['data1-data' num2str(size(data,2)) '~X'],'WithinDesign',W);

% run ANOVA
pAll = ranova(rm,'WithinModel','Time');
pAnova = pAll{'X','pValue'};

% run ttest with multiple corrections
[~,MC] = ttest2(data1,data2);
[~,pMC] = bonferroni_holm(MC');


