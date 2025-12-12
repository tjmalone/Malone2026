function [pAnova,pLabel] = anovaRM3W(data,varLabels)
% Performs repeated measure 3-Way ANOVA assuming sphericity.Has not been
% tested against Prism outputs.
%
% Inputs: (must have same number of elemments in the within vars variables)
%       data - cell array containing data. Each element contains the data
%       matrix for one group (V1 x V2 x X).
%       withinVars - labels for the repeated measures variables that are 
%
% Outputs:
%       pANOVA - p value for ANOVA for X, T, X-T, X-T alternate. Matlab has
%           two ways to produce an interaction p-value. The first of these
%           matches the output of Prism.
%       pLabel - labels for the Anova p values.
%

%% Define data

catGroup = {'A','B'};

% generate predictor categories (X)
allData = [];
allX = {};

newN = size(data{1},1)*size(data{1},2);

for ii = 1:2
    % reshape data
    curData = data{ii};
    curReshape = reshape(curData,newN,[])';

    % define groups type
    curX = repmat(catGroup(ii),size(curReshape,1),1);
    
    % combine variables
    allData = cat(1,allData,curReshape);
    allX = cat(1,allX,curX);
end

% define variable 1 groups
baseV1 = repmat((1:size(data{1},1))',1,size(data{1},2));
allV1 = reshape(baseV1,1,newN)';

% define variable 1 groups
baseV2 = repmat(1:size(data{1},2),size(data{1},1),1);
allV2 = reshape(baseV2,1,newN)';


%% Generate model

% generate data table
Q = [table(allX) array2table(allData)];

% generate repeated measures table 
if nargin<2
    factorNames = {'V1','V2'};
else
    factorNames = varLabels;
end
W = table(allV1,allV2,'VariableNames',factorNames);

% create repeated measures model
rm = fitrm(Q,['allData1-allData' num2str(size(allData,2)) '~allX'],'WithinDesign',W);

% run ANOVA
pAll = anova(rm,'WithinModel',[factorNames{1} '*' factorNames{2}]);

% % run my repeated measures anova here
% [ranovatbl] = ranova(rm,'WithinModel','V1*V2');

% save output p values
pAnova = pAll{[2 5 8 11],'pValue'};
pLabel = {'X',['X*' factorNames{1}],['X*' factorNames{2}],['X*' factorNames{1} '*' factorNames{2}]};

end
