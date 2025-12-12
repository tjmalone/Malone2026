function [lossAll,mdlAll] = fitlme_crossval(tbl,formula,cv,errorType)
%% fitlme_crossval
% Peforms k-way cross-validation of a fitlme model. Output the loss and the
% individual submodels.
%
% Inputs:
%       tbl - data in a table format (required)
%       formula - model formula (required)
%       cv - 
%           if not cell: cross-validation partition
%           if cell: defines cross-validaiton partition (default)
%               k - number of validation folds (default = 5)
%               stratGroups - stratification variable. Length must equal
%               number of rows in table. (default is no stratification)
%       errorType - loss function for cross-validaiton. MSE (default) or
%           adjusted MAE
%


%% Process inputs

% define partition settings if not set
if nargin<3 || isempty(cv)
    k = 5;
    stratGroup = size(tbl,1);
    cv = {k,stratGroup};
end

% create partition if cv is cell (k, stratGroups)
if iscell(cv)
    cv = cvpartition(cv{2},'KFold',cv{1});
else
    k = cv.NumTestSets;
end

if nargin<4 || isempty(errorType)
    errorType = 'MSE';
end

% ititialize model and loss
mdlAll = cell(k,1);
lossAll = zeros(k,1);


%% Perform cross-validation

for kk = 1:k
    % define training and test indices
    idxTrain = training(cv,kk);
    idxTest = test(cv,kk);
    
    % define training and test table
    tblTrain = tbl(idxTrain,:);
    tblTest = tbl(idxTest,:);

    % perform linear fit
    mdl = fitlme(tblTrain,formula,'FitMethod','ML');
    yTest = tblTest.(mdl.ResponseName);
    yFit = predict(mdl,tblTest);

    % calculate MAE and MSE
    MAE = mean(abs(yTest-yFit),'omitnan');
    adjMAE = MAE/range(yFit);
    MSE = mean((yTest - yFit).^2,'omitnan');

    % store loss variable
    if strcmp(errorType,'MAE')
        lossAll(k) = adjMAE;
    elseif strcmp(errorType,'MSE')
        lossAll(kk) = MSE;
    else
        error('Invalid loss function')
    end

    % store model
    mdlAll{kk} = mdl;
end

end

