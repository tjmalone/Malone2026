function [domDataTrain,yPredTestTier] = dominanceAnalysis(tblTrain,factors,fitParam,modelType,tblTest)
%% dominance_analysis
% performs dominance analsysis on a linear mixed effects model or a general
% linear model based on the algorithm described at
% https://github.com/dominance-analysis/dominance-analysis. R2 values are
% calculated directly based on training data. If a test data is provided,
% predictions are returned for test data.
%
% Inputs:
%   tblTrain - table containing all data. Column names should be the factor
%       names. y must be the name of the output variable. There can be
%       extra columns as long as all input factors are included.
%   factors - factors for defining full model. Can be a cell array
%       containing the fixed and random factors or can be a string/char
%       array defining the full model formula.
%   fitParam - parameter for model run.
%       If modelType is lme: which mixed effects fit method will be used
%           for models. 'ML' (default) or 'REML'. If 'REML' the pLRT is
%           only valid for random factors.
%       If modelType is glm: the number of trials that response variable is
%           based on. Default = 19
%   tblTest - table containing test data.
%
% Outputs:
%   domDataTrain - struct containing the various R2 dominance calculations.
%       tier: contains the R2 contributions for each factor (columns) to
%           each tier of models (rows).
%       domInd: individual dominance, variance explained by each factor
%           alone. Contribution of random factor will be NaN. Analogous to
%           pariwise correlation.
%       domPart: average partial dominance, average R2 contribution to all
%           models except full model and individual model (least useful)
%       domInter: interactional dominiance: R2 contribution to the full
%           model (semi-partial R2). Measures "unique" contribution to full
%           model.
%       domTot: total dominance, summarizes R2 contribution of each factor
%           to all subset models.
%       pLRT: significance of the unique contribution of each factor based
%           on the likelihood ratio test. Should primarily be used to
%           interpret random factors rather than fixed factors.
%   yPredTier - test predictions for each tier of dominance analysis.
%       Analagous to the "tier" variable in "domDataTrain"
%


%% Process inputs

if ~iscell(factors)
    splitForm = strsplit(char(factors),{'+','~'});
    factors = splitForm(2:end);
end

% define model and model parameters
if nargin<4 || isempty(modelType)
    modelType = 'lme';
end

% define model parameter
if nargin<3 || isempty(fitParam)
    if strcmp(modelType,'lme')
        fitParam = 'ML';
    elseif strcmp(modelType,'glm')
        fitParam = 19;
    else
        error('Invalid model type')
    end
end

if strcmp(modelType,'lme')
    predictScale = 1;
elseif strcmp(modelType,'glm')
    predictScale = fitParam;
end

if nargin>=5 && ~isempty(tblTest)
    nTest = size(tblTest,1);
else
    nTest = 0;
end


%% Determine factor/row matching

% generate factor matrix
nFactors = length(factors);
factMat = logical(dec2bin(0:2^nFactors-1) - '0');
nMods = size(factMat,1);

factorData = struct();
for ff = 1:nFactors
    % current factor
    factorData(ff).factor = factors{ff};

    % current include rows
    rowInc = find(factMat(:,ff)==1);
    nRows = length(rowInc);

    % current match rows
    rowMatch = zeros(nRows,1);
    for rr = 1:nRows
        rowSub = factMat(rowInc(rr),:);
        rowSub(ff) = 0;
        rowMatch(rr) = find(ismember(factMat,rowSub,'rows'),1);
    end

    % get row tier
    rowTier = sum(factMat(rowInc,:),2);

    % store row info
    factorData(ff).rowInc = rowInc;
    factorData(ff).rowMatch = rowMatch;
    factorData(ff).rowTier = rowTier;
end


%% Generate all models

modelsAll = cell(nMods,1);
R2 = zeros(nMods,1);
yPredTest = cell(nMods,1);

for mm = 1:nMods
    % create formula
    formCur = ['y ~' strjoin(factors(factMat(mm,:)), '+')];
    try
        % run model
        if strcmp(modelType,'lme')
            modelsAll{mm} = fitlme(tblTrain,formCur,'FitMethod',fitParam);
        elseif strcmp(modelType,'glm')
            modelsAll{mm} = fitglm(tblTrain,formCur,'Distribution','binomial','BinomialSize',fitParam);
        end

        % store fit quality
        R2(mm) = modelsAll{mm}.Rsquared.Ordinary;

        % predict test data model
        if nTest~=0
            yPredTest{mm} = predict(modelsAll{mm},tblTest)*predictScale;
        else
            yPredTest{mm} = nan(nTest,1);
        end
    catch
        % define null and error model outputs
        modelsAll{mm} = NaN;
        if strcmp(formCur,'y ~')
            R2(mm) = 0;
            yPredTest{mm} = zeros(nTest,1);
        else
            R2(mm) = NaN;
            yPredTest{mm} = nan(nTest,1);
        end
    end
end


%% Process R2 values

% intialize R2 difference matrix
R2_tier = zeros(nFactors);

% initialize unique contribution p-value
pLRT = zeros(1,nFactors);

for ff = 1:nFactors
    % extract R2 values
    R2_inc = R2(factorData(ff).rowInc);
    R2_match = R2(factorData(ff).rowMatch);

    % calculate R2 matched differences
    R2_diff = R2_inc-R2_match;
    factorData(ff).R2_diff = R2_diff;

    % calculate tier averages
    R2_tierCur = zeros(nFactors,1);
    for gg = 1:nFactors
        R2_tierCur(gg) = mean(R2_diff(factorData(ff).rowTier==gg),'omitnan');
    end

    % store tier averages
    factorData(ff).R2_tier = R2_tierCur;
    R2_tier(:,ff) = R2_tierCur;

    if strcmp(modelType,'lme')
        % identify full model and match
        idxFull = find(factorData(ff).rowTier==gg);
        mdlFull = modelsAll{factorData(ff).rowInc(idxFull)};
        mdlSub = modelsAll{factorData(ff).rowMatch(idxFull)};

        % compare sub model p-values
        LRT = compare(mdlSub,mdlFull);
        pLRT(ff) = LRT.pValue(2);
    elseif strcmp(modelType,'glm')
        pLRT(ff) = NaN;
    end
end

% calculate individual dominance
R2_DomInd = R2_tier(1,:);

% calculate partial dominance
R2_DomPart = mean(R2_tier(2:end-1,:),1,'omitnan');

% calculate interactional dominance
R2_DomInter = mean(R2_tier(end,:),1,'omitnan');

% calculate total dominance
R2_DomTot = mean(R2_tier,1,'omitnan');

% store R2 values
domDataTrain = struct();
domDataTrain.tier = R2_tier;
domDataTrain.domInd = R2_DomInd;
domDataTrain.domPart = R2_DomPart;
domDataTrain.domInter = R2_DomInter;
domDataTrain.domTot = R2_DomTot;
domDataTrain.pLRT = pLRT;


%% Process test predictions

% intialize yPred difference cell
yPredTestTier = cell(nFactors);

if nTest~=0
    for ff = 1:nFactors
        % extract R2 values
        yPred_inc = yPredTest(factorData(ff).rowInc);
        yPred_match = yPredTest(factorData(ff).rowMatch);

        % calculate R2 matched differences
        yPred_cat = [yPred_inc,yPred_match];
        factorData(ff).yPred_cat = yPred_cat;

        % calculate tier averages
        yPred_tierCur = cell(nFactors,1);
        for gg = 1:nFactors
            yPred_tierCur{gg} = yPred_cat(factorData(ff).rowTier==gg,:);
        end

        % store tier averages
        factorData(ff).yPred_tier = yPred_tierCur;
        yPredTestTier(:,ff) = yPred_tierCur;
    end
end

end

