function mdlFinal = stepwiseModel_post(tblTrain,tblTest,idxFinal,idxOther,idxY,modelType,nTrials,verbose,skipVal)
%% stepwiselme_postTest
% Run after determining a final model to evaluate factor interaction
% contributions. Uses training and test data to define r/R2. Can work with
% both 'lme' and 'glm' models.

%% Initialize inputs

% define default model
if nargin<6 || isempty(modelType)
    modelType = 'lme';
end

% define default trial number
if nargin<7 || isempty(nTrials)
    nTrials = 19;
end

% reset trial number for lme models
if strcmp(modelType,'lme')
    nTrials = 1;
end

% define whether to print results
if nargin<8 || isempty(verbose)
    verbose = 1;
end

% define whether to skip sub model and expansion model
if nargin<9 || isempty(skipVal)
    skipVal = 0;
end

% get variable names
labsFull = string(tblTrain.Properties.VariableNames(1:end-1));
labsFinal = string(tblTrain.Properties.VariableNames(idxFinal));
labsOther = string(tblTrain.Properties.VariableNames(idxOther));
labsY = string(tblTrain.Properties.VariableNames(idxY));

% generate final model formula
formFinal = labsY + ' ~ ' + strjoin(labsFinal,' + ');

% run full model
if strcmp(modelType,'lme')
    mdl = fitlme(tblTrain,formFinal);
elseif strcmp(modelType,'glm')
    mdl = fitglm(tblTrain,formFinal,'Distribution','binomial','BinomialSize',nTrials);
else
    error('Invalid model type')
end

% find base model correlation
yPredTrain = predict(mdl,tblTrain)*nTrials;
yPredTest = predict(mdl,tblTest)*nTrials;
rBase = corr(tblTest.y,yPredTest);

% store final model
mdlFinal = struct();
mdlFinal.mdl = mdl;
mdlFinal.formula = formFinal;
mdlFinal.tblTrain = tblTrain;
mdlFinal.tblTrain = tblTest;
mdlFinal.yTrain = tblTrain.y;
mdlFinal.yPredTrain = yPredTrain;
mdlFinal.yTest = tblTest.y;
mdlFinal.yPredTest = yPredTest;
mdlFinal.r = rBase;
mdlFinal.idxGlob = idxFinal;

% skip end of function for a minimal run
if skipVal==1
    return
end


%% Run sub models for final factors

% print sub modeling message
if verbose
    fprintf('\nFinal factor sub-models: r_base = %.2f\n',mdlFinal.r);
end

% split formula for sub-modeling
nFactors = length(labsFinal);

% identify sub-model factors
factMat = logical(~eye(nFactors));

% initialize p-value array
nTrain = size(tblTrain,1);
nTest = size(tblTest,1);
factorRem = cell(1,nFactors);
yPredRemTrain = zeros(nTrain,nFactors);
rRemTrain = zeros(1,nFactors);
yPredRemTest = zeros(nTest,nFactors);
rRemTest = zeros(1,nFactors);

for mm = 1:nFactors
    % create formula
    formSub = 'y ~ ' + strjoin(labsFinal(factMat(mm,:)), ' + ');
    factorRem{mm} = labsFinal{mm};

    % run sub model
    if strcmp(modelType,'lme')
        mdlSub = fitlme(tblTrain,formSub);
    elseif strcmp(modelType,'glm')
        mdlSub = fitglm(tblTrain,formSub,'Distribution','binomial','BinomialSize',nTrials);
    end

    % calculate prediction and R for training data
    curYPredTrain = predict(mdlSub,tblTrain)*nTrials;
    rCurTrain = corr(tblTrain.y,curYPredTrain);
    yPredRemTrain(:,mm) = curYPredTrain;
    rRemTrain(mm) = rCurTrain;

    % calculate prediction and R for test data
    curYPredTest = predict(mdlSub,tblTest)*nTrials;
    rCurTest = corr(tblTest.y,curYPredTest);
    yPredRemTest(:,mm) = curYPredTest;
    rRemTest(mm) = rCurTest;

    % print sub-model info
    if verbose
        fprintf('%d. Factor %s, r = %.2f\n',mm,erase(labsFinal{mm},' '),rCurTest);
    end
end

% store factor removal statistics
mdlFinal.factorRem = string(factorRem);
mdlFinal.yPredRemTrain = yPredRemTrain;
mdlFinal.rRemTrain = rRemTrain;
mdlFinal.yPredRemTest = yPredRemTest;
mdlFinal.rRemTest = rRemTest;


%% Run expansion models for removed factors

% print expansion modeling message
if verbose
    fprintf('\nExtra factor expansion-models: r_base = %.2f\n',mdlFinal.r);
end

% identify expansion-model factors
factMat = logical(eye(length(labsFull)));
factMat(:,idxFinal) = 1;
factMat(idxFinal,:) = [];
nFactors = size(factMat,1);

% initialize p-value array
yPredExpTrain = zeros(nTrain,nFactors);
rExpTrain = zeros(nFactors,1);
yPredExpTest = zeros(nTest,nFactors);
rExpTest = zeros(nFactors,1);

for mm = 1:nFactors
    % create formula
    formExp = 'y ~ ' + strjoin(labsFull(factMat(mm,:)), ' + ');

    % run expansion model
    if strcmp(modelType,'lme')
        mdlExp = fitlme(tblTrain,formExp);
    elseif strcmp(modelType,'glm')
        mdlExp = fitglm(tblTrain,formExp,'Distribution','binomial','BinomialSize',nTrials);
    end

    % calculate prediction and R for training data
    curYPredTrain = predict(mdlExp,tblTrain)*nTrials;
    rCurTrain = corr(tblTrain.y,curYPredTrain);
    yPredExpTrain(:,mm) = curYPredTrain;
    rExpTrain(mm) = rCurTrain;

    % calculate prediction and R for test data
    curYPredTest = predict(mdlExp,tblTest)*nTrials;
    rCurTest = corr(tblTest.y,curYPredTest);
    yPredExpTest(:,mm) = curYPredTest;
    rExpTest(mm) = rCurTest;

    % print expansion-model info
    if verbose
        fprintf('%d. Factor %s, r = %.2f\n',mm,erase(labsOther{mm},' '),rCurTest);
    end
end

if verbose
    fprintf('\n')
end

% store factor expansion statistics
mdlFinal.yPredExpTrain = yPredExpTrain;
mdlFinal.rExpTrain = rExpTrain;
mdlFinal.yPredExpTest = yPredExpTest;
mdlFinal.rExpTest = rExpTest;
mdlFinal.factorExp = labsOther';

end

