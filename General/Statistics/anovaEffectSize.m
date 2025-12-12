function [anovaStats] = anovaEffectSize(data1,data2,testName,nUnits,isDep)
% Performs repeated measure 2-Way ANOVA assuming sphericity. Outputs all
% statistics in a cell array. Can be used to manually generate ANOVA
% information for compliance with journal manuscript requirements.
% Currently only the group difference is output.
%
% Inputs: (must have same number of time points)
%       data1 - first data matrix (subject x time)
%       data2 - second data matrix (subject x time)
%       testName - name of statistical test
%       nUnits - unit of comparsion
%       isDep - repeated measures variables. If input data are cells (i.e.
%           unequal numbers of rows), default is no repeated measure. If
%           input data are numerical (i.e. equal number of rows) default is
%           repeated measure across columns (i.e. time).
%
% Output:
%       anovaStats - cell array containing the following ANOVA information
%           {p-value, test name, test statistic, sample size, effect
%           size, confidence interval, degrees of freedom}
%

%% Perform ANOVA

if iscell(data1)
    % make data arrays
    nTime = size(data1,2);
    day1 = cell(nTime,1);
    day2 = cell(nTime,1);
    dataM1 = cell(nTime,1);
    dataM2 = cell(nTime,1);

    % concatenate by time point
    for ii = 1:nTime
        dataM1{ii} = cat(1,data1{:,ii});
        dataM2{ii} = cat(1,data2{:,ii});

        dataM1{ii}(isnan(dataM1{ii})) = [];
        dataM2{ii}(isnan(dataM2{ii})) = [];

        day1{ii} = ii*ones(size(dataM1{ii}));
        day2{ii} = ii*ones(size(dataM2{ii}));
    end

    % concatenate data
    cat1 = cat(1,dataM1{:});
    cat2 = cat(1,dataM2{:});

    % concatenate group ids
    num1 = 1*ones(size(cat1));
    num2 = 2*ones(size(cat2));
    dayCat1 = cat(1,day1{:});
    dayCat2 = cat(1,day2{:});

    % combine concatenations
    dataL = [cat1;cat2];
    numCat = [num1;num2];
    dayCatAll = [dayCat1;dayCat2];
    G = cat(2,numCat,dayCatAll);

    % set repeated measure
    if nargin<5 || isempty(isDep)
        isDep = [0,0];
    end

else
    data1(any(isnan(data1),2),:) = [];
    data2(any(isnan(data2),2),:) = [];

    % put all data in 1D array
    data = cat(1,data1,data2);
    dataL = data(:);

    % generate group ID in 1D array
    G1 = cat(1,ones(size(data1,1),1),repmat(2,size(data2,1),1));
    G1L = repmat(G1,size(data,2),1);

    % generate time point ID in 1D array
    G2L = zeros(length(dataL),1);
    for ii = 1:size(data,2)
        G2L(size(data,1)*(ii-1)+1:size(data,1)*(ii)) = ii;
    end

    % generate group matrix
    G = [G1L,G2L];

    % set repeated measure
    if nargin<5 || isempty(isDep)
        isDep = [0,1];
    end
end

% perform ANOVA
[statsAll,tblAll] = mes2way(dataL,G,'eta2','isDep',isDep,'nBoot',10000);


%% Generate output array

% set p-value
pValue = ['p=' num2str(tblAll{3,6},'%.2g')];

% set test name
if nargin<3 || isempty(testName)
    testName = '';
end

% set test statistic
testStatistic = ['F=' num2str(tblAll{3,5},2)];

% set sample number
if iscell(data1)
    szNum = [];
    for ii = 1:nTime
        szNum = [szNum num2str(length(dataM1{ii})) ',' num2str(length(dataM2{ii}))];

        if ii~=nTime
            szNum = [szNum '; '];
        end
    end
else
    szAll = unique([size(data1,1),size(data2,1)],'stable');
    if length(szAll)==1
        szNum = num2str(szAll);
    else
        szNum = [num2str(szAll(1)) ',' num2str(szAll(2))];
    end
end

% set sample type
if nargin<4 || isempty(nUnits)
    nUnits = 'biologically independent samples';
end

% set sample size
sampleSize = ['n=' szNum ' ' nUnits];

% set effect size
effectSize = ['eta=' num2str(statsAll.eta2(1),2)];

% set confidence interval
CI = statsAll.eta2Ci(1,:);
confidenceInterval = ['CI=' num2str(CI(1),2) ',' num2str(CI(2),2)];

% set degrees of freedom
degreesOfFreedom = ['df=' num2str(tblAll{3,3},2)];

% generate output array
anovaStats = {pValue,testName,testStatistic,sampleSize,effectSize,...
    confidenceInterval,degreesOfFreedom};

end

