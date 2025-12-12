function [pAnova,pMC,pLabel] = anovaO2W_BH(data,mcON)
% Performs an ordinary 2-Way ANOVA assuming sphericity. Calculates
% multiple comparisons using a non-paired ttest with Bonferonni correction.
% Anova p-vaule has not been check against Prism values. Intended for
% situations where the time is not matched across samples.
%
% Inputs:
%       data - cell array of data. Each cell should be a 1D array of values
%       for a given group and time point (2 groups, n times)
%       mcON - whether multiple comparisons p values should be corrected
%       for number of comparisons
%

% calculate condition numbers
nG = size(data,1);
if nG~=2
    error('Only two groups can be compared')
end
nT = size(data,2);

% define conditions
valG = 1:nG;
valT = 1:nT;

% initialize anova arrays
catData = [];
grG = [];
grT = [];

% initialize multiple comparisons
MC = zeros(1,nT);


%% Create anova arrays

for ii = 1:nT
    for jj = 1:nG
        curData = data{jj,ii};
        curLen = length(curData);
        
        catData = [catData; curData];
        grG = [grG; valG(jj)*ones(curLen,1)];
        grT = [grT; valT(ii)*ones(curLen,1)];
    end
    
    % calcualte multiple comparisons p value
    [~,MC(ii)] = ttest2(data{1,ii},data{2,ii});
end

% correct for NaN values
if any(isnan(MC))
    MC(isnan(MC)) = 1;
    disp('NaN p-values set to 1')
end


%% Run ANOVA

[pAnova,~,~] = anovan(catData,{grG,grT},'model','interaction',...
    'varnames',{'geno','time'},'display','off');

pLabel = {'X','T','X:T'};

% run multiple corrections
if nargin<2 || mcON~=0
    [~,pMC] = bonferroni_holm(MC);
else
    pMC = MC;
end

end

