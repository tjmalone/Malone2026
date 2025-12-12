function [pAnova,pMC] = errorSig(X,data,colors,lgnd,dayCats,indON,anovaType)
%% errorSig
% plots line graph with errorbars between two groups and calculates
% significance values. Uses plotErrorSig for actual plotting.
%

nSess = size(data,2);

if isempty(X)
    X = 1:nSess;
end

if nargin<3 || isempty(colors)
    colors = {[0 0 1];[1 0 0]};
end

if nargin<4 || isempty(lgnd)
    lgnd = [];
end

if nargin<5 || isempty(dayCats)
    dayCats = {1:nSess};
end

if nargin<6 || isempty(indON)
    indON = 0;
end

if nargin<7 || isempty(anovaType)
    anovaType = 'O2W';
end

nCats = length(dayCats);


%% Calculate p values 

pMC = nan(nSess,1);
pAnova = nan(nSess,2);

for ii = 1:nCats
    lenCur = length(dayCats{ii});
    if lenCur>1
        % calculate anova p values with multiple comparisons
        if strcmp(anovaType,'O2W')
            % perform ordinary 2-way anova
            [pACur,pMCCur] = anovaO2W_BH(data(:,dayCats{ii}),1);
        elseif strcmp(anovaType,'RM2W')
            % perform repeated measures 2-way anova
            curData = data(:,dayCats{ii});
            AData = cell(2,1);
            for aa = 1:2
                AData{aa} = cat(2,curData{aa,:});
            end
            [pACur,pMCCur] = anovaRM2W_full_BH(AData{1},AData{2},1);
        else
            error('Invalid ANOVA type')
        end
        pAnova(dayCats{ii},:) = pACur([1 3])'.*ones(lenCur,2);
    else
        % calculate single day p values
        [~,pMCCur] = ttest2(data{1,dayCats{ii}},data{2,dayCats{ii}});
    end
    
    pMC(dayCats{ii}) = pMCCur;
end


%% Plot data

hold on

plotErrorSigCell(1:nSess,data,lgnd,pAnova,pMC,colors,indON)

xticks(1:nSess)
set(gca,'XTickLabels',X)

end

