%% validationPS_forFigs
% Plots histology validation for use in figures


%% Plot data

clear; close all; clc

p1 = 'D:\AD_Project\imagingData\analysis_ManualSelection\HistologyValidation';
cd(p1)

% load data
load('valDataFinal.mat','valDataFinal')
lens = cellfun(@length,valDataFinal);
nData = length(valDataFinal);

% define ks plot parameters
width = 8;
pts = 0:1:400;
types = {'pyramidal','stellate','combined'};
legs = types;

% define separation parameters
center = 151;
cRange = 20;


%% Calculate accuracy

% calculate all distances
dists = 0:100;

% intitialize error rates
correctRates = zeros(2,length(dists));

% calculate correct rates
for ii = 1:2
    isCal = valDataFinal{ii} < center-dists;
    isReel = valDataFinal{ii} > center+dists;

    numCal = sum(isCal,1);
    numReel = sum(isReel,1);

    if ii==1
        correctRates(ii,:) = numCal./(numCal+numReel)*100;
    elseif ii==2
        correctRates(ii,:) = numReel./(numReel+numCal)*100;
    end

    % add correct rates to legend
    legs{ii} = [legs{ii} ' (' num2str(correctRates(ii,cRange),2) '% correct)'];
end




%% Plot size curves

figure; hold on

% plot individaual ks density
for ii = 1:nData
    [k,xi] = ksdensity(valDataFinal{ii},pts,'width',width);
    plot(xi,k*lens(ii))
end


% plot x lines
xline(center)
xline(center+cRange)
xline(center-cRange)

% set labels
legend(legs)
ylabel('Cell count density')
xlabel('Cell area (\mum^2)')

sgtitle('PS Validation')
savefig('valPS_forFig.fig')



