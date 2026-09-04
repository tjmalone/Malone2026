%% validationPS_forFigs
% Plots histology validation for use in figures


%% Plot data

clear; close all; clc

p1 = 'D:\AD_Project\imagingData\analysis_ManualSelection\HistologyValidation';
cd(p1)

% load data
load('allData.mat','allData')
useData = allData.area;
lens = cellfun(@length,useData);

% set use FOVs
useFOV = 1:8;
% useFOV = [3 4 6 7 8];
nFOV = length(useFOV);

svFolder = 'D:\AD_Project\Revision\PyrSteAnalysis\Results\Validation\';


%% Calculate threshold splits

figure; tiledlayout(1,length(useFOV))
for ff = 1:nFOV
    nexttile(); hold on
    curFOV = useFOV(ff);
    curData = cat(1,useData{curFOV,:});

    ksdensity(curData,0:1:300,'Bandwidth',10);
    ksdensity(useData{curFOV,1},'Bandwidth',10);
    ksdensity(useData{curFOV,2},'Bandwidth',10);

    xline(151)
end


%% Calculate accuracy

% define settings
useFOV = [3 4 6 7 8];
indThresh = [146,161,112,165,151];
nFOV = length(useFOV);
globThresh = 151;

% define ks plot parameters
width = 8;
pts = 0:1:400;
thTypes = {'ind','global'};
cTypes = {'pyr','ste'};
nTh = length(thTypes);
nC = length(cTypes);

% multiplier to deal with comparison direction
cMult = [1 -1];

figure; tiledlayout(1,nTh)

% calculate all distances
dists = 0:100;
yLims = [70 100];
plotData = struct();

for tt = 1:nTh
    % intitialize error rates
    correctRates = zeros(2,length(dists));

    % loop through true cell types
    for cc = 1:nC
        % apply size threshold
        curCell = cell(nFOV,2);
        for ff = 1:nFOV
            curData = useData{useFOV(ff),cc};

            % loop through test cell types
            for mm = 1:2
                cM = cMult(mm);
                if strcmp(thTypes{tt},'ind')
                    curCell{ff,mm} = cM*curData < cM*indThresh(ff)-dists;
                elseif strcmp(thTypes{tt},'global')
                    curCell{ff,mm} = cM*curData < cM*globThresh-dists;
                else
                    error('Invalid threshold type')
                end
            end
        end

        curCat = cell(1,2);
        curNum = cell(1,2);
        for mm = 1:2
            % concatenate identities
            curCat{mm} = cat(1,curCell{:,mm});

            % calculate total cell numbers
            curNum{mm} = sum(curCat{mm},1);
        end

        % calculate true positive rate
        correctRates(cc,:) = curNum{cc}./(curNum{1}+curNum{2})*100;

        % store cell type info
        plotData.(thTypes{tt}).(cTypes{cc}).cell = curCell;
        plotData.(thTypes{tt}).(cTypes{cc}).cat = curCat;
        plotData.(thTypes{tt}).(cTypes{cc}).num = curNum;
    end


    %% Plot error rates

    % plot error rates
    nexttile(tt)
    plot(correctRates')

    % set labels
    title(['True Positive Rates (' thTypes{tt} ')'])
    legend(cTypes,'Location','best')
    xlabel('distance threshold')
    ylabel('True Positive Rate')
    ylim(yLims)
end

savefig([svFolder 'identAccuracy'])

