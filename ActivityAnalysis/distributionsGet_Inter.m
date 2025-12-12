%% distributionsGet_Inter.m
% Calculate consistency distributions as a function of track location for
% learning phase. Specifically, caluclates inter-day consistency and
% prepares intra-day consistency based on calculations in
% distributionsGet_master.
%


%% Set input parameters

clear; clc; close all

p1 = 'D:\AD_Project\imagingData';
cd(p1)

% load learning alignments
load('foldersLearning.mat','foldersLearning')

% load recall alignments
load('foldersRecall.mat','foldersRecall')

% define learning type parameters
LTypes = {'learning','recall'};
nFOV = [length(foldersLearning) length(foldersRecall)];
nDays = [length(foldersLearning{1}) length(foldersRecall{1})];

% set rolling average parameter
rollN = 5;
nBins = 120;
rollBinN = nBins-(rollN-1);

% load inter-day learning data
load('D:\AD_Project\imagingData\data\interDayData_dfof.mat','dataAllInter')
interDataLearn = dataAllInter(1).meanDfof;

% load inter-day recall data
load('data\activityMatGlobalRecall_dfof.mat','activityMatGlobal')
interDataRecall = cellfun(@(x) permute(cat(3,x{:}),[3 2 1]),activityMatGlobal,'UniformOutput',false);


%% Calculate distributions

distActivityInter = struct();

for LT = 1:2
    nFOVCur = nFOV(LT);
    nDaysCur = nDays(LT);

    distActivityInterCur = cell(nFOVCur,nDaysCur);

    % loop through FOV
    for ff = 1:nFOVCur
        disp(ff)

        % concatenate FOV data
        if LT==1
            curInterSource = interDataLearn(ff,:);
            curInterCat = cat(1,curInterSource{:});
        else
            curInterCat = interDataRecall{ff};
        end
        nCells = size(curInterCat,3);

        % initialize output array
        curInterCorr = zeros(rollBinN,nDaysCur-1,nCells);

        % caculate inter-day consistency by cell
        for cc = 1:nCells
            cellInter = curInterCat(:,:,cc);

            % calculate correlation for rolling average
            for m=1:size(cellInter,2)-(rollN-1)
                curRange = m:m+rollN-1;

                % lap to lap correlations
                c = corr(cellInter(:,curRange)',cellInter(:,curRange)');

                % adjacent day correlation
                corrInter = diag(c,-1);
                curInterCorr(m,:,cc) = corrInter;
            end
        end

        for dd = 1:nDaysCur-1
            distActivityInterCur{ff,dd} = squeeze(curInterCorr(:,dd,:))';
        end
        distActivityInterCur{ff,end} = nan(nCells,rollBinN);
    end

    distActivityInter.(LTypes{LT}) = distActivityInterCur;
end

save('data\distActivityInter_dfof.mat','distActivityInter')


%% Prep Intra-day consistency

% load distribution data
load('data\distActivity.mat','distActivity')

% select intra-day consistency
distActivityIntra = distActivity.learn.allRuns.oRBR.all;

% save intra-day consistency
save('data\distActivityIntra.mat','distActivityIntra')

