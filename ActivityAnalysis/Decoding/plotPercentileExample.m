% plotPercentileExample

clear; close all;clc

p1 = 'D:\AD_Project\imagingData\analysis_Decoding';
cd(p1)

load('D:\AD_Project\imagingData\foldersLearning.mat')
load('D:\AD_Project\imagingData\groupIDs.mat')
load('D:\AD_Project\imagingData\data\cueIdent\cueIdentityDataCat_common.mat')

dataCue = cueIdentityDataCat.ampDiffMean.all.cue;
dataAnti = cueIdentityDataCat.ampDiffMean.all.anticue;

nFOV = 42;
nDays = 10;
nGeno = 2;

% define example info
FOV = 31;
day = 10;
cellID = 29;

dataMeanCue = cellfun(@(x) mean(x,2,'omitnan'),dataCue,'UniformOutput',false);

idxs = [0 0];
% loop through FOV
for ii = 1:nFOV

    % loop through genotypes
    for jj = 1:nGeno
        %% Calculate amplitude difference
        if ~ismember(ii,groups{jj}); continue; end

        % get current size and indices
        curN = size(alignsLearning{ii},1);
        curIdx = idxs(jj)+1:idxs(jj)+curN;
        idxs(jj) = idxs(jj)+curN;

        %% Calculate percentile of anti-cue

        for kk = 1:nDays
            % select current data
            curDataCue = dataMeanCue{jj,kk}(curIdx);
            curDataAnti = dataAnti{jj,kk}(curIdx,:);

            for dd = 1:size(curDataCue,1)
                if ismember(ii,FOV) && kk==day && ismember(dd,cellID)
                    %% Plot example

                    % calculate percentile
                    sortAnti = sort(curDataAnti(dd,:));
                    sortAnti(isnan(sortAnti)) = [];
                    curLess = sum(sortAnti<curDataCue(dd));
                    curEqual = sum(sortAnti==curDataCue(dd));
                    curTile = (curLess+curEqual/2)/length(sortAnti)*100;

                    % skip full nans
                    if all(isnan(sortAnti))
                        continue
                    end

                    % plot histogram
                    figure; hold on
                    ksdensity(sortAnti,0:0.01:1,'bandwidth',0.05)
                    xline(curDataCue(dd),'b')

                    % plot 50 percentile
                    prc = prctile(sortAnti,50);
                    xline(prc,'k')

                    title(['Percentile = ' num2str(curTile,2)])

                    % save example data
                    exampleCue = curDataCue(dd);
                    exampleAnti = sortAnti;
                    save('examplePercentile.mat','exampleCue','exampleAnti')

                    % ylim([0 2])
                    savefig('examplePercentile.fig')
                end
            end
        end
    end
end

