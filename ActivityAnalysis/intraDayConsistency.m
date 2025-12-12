%% intraDayConsistency.m
% Calculate and plot RBR consistency across days.

clear; clc; close all;


%% Set input parameters

p1 = '';
cd(p1)

% load alignmeft data
load('foldersLearning.mat')

useDays = 1:11;
nDays = length(useDays);
nFOV = length(foldersLearning);

% define region names
regionTypes = {'full','start','stop'};
nRegions = length(regionTypes);

% define reward zones
rewLocs = [240 290;510 560];
useRew = [1 2*ones(1,10)];
trackLength = 600;
binWidth = 5;
maxLaps = 20;


%% Find common cell activity for all mice

dataActIntra = struct();

for rr = 1:nRegions
    dataActIntra.(regionTypes{rr}).RBR_sigMean = cell(nFOV,nDays);
    dataActIntra.(regionTypes{rr}).RBR_sigNext = cell(nFOV,nDays);
    dataActIntra.(regionTypes{rr}).RBR_sigOthers = cell(nFOV,nDays);

    % cycle through FOV
    for ii = 1:nFOV
        disp(ii)

        % cycle through days
        for jj = 1:nDays
            cd(foldersLearning{ii}{jj})
            
            % get lap numbers
            load('imageLapIdx.mat')
            load('md.mat')
            runIdx = getLapIdx(md);
            if size(runIdx,1)>21
                disp(foldersLearning{ii}{jj})
            end

            % define image index range
            lastLapRun = find(imageLapIdx.lapRun<=maxLaps,1,'last');
            lastLapIdx = imageLapIdx.lapIdx(lastLapRun);
            endIdx = runIdx(lastLapIdx,2);

            % get dfofaveragesmooth_interp
            % d1 = dir('dfof_sig*');
            % load(d1(1).name,'dfof_sig')
            % load('abf.mat','abf')
            % dfofaveragesmooth_interp = get_dfofaveragesmoothSubset(dfof_sig,abf,endIdx,trackLength,binWidth);
            load('dfofaveragesmooth_sig_interp.mat')
            dfofaveragesmooth_interp = dfofaveragesmooth_sig_interp;

            % define use bins
            if rr==1
                binStart = 1;
                binEnd = size(dfofaveragesmooth_interp,1);
            elseif rr==2
                binStart = 1;
                binEnd = rewLocs(useRew(jj),1)/binWidth;
            elseif rr==3
                binStart = rewLocs(useRew(jj),2)/binWidth;
                binEnd = size(dfofaveragesmooth_interp,1);
            end

            load('RunByRun_sig\dfofMInterpM_sig.mat')
            dfofMInterpMUse = cellfun(@(x) x(1:lastLapRun,:),dfofMInterpM_sig,'UniformOutput',0);
            % dfofMInterpMUse = dfofMInterpM_sig;

            % try
            % catch
            %     imageLapIdx = findImageLaps([],1);
            %     lastLapRun = find(imageLapIdx.lapRun<=maxLaps,1,'last');
            %     dfofMInterpMUse = cellfun(@(x) x(1:lastLapRun,:),dfofMInterpM_sig,'UniformOutput',0);
            % end

            % get correlation info
            corrInfo = dfofMCorrelationSubset(dfofMInterpMUse,...
                dfofaveragesmooth_interp,binStart,binEnd);

            % extract RBR consistency
            meantoMean = corrInfo.meantoMean(:)';
            meantoNext = corrInfo.meantoNext(:)';
            meantoOthers = corrInfo.meantoOthers(:)';

            % identify true day index
            trueIdx = trueDays(ii,jj);
            if trueIdx>jj
                % set skipped day data to nan
                if trueDays(ii,jj-1)==trueIdx-2
                    dataActIntra.(regionTypes{rr}).RBR_sigMean{ii,jj} = nan(size(meantoMean));
                    dataActIntra.(regionTypes{rr}).RBR_sigNext{ii,jj} = nan(size(meantoNext));
                    dataActIntra.(regionTypes{rr}).RBR_sigOthers{ii,jj} = nan(size(meantoOthers));
                end

                % skip days past day limits
                if trueIdx>max(useDays)
                    continue
                end
            end

            % save data
            dataActIntra.(regionTypes{rr}).RBR_sigMean{ii,trueIdx} = meantoMean;
            dataActIntra.(regionTypes{rr}).RBR_sigNext{ii,trueIdx} = meantoNext;
            dataActIntra.(regionTypes{rr}).RBR_sigOthers{ii,trueIdx} = meantoOthers;
        end
    end
end

cd(p1)
save('data\intraDayData.mat','dataActIntra')


%% Initialize figure parameters

clear; close all; clc

p1 = 'D:\AD_Project\imagingData';
cd(p1)

% load data
load('data\intraDayData.mat','dataActIntra')
load('foldersLearning.mat')
load('data\globalCellTypes.mat','globalCellTypeIdx')

useDays = 1:11;
nDays = length(useDays);
nFOV = length(foldersLearning);

% load genotypes {[WT],[AD]}
load('groupIDs.mat')
nGroups = length(groups);

% manually set genotypes
% groups = {[1:6 9:14 23:24],[7:8 15:22 25:28]};

% define legend
legs = groupIDs;

% define plot colors
colors = {[0 0 1],[1 0 0]};

% set save folder name
svFile = [p1 '\Figures\Figures' char(datetime('today','format','yyyyMMdd')) '\consistency'];
mkdir(svFile);
mkdir([svFile '\data']);


%% Plot RBR consistency

fields = {'RBR_sigMean','RBR_sigNext','RBR_sigOthers','RBR_sigNext','RBR_sigNext'};
regionTypes = {'full','full','full','start','stop'};
ylabs = 'RBR consistency';
ttls = {'RBR consistency: to mean','RBR consistency: to next','RBR consistency: to others',...
    'RBR consistency: pre-reward','RBR consistency: post-reward'};
ylimits = {[0.5 0.725],[0.3 0.6],[0.25 0.55],[0.25 0.55],[0.25 0.55]};

% fields = {'RBR_sigOthers','RBR_sigOthers'};
% regionTypes = {'start','stop'};
% ylabs = 'RBR consistency';
% ttls = {'RBR consistency: pre-reward','RBR consistency: post-reward'};
% ylimits = {[0.25 0.55],[0.25 0.55]};

% set plotting paramters
X = ['FE' string(1:10)];
dayCats = {1,2:11};
sfx ={'cell','fov'};

% close all figures
clAll = 0;

% loop through all fields
for ff = 1:length(fields)
    
    % define current data
    curData = dataActIntra.(regionTypes{ff}).(fields{ff});
    for ii = 1:nFOV
        for jj = 1:nDays
            % identify true day index
            trueIdx = trueDays(ii,jj);
            if trueIdx>jj
                % set skipped day data to nan
                if trueDays(ii,jj-1)==trueIdx-2
                    curData{ii,jj} = nan(size(curData{ii,jj}(alignsLearning{ii}(:,jj))));
                end

                % skip days past day limits
                if trueIdx>max(useDays)
                    continue
                end
            end

            % select cell type (currently hard set to common cells)
            curData{ii,trueIdx} = curData{ii,trueIdx}(alignsLearning{ii}(:,jj));
            % curData{ii,jj} = curData{ii,jj};
        end
    end
    

    %% Plot data

    for mm = 1:2
        figure; hold on

        % separate data by genotype
        useDataCat = cell(nGroups,nDays);

        for ii = 1:length(groups)
            useData = curData(groups{ii},:);

            if mm==1
                % average by cell
                for jj = 1:nDays
                    useDataCat{ii,jj} = cat(2,useData{:,jj})';
                end

            elseif mm==2
                % average by FOV
                useDataFOV = cellfun(@(x) mean(x,'omitnan'),useData);
                for jj = 1:nDays
                    useDataCat{ii,jj} = useDataFOV(:,jj);
                end
            end
        end

        % plot data
        [pAnova] = errorSig(X,useDataCat,colors,legs,dayCats);

        % data1 = cat(2,useDataCat{1,:});
        % data2 = cat(2,useDataCat{2,:});
        % 
        % mean(data1,'all','omitnan')
        % mean(data2,'all','omitnan')

        % set labals
        xlabel('Session')
        ylabel(ylabs)
        title([ttls{ff} ' (by ' sfx{mm} ')'])
        set(gca,'FontSize',12)

        % set axes
        xlim([useDays(1)-0.5 useDays(end)+0.5])
        ylim(ylimits{ff})
        axis('square')

        % save figure
        savefig([svFile '\' regionTypes{ff} '-' fields{ff} '_' sfx{mm} '.fig'])

        % save data
        save([svFile '\data\' regionTypes{ff} '-' fields{ff} '_' sfx{mm} '.mat'],'X','useDataCat')

        % close figure
        if clAll
            close
        end
    end
end
