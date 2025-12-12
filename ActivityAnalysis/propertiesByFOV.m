%% propertiesByFOV.m
% Calculates/collects properties of neural acitivity on a per FOV basis
% for analysis elsewhere. Currently only collects speed score.

clear; clc; close all;


%% Set input parameters

p1 = 'D:\AD_Project\imagingData';
cd(p1)

% load alignment data
load('foldersLearning.mat')

useDays = 1:11;
nDays = length(useDays);
nFOV = length(foldersLearning);
refDay1 = 2;
refDay2 = 11;

subTypes = {'full','even','odd','shuffle'};
nType = length(subTypes);


%% Collect properties for all mice

types = {'thetaEdge','thetaSub'};
dataDiss.thetaEdge = zeros(nDays,nDays,nFOV);
dataDiss.thetaSub = zeros(nDays,nDays,nFOV);

% cycle through FOV
for ii = 1:nFOV
    disp(ii)

    % loop through analysis types
    for fType = 1:2
        %% Collect data for current FOV/analysis type
         
        curAll = struct();
        for st = 1:nType
            curAll.(subTypes{st}) = cell(1,nDays);
        end

        % cycle through days
        for jj = 1:nDays
            % identify true day index
            trueDay = find(trueDays(ii,:)==jj);

            % skip days past day limits
            if isempty(trueDay)
                continue
            end

            cd(foldersLearning{ii}{trueDay})

            % load data
            if fType==1
                load('populationGeometry\sig\thetaEdge.mat')
                curData = thetaEdge;
            else
                load('populationGeometry\sig\thetaSub.mat')
                curData = thetaSub;
            end

            for st = 1:nType
                curAll.(subTypes{st}){jj} = curData.(subTypes{st});
            end

        end

        cd(p1)

        %% Process current FOV/analysis type to get normalized dissimilarity matrix
        
        % calculate initial dissimilarity matrix for all day combinations
        dissInit = zeros(nDays,nDays);
        dissShuff = zeros(nDays,nDays);
        for jj = 1:nDays
            for kk = 1:nDays
                % define analysis pair
                if jj==kk
                    data1 = curAll.odd{jj};
                    data2 = curAll.even{jj};
                else
                    data1 = curAll.full{jj};
                    data2 = curAll.full{kk};
                end

                % calculate dissimilarity value
                if isempty(data1) || isempty(data2)
                    curDiss = nan;
                else
                curDiss = sum((data1-data2).^2,'all','omitnan')^0.5;
                end

                % store dissimiliarity value
                dissInit(jj,kk) = curDiss;

                % calculate shuffle dissimilarity
                data1 = curAll.shuffle{jj};
                data2 = curAll.shuffle{kk};
                if (isempty(data1) || isempty(data2)) || jj==kk
                    curShuff = nan;
                else
                    curShuff = sum((data1-data2).^2,'all','omitnan')^0.5;
                end
                dissShuff(jj,kk) = curShuff;
            end
        end

        % calculate within-day average
        dissWithinMean = mean(diag(dissInit),'omitnan');

        % calculate shuffle average
        dissShuffMean = mean(dissShuff,'all','omitnan');

        % calculate normalized dissimilarity matrix
        dissNorm = (dissInit-dissWithinMean)/(dissShuffMean-dissWithinMean);
        
        % store normalized matrix
        dataDiss.(types{fType})(:,:,ii) = dissNorm;
    end
end

% save data
cd(p1)
save('data\dataByFOV.mat','dataDiss')


%% Plot data

clear; close all; clc

load('groupIDs.mat');
load('data\dataByFOV.mat','dataDiss')
fields = fieldnames(dataDiss);
useDays = 2:11;

for fType = 1:2

    figure
    sgtitle(fields{fType})
    pixMax = nan;
    pixMin = nan;

    for gg = 1:2
        subplot(1,2,gg); hold on

        curPlot = mean(dataDiss.(fields{fType})(:,:,groups{gg}),3,'omitnan');
        curPlot = curPlot(useDays,useDays);
        imagesc(curPlot)

        title(groupIDs{gg})

        pixMax = max(pixMax,max(curPlot,[],'all'));
        pixMin = min(pixMin,min(curPlot,[],'all'));
    end

    for gg = 1:2
        subplot(1,2,gg); hold on

        clim([pixMin pixMax])
        colorbar
        ylim([0 length(useDays)]+0.5)
        xlim([0 length(useDays)]+0.5)
        axis('square')
    end


end


