%% BAmodel_postProcess
% Process data using a final model determined elsewhere. This
% generic version processes both training and test data, using lme of
% glm, and using LOOCV or LTOCV. BAmodel_plotGen should be run after to
% plot results.
%


%% Perform genotype classification

clear; close all; clc

p1 = 'D:\AD_Project\Behavior';
cd(p1)

% load data
load('data/model_vB/LOOCVpredictGeno_noDS_vB.mat','LOOData')

% load genotype and sex information
load('groupIDs.mat')
useSex = 1:3;
nSex = length(useSex);
nGeno = length(groups);
nDays = 10;

colors3 = {[0 0 0];[255 102 102]/255;[0 0 175]/255};


%%

% define cell types
cellTypes = fieldnames(LOOData);
nCellTypes = length(cellTypes);

% initialize figure
figure; tiledlayout(1,nCellTypes)
sgtitle('WT vs. PS19 Classification: ')

for cc = 1:nCellTypes
    % get current data
    LOODataCur = LOOData.(cellTypes{cc});

    % define number of samples
    % initialize classification output
    nShuff = size(LOODataCur.yTest,2);
    perCorrect = zeros(nShuff,nSex);

    % define use data
    dataTrue = LOODataCur.yTest;
    dataPred = LOODataCur.yPredTest;
    useMice = LOODataCur.idxMouse;
    nMiceAll = length(useMice);

    % reconstruct mouse table
    tblM = LOODataCur.tblM;
    tblMFull = zeros(nMiceAll,nDays);
    curCol = 1;
    tblMFull(ismember(useMice,tblM(1)),curCol) = 1;
    for ii = 2:length(tblM)
        if tblM(ii)<tblM(ii-1)
            curCol = curCol+1;
        end
        tblMFull(ismember(useMice,tblM(ii)),curCol) = 1;
    end

    % define use days as the last N days. This is only valid up to 4 because of missing days
    useDaysN = 4;

    % construct use days cell array
    daysArray = cell(nMiceAll,1);
    for mm = 1:nMiceAll
        % find first day in range with data
        curMRow = tblMFull(mm,:);
        dayStart = find(curMRow(end-useDaysN+1:end)==1,1)+nDays-useDaysN;

        % skip empty days
        if ~isempty(dayStart)
            daysArray{mm} = sum(curMRow(1:dayStart)):sum(curMRow);
        end
    end
    daysArrayRep = repmat(daysArray,1,nShuff);

    meanTrue = cellfun(@(x,y) mean(x(y),'omitnan'),dataTrue,daysArrayRep);
    meanPred = cellfun(@(x,y) mean(x(y),'omitnan'),dataPred,daysArrayRep);

    for ss = 1:nSex
        % define current mice
        curMice = intersect(useMice,sexes{ss});
        nMice = length(curMice);

        % initialize classification array (1=WT, 2=AD)
        catTrue = zeros(nMice,nShuff);
        catPred = zeros(nMice,nShuff);

        for mm = 1:nMice
            curMouse = curMice(mm);
            if ismember(curMouse,groups{1})
                catTrue(mm,:) = 1;
            elseif ismember(curMouse,groups{2})
                catTrue(mm,:) = 2;
            end

            % calculate current mouse excluded mean of true behavior
            behMeanExc = zeros(2,nShuff);
            for gg = 1:nGeno
                % define current genotype mice (excluding current mouse)
                excMice = intersect(curMice,groups{gg});

                % exclude current mouse except for pyramidal 
                if cc~=5
                    excMice(excMice==curMouse) = [];
                end

                % calculate current mouse excluded mean
                behMeanExc(gg,:) = mean(meanTrue(ismember(useMice,excMice),:),1,'omitnan');
            end

            % get current mouse prediction
            curPred = meanPred(ismember(useMice,curMouse),:);

            % classify mouse
            predDif = abs(behMeanExc-curPred);
            [~,idx] = min(predDif,[],1);
            catPred(mm,:) = idx;

            % remove ties and nan (code as -1 and -2)
            catPred(mm,predDif(1,:)==predDif(2,:)) = -1;
            catPred(mm,all(isnan(predDif),1)) = -2;
        end

        % find correct classification prediction
        nanPred = catPred<0;
        curCorrect = nan(size(catTrue));
        curCorrect(~nanPred) = catTrue(~nanPred)==catPred(~nanPred);
        curCorrect(catPred==-1) = 0.5;
        perCorrect(:,ss) = mean(curCorrect,1,'omitnan');
    end


    %% Plot classification accuracy versus shuffle

    nexttile(cc); hold on
    ttlSfx = ': prctile = ';

    h = zeros(1,nSex);
    for ss = 1:nSex
        % get current shuffle and true data
        plotShuff = perCorrect(1:end-1,ss);
        plotTrue = perCorrect(end,ss);

        % calculate percentile of shuffle
        curPrc = invPercentile(plotTrue,plotShuff);
        if ss>1
            ttlSfx = [ttlSfx ', '];
        end
        ttlSfx = [ttlSfx num2str(curPrc,'%.3g')];

        % plot shuffle distribution
        [yShuff,xShuff] = ksdensity(plotShuff,0:0.01:1);
        plot(xShuff,yShuff,'Color',colors3{ss});

        % plot true value
        h(ss) = xline(plotTrue,'Color',colors3{ss});
    end

    title([cellTypes{cc} ttlSfx])
    legend(h,sexIDs,'Location','best')
    set(gca,'FontSize',12)
end


