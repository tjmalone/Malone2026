%% propertiesByCell.m
% Calculates/collects properties of neural acitivity on a per cell basis
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

% define parameters
patterns = [1 1; 1 0; 0 0; 0 1];
nPatt = size(patterns,1);
% selfBins = [1:50;121:170];
selfBins = [1:102;121:222];

maxBin = 240;


%% Find properties for all mice

dataAll = struct();
dataAll.speedScore = cell(nFOV,nDays);
dataAll.speedTile = cell(nFOV,nDays);
dataAll.spatialInfo = cell(nFOV,nDays);
dataAll.spatialSelectivity = cell(nFOV,nDays);
dataAll.dfofInField = cell(nFOV,nDays);
dataAll.dfofNonField = cell(nFOV,nDays);
dataAll.FWidthNum = cell(nFOV,nDays);
dataAll.FWidthSum = cell(nFOV,nDays);
dataAll.FWidthMean = cell(nFOV,nDays);
dataAll.FWidthMax = cell(nFOV,nDays);
dataAll.corrNear = cell(nFOV,nDays);
dataAll.corrFar = cell(nFOV,nDays);
dataAll.corrDiff = cell(nFOV,nDays);
dataAll.ratioSS = cell(nFOV,nDays);
dataAll.ratioFF = cell(nFOV,nDays);
dataAll.RBR_S = cell(nFOV,nDays);
dataAll.RBR_F = cell(nFOV,nDays);
dataAll.selfRun_S = cell(nFOV,nDays);
dataAll.selfRun_F = cell(nFOV,nDays);

% cycle through FOV
for ii = 1:nFOV
    disp(ii)

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
        load('speed_dfof_sig\speed.mat')
        nCells = length(speed.scores);

        % process speed scores
        curScores = speed.scores;
        sortShuffs = sort(speed.shuffleScores,1);
        speedTile = invPercentile(curScores,sortShuffs);

        % store speed scores
        dataAll.speedScore{ii,jj} = curScores';
        dataAll.speedTile{ii,jj} = speedTile';

        % load spatial information
        load('spatialInfo_sig\SInew.mat')
        dataAll.spatialInfo{ii,jj} = SI;

        % load spatial selectivity information
        load('spatialInfo_sig\spatialSelectivity.mat','SS','inFields','nonFields')
        dataAll.spatialSelectivity{ii,jj} = SS;
        dataAll.dfofInField{ii,jj} = inFields;
        dataAll.dfofNonField{ii,jj} = nonFields;

        % load field information
        if isfile('gridAnalysis_sig\allCellsWidths.mat')
            load('gridAnalysis_sig\allCellsWidths.mat','fieldWidths')
        else
            load('gridAnalysis_sig\allCellsCorrected.mat','allCellsCorrected')
            fieldWidths = allCellsCorrected.fieldWidths;
            save('gridAnalysis_sig\allCellsWidths.mat','fieldWidths')
        end

        % calculate width information
        widthNum = cellfun(@length,fieldWidths);
        widthSum = cellfun(@sum,fieldWidths);
        widthMean = widthSum./widthNum;
        widthMax = cellfun(@max, fieldWidths);

        % store width information
        dataAll.FWidthNum{ii,jj} = widthNum';
        dataAll.FWidthSum{ii,jj} = widthSum';
        dataAll.FWidthMean{ii,jj} = widthMean';
        dataAll.FWidthMax{ii,jj} = widthMax';


        %% Calculate near-/far- inter-lap correlation

        if isfile('RunByRun_sig\corrInfoLapDistance.mat')
            load('corrInfoLapDistance.mat','corrInfoLapDistance')
        else
            % define last run
            maxLaps = 20;
            load('imageLapIdx.mat','imageLapIdx')
            lastLapRun = find(imageLapIdx.lapRun<=maxLaps,1,'last');

            % load dfofsmooth
            load('dfofaveragesmooth_sig_interp.mat','dfofaveragesmooth_sig_interp')

            % load and process spatial dfof
            load('RunByRun_sig\dfofMInterpM_sig.mat','dfofMInterpM_sig')
            dfofMInterpMUse = cellfun(@(x) x(1:lastLapRun,:),dfofMInterpM_sig,'UniformOutput',0);

            % get correlation info
            corrInfo = dfofMCorrelationSubset(dfofMInterpMUse,dfofaveragesmooth_sig_interp);
            mask = repmat(tril(true(lastLapRun), -1),[1,1,nCells]);
            corrInfoLapDistance = nan(lastLapRun,lastLapRun,nCells);
            corrInfoLapDistance(mask) = corrInfo.toOthers;

            % save correlation info
            save('corrInfoLapDistance.mat','corrInfoLapDistance')
        end

        % define near and far distances; close is set to 3 based on data,
        % far is set to 9 to match number of elements. About one third each
        distNear = 3;
        distFar = 9;

        % define near and far indices
        idxNear = tril(true(lastLapRun),-1) & ~tril(true(lastLapRun),-1-distNear);
        idxFar = tril(true(lastLapRun),-distFar);

        corrNear = zeros(nCells,1);
        corrFar = zeros(nCells,1);
        for cc = 1:nCells
            curCorr = corrInfoLapDistance(:,:,cc);
            corrNear(cc) = mean(curCorr(idxNear),'omitnan');
            corrFar(cc) = mean(curCorr(idxFar),'omitnan');
        end

        dataAll.corrNear{ii,jj} = corrNear;
        dataAll.corrFar{ii,jj} = corrFar;
        dataAll.corrDiff{ii,jj} = corrNear-corrFar;


        %% Calculate success order ratios

        % % load success/fail runs
        % load('successFailAnalysis\Success\useRuns.mat','useRunIdxAll')
        % idxS = useRunIdxAll;
        % load('successFailAnalysis\Fail\useRuns.mat','useRunIdxAll')
        % idxF = useRunIdxAll;
        % 
        % % merge success/fail runs
        % idxAll = sort([idxS;idxF]);
        % runSF = false(size(idxAll));
        % runSF(ismember(idxAll,idxS)) = true;
        % 
        % % find success/fail matches
        % matches = zeros(nPatt,1);
        % for kk = 1:nPatt
        %     matches(kk) = length(strfind(runSF',patterns(kk,:)));
        % end
        % matches = reshape(matches,2,2);
        % 
        % % determine success fail ratio
        % sfRatio = matches(1,:)./sum(matches,1);
        % 
        % dataAll.ratioSS{ii,jj} = sfRatio(1)*ones(nCells,1);
        % dataAll.ratioFF{ii,jj} = sfRatio(2)*ones(nCells,1);


        %% Calculate success/fail RBR subset

        % % load success/fail RBR
        % if isfile('successFailAnalysis\Success\allCellsOvrRBR.mat')
        %     load('successFailAnalysis\Success\allCellsOvrRBR.mat')
        %     RBR_S = allCellsSubRBR.meantoNext';
        % else
        %     RBR_S = nan(nCells,1);
        % end
        % 
        % if isfile('successFailAnalysis\Fail\allCellsOvrRBR.mat')
        %     load('successFailAnalysis\Fail\allCellsOvrRBR.mat')
        %     RBR_F = allCellsSubRBR.meantoNext';
        % else
        %     RBR_F = nan(nCells,1);
        % end
        % 
        % dataAll.RBR_S{ii,jj} = RBR_S;
        % dataAll.RBR_F{ii,jj} = RBR_F;


        %% Calculate joint lap correlation

        % load success/fail indices
        load('successFailAnalysis\Success\useRuns.mat')
        idxS = useRunIdxAll;

        load('successFailAnalysis\Fail\useRuns.mat')
        idxF = useRunIdxAll;

        % load success/fail RBR
        for sf = 1:2
            % define run class
            if sf==1
                subFolder = 'Success';
            else
                subFolder = 'Fail';
            end

            loadFile = ['successFailAnalysis\' subFolder '\jointActivity.mat'];
            if isfile(loadFile)

                % load joint activity
                load(loadFile)

                if ~isempty(jointActivity)
                    %% Calculate dfofM_interp

                    rawDfofM = jointActivity.dfofM;

                    if ~isfield(jointActivity,'dfofMInterpM')
                        dfofMInterpM = cell(size(rawDfofM));
                        for n=1:length(rawDfofM)
                            % define current cell
                            thisdfofM = rawDfofM{n};
                            if isempty(thisdfofM); continue; end
                            dfofMInterpM{n} = zeros(size(thisdfofM));

                            % interpolate cell
                            for m = 1:size(thisdfofM,1)
                                %if there are less than 2 non-nan numbers, keep nan
                                if sum(~isnan(thisdfofM(m,:)))>1
                                    dfofMInterpM{n}(m,:) = naninterp(thisdfofM(m,:));
                                else
                                    dfofMInterpM{n}(m,:) = thisdfofM(m,:);
                                end
                            end
                        end
                        jointActivity.dfofMInterpM = dfofMInterpM;

                        save(loadFile,'jointActivity')
                    end

                    curDfofM = jointActivity.dfofMInterpM;


                    %% Calculate correlation

                    % identify next success/fail subsets
                    load(['successFailAnalysis\' subFolder '\abfJoint.mat'],'useRunIdxAll');
                    idxNext = useRunIdxAll+1;
                    NS = ismember(idxNext,idxS);
                    NF = ismember(idxNext,idxF);

                    % fill missing cell data
                    curDfofM(cellfun(@isempty,curDfofM)) = {nan(length(idxNext),maxBin)};

                    % perform run-by-run correlation for all runs
                    curRunCorrAll = cellfun(@(x) corr(x(:,selfBins(1,:))',x(:,selfBins(2,:))',...
                        'rows','pairwise'),curDfofM,'UniformOutput',false);

                    % take mean of self correlations
                    curRunCorrMean = cellfun(@(x) mean(diag(x),'omitnan'),curRunCorrAll)';

                    % take mean of next success/fail subsets
                    curRunCorrMeanNS = cellfun(@(x) mean(diag(x(NS,NS)),'omitnan'),curRunCorrAll)';
                    curRunCorrMeanNF = cellfun(@(x) mean(diag(x(NF,NF)),'omitnan'),curRunCorrAll)';

                else
                    % fill missing session
                    curRunCorrMean = nan(nCells,1);
                    curRunCorrMeanNS = nan(nCells,1);
                    curRunCorrMeanNF = nan(nCells,1);
                end
            else
                % fill missing session
                curRunCorrMean = nan(nCells,1);
                curRunCorrMeanNS = nan(nCells,1);
                curRunCorrMeanNF = nan(nCells,1);
            end

            if length(curRunCorrMean)~=length(curRunCorrMeanNS)
                error('Mismatch')
            end
            if sf==1
                dataAll.selfRun_S{ii,jj} = curRunCorrMean;
                dataAll.selfRun_SS{ii,jj} = curRunCorrMeanNS;
                dataAll.selfRun_SF{ii,jj} = curRunCorrMeanNF;
            else
                dataAll.selfRun_F{ii,jj} = curRunCorrMean;
                dataAll.selfRun_FS{ii,jj} = curRunCorrMeanNS;
                dataAll.selfRun_FF{ii,jj} = curRunCorrMeanNF;
            end
        end

    end
end


cd(p1)
save('data\dataByCell.mat','dataAll')

