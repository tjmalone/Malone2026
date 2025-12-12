%% distributionsGet_master.m
% Calculate neuronal activity distributions as a function of track
% location fo all learning types (learn, recall), run types (all runs, true
% succes/fail, conditional success/fail), and distribution types (dfof,
% field, to mean RBR, to next RBR, to others RBR). Distribution suffix type
% (dfof vs. sig) can currently be set, but eill cause issues if it is set
% to dfof.
%


%% Set input parameters

clear; clc; close all

p1 = 'D:\AD_Project\imagingData';
cd(p1)

% load learning alignments
load('foldersLearning.mat','foldersLearning','trueDays')

% load recall alignments
load('foldersRecall.mat','foldersRecall')

% set distribution suffix type (_dfof, _sig)
suffix = '_sig';

% set maximum lap number
maxLaps = 20;

% set RBR consistency parameters
rbrStart = 1;
% rbrEnd = 50;
rbrEnd = 102;


%% Calculate distributions

% loop key:
% learning type (lt), run type (rp), distribution type, (dt), success/fail
% (sf), FOV (ff), day (dd)

learningTypes = {'learn','recall'};
runTypes = {'allRuns','trueSF','condSF'};
distTypes = {'dfof','field','fieldAmp','mRBR','nRBR','oRBR','oEMD'};
sfTypes = {'Success','Fail'};

% initialize data struct
distActivity = struct();

% loop through learning types
for lt = learningTypes
    % set current folder and alignment info
    if strcmp(lt,'learn')
        curFolders = foldersLearning;
    elseif strcmp(lt,'recall')
        curFolders = foldersRecall;
    else
        error('Fix learing type list')
    end

    % set FOV and day number
    nFOV = length(curFolders);
    nDays = length(curFolders{1});

    % loop through run types
    for rt = runTypes
        % loop through distribution types
        for dt = distTypes

            % define field identification folder
            if strcmp(rt,'allRuns')
                fieldSubFolds = {'all'};
            else
                fieldSubFolds = sfTypes;
            end

            % loop through success/fail types
            for sf = fieldSubFolds
                distActivity.(lt{:}).(rt{:}).(dt{:}).(sf{:}) = cell(nFOV,nDays);
            end
        end
    end
end


%% Calculate distributions

% loop through learning types
for lt = learningTypes
    % set current folder and alignment info
    if strcmp(lt,'learn')
        curFolders = foldersLearning;
    elseif strcmp(lt,'recall')
        curFolders = foldersRecall;
    else
        error('Fix learing type list')
    end

    % set FOV and day number
    nFOV = length(curFolders);
    nDays = length(curFolders{1});

    % loop through FOV
    for ff = 1:nFOV
        % loop through days
        for dd = 1:nDays
            disp([num2str(ff) '-' num2str(dd)])

            % correct for nan days
            if strcmp(lt,'learn')
                % identify true day index
                trueDay = trueDays(ff,dd);

                % skip days past day limits
                if trueDay>nDays
                    continue
                end
            else
                trueDay = dd;
            end

            % move to current folder
            cd(curFolders{ff}{dd})

            % load interpolated dfof per lap
            if strcmp(suffix,'_dfof')
                load('RunByRun_dfof\dfofMInterpM.mat','dfofMInterpM')
                dfofM = dfofMInterpM;
            else
                load('RunByRun_sig\dfofMInterpM_sig.mat','dfofMInterpM_sig')
                dfofM = dfofMInterpM_sig;
            end


            %% loop through run types
            for rt = runTypes
                
                % define field identification folder
                if strcmp(rt,'allRuns')
                    fieldFold = ['gridAnalysis' suffix];
                    fieldSubFolds = {'all'};
                else
                    if strcmp(rt,'trueSF')
                        fieldFold = 'successFailAnalysis\';
                    else
                        fieldFold = 'successFailAnalysis\post';
                    end
                    fieldSubFolds = sfTypes;
                end
                

                %% loop through success/fail types
                for sf = fieldSubFolds
                    % find last lap index
                    load('imageLapIdx.mat','imageLapIdx')
                    lastLapIdx = find(imageLapIdx.lapRun<=maxLaps,1,'last');

                    if strcmp(rt,'allRuns')
                        % define full field identification folder
                        fieldFoldFull = fieldFold;

                        % set use laps
                        useRuns = 1:lastLapIdx;
                    else
                        % define full field identification folder
                        fieldFoldFull = [fieldFold sf{:}];

                        % set use laps
                        load([fieldFoldFull '\useRuns.mat'],'useRunIdxAll')
                        useRuns = useRunIdxAll;
                        useRuns(useRuns>lastLapIdx) = [];
                    end

                    if isempty(useRuns)
                        continue
                    end

                    % load field information
                    if isfile([fieldFoldFull '\allCellsSub.mat'])
                        % load field locations
                        load([fieldFoldFull '\allCellsSub.mat'])
                    else
                        % find field locations
                        load([fieldFoldFull '\allCellsCorrected.mat'])

                        % save field locations
                        allCellsSub = struct();
                        allCellsSub.dfof = allCellsCorrected.dfofaveragesmooth';
                        allCellsSub.fieldsAmp = allCellsCorrected.dfofaveragesmoothFields';
                        allCellsSub.fields = allCellsSub.fieldsAmp;

                        % threshold fields
                        allCellsSub.fields(allCellsSub.fields~=0) = 1;

                        % interpolate dfof
                        allCellsSub.dfofInterp = zeros(size(allCellsSub.dfof));
                        for n=1:size(allCellsSub.dfof,1)
                            allCellsSub.dfofInterp(n,:) = naninterp(allCellsSub.dfof(n,:));
                        end

                        save([fieldFoldFull '\allCellsSub.mat'],'allCellsSub')
                    end

                    % store dfof and field distribution
                    distActivity.(lt{:}).(rt{:}).dfof.(sf{:}){ff,trueDay}...
                        = allCellsSub.dfofInterp;
                    distActivity.(lt{:}).(rt{:}).field.(sf{:}){ff,trueDay}...
                        = logical(allCellsSub.fields);
                    distActivity.(lt{:}).(rt{:}).fieldAmp.(sf{:}){ff,trueDay}...
                        = allCellsSub.fieldsAmp;

                    % skip run-by-run consistency for single runs
                    if length(useRuns)==1; continue; end

                    % calculate run-by-run consistency
                    if isfile([fieldFoldFull '\allCellsRBR.mat'])
                        load([fieldFoldFull '\allCellsRBR.mat'])
                    else
                        dfofMCur = cellfun(@(x) x(useRuns,:),dfofM,'UniformOutput',false);
                        allCellsRBR = dfofM_RBRbyLoc_Subset(dfofMCur,allCellsSub.dfofInterp');
                        save([fieldFoldFull '\allCellsRBR.mat'],'allCellsRBR')
                    end

                    % calculate overall run-by-run consistency
                    if ~isfile([fieldFoldFull '\allCellsOvrRBR.mat'])
                        dfofMCur = cellfun(@(x) x(useRuns,:),dfofM,'UniformOutput',false);
                        allCellsOvrRBR = dfofMCorrelationSubset(dfofMCur,allCellsSub.dfofInterp');
                        allCellsSubRBR = dfofMCorrelationSubset(...
                            dfofMCur,allCellsSub.dfofInterp',rbrStart,rbrEnd);
                        save([fieldFoldFull '\allCellsOvrRBR.mat'],...
                            'allCellsOvrRBR','allCellsSubRBR','rbrStart','rbrEnd')
                    end

                    % store dfof and field distribution
                    distActivity.(lt{:}).(rt{:}).mRBR.(sf{:}){ff,trueDay}...
                        = allCellsRBR.toMeanMean;
                    distActivity.(lt{:}).(rt{:}).nRBR.(sf{:}){ff,trueDay}...
                        = allCellsRBR.toNextMean;
                    distActivity.(lt{:}).(rt{:}).oRBR.(sf{:}){ff,trueDay}...
                        = allCellsRBR.toOthersMean;
                    distActivity.(lt{:}).(rt{:}).oEMD.(sf{:}){ff,trueDay}...
                        = allCellsRBR.EMDtoOthersMean;

                end
            end
        end
    end
end

cd(p1)


%% Save distribution data


save('data\distActivity.mat','distActivity','-v7.3')


