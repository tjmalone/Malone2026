%%
% Prepares pre-calculated data for use in timecoursePlot_master.
% Specifically, process inter-day consistency, cue amplitude difference,
% and cue template correlations. Some calculations are based on previously
% run scripts, wihle some are directly calculated.
%

%% Correct inter-day consistency

clear
load('data/interDayData.mat','dataAllInter')
load('foldersLearning.mat')
load('groupIDs.mat')

interDay = dataAllInter(1).corrMInd;

nFOV = 42;
nDays = 10;
nGeno = 2;

excludeFOV = find(trueDays(:,end)~=size(trueDays,2));
excludeFOV = [31 32 17 18 20 30];
excludeFOV = [];

dataOut = cell(nFOV,nDays);

idxs = [0 0];
% loop through FOV
for ii = 1:nFOV

    % loop through genotypes
    for jj = 1:nGeno
        if ~ismember(ii,groups{jj}); continue; end

        % get current size and indices
        curN = size(alignsLearning{ii},1);
        curIdx = idxs(jj)+1:idxs(jj)+curN;
        idxs(jj) = idxs(jj)+curN;

        % skip excluded mice
        if ismember(ii,excludeFOV); continue; end

        % extract day-to-day correlation
        curData = interDay{jj}(:,:,curIdx);
        adjDiag = zeros(size(curData,3),size(curData,1)-1);
        for kk = 1:size(curData,3)
            adjDiag(kk,:) = diag(curData(:,:,kk),1);
        end

        for kk = 1:nDays
            dataOut{ii,kk} = adjDiag(:,kk);
        end
    end
end

save('data/interDaySubset.mat','dataOut')


%% Correct amplitude difference

clear

p1 = 'D:/AD_Project/imagingData';
cd(p1)

load('foldersLearning.mat')
load('groupIDs.mat')
load('D:/AD_Project/imagingData/data/cueIdent/cueIdentityDataCatNL_common.mat')

dataCue = cueIdentityDataCat.ampDiffMean.all.cue;
dataAnti = cueIdentityDataCat.ampDiffMean.all.anticue;

nFOV = 42;
nDays = 10;
nGeno = 2;
nAntiPairs = 200;
% cueUsePairs = [1 4 5 7 8 9];
cueUsePairs = 1:15;

dataMeanCue = cellfun(@(x) mean(x(:,cueUsePairs),2,'omitnan'),dataCue,'UniformOutput',false);
dataMeanAnti = cellfun(@(x) mean(x,2,'omitnan'),dataAnti,'UniformOutput',false);

dataOut = struct();
dataOut.cue = cell(nFOV,nDays+1);
dataOut.anti = cell(nFOV,nDays+1);
dataOut.cuePer = cell(nFOV,nDays+1);
dataOut.IOraw = cell(nFOV,nDays+1);
dataOut.IOnorm = cell(nFOV,nDays+1);

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

        for kk = 1:nDays
            % extract day-to-day correlation
            dataOut.cue{ii,kk+1} = dataMeanCue{jj,kk}(curIdx);
            dataOut.anti{ii,kk+1} = dataMeanAnti{jj,kk}(curIdx);
        end

        dataOut.cue{ii,1} = nan(length(curIdx),1);
        dataOut.anti{ii,1} = nan(length(curIdx),1);


        %% Calculate percentile of anti-cue

        for kk = 1:nDays
            % select current data
            curDataCue = dataMeanCue{jj,kk}(curIdx);
            curDataAnti = dataAnti{jj,kk}(curIdx,:);

            % % remove nans
            % curNans = isnan(curDataCue);
            % curCueData = curDataCue(curNans==0);
            % curAntiData = curDataAnti(curNans==0,:);

            tileData = zeros(size(curDataCue,1),1);

            for dd = 1:size(curDataCue,1)
                if isnan(curDataCue(dd))
                    tileData(dd) = NaN;
                else
                    sortAnti = sort(curDataAnti(dd,:));
                    curLess = sum(sortAnti<curDataCue(dd));
                    curEqual = sum(sortAnti==curDataCue(dd));

                    curTile = (curLess+curEqual/2)/sum(~isnan(sortAnti))*100;
                    if curTile>100
                        return
                    end
                    tileData(dd) = curTile;
                end
            end
            % store percentile
            dataOut.cuePer{ii,kk+1} = tileData;

            % calculate IO ratio
            I = curDataCue;
            O = mean(curDataAnti,2,'omitnan');
            dataOut.IOraw{ii,kk+1} = I./O;
            dataOut.IOnorm{ii,kk+1} = (I-O)./(I+O);
        end

        dataOut.cuePer{ii,1} = nan(length(curIdx),1);

    end
end

save('data/cueIdentitySubset.mat','dataOut')



%% Generate cue template correlations

clear; clc;
p1 = '/MATLAB Drive/FY2025/imagingData';
cd(p1)

% load folder and index data
load('foldersLearning.mat')
load('data/cellSelect.mat','cellSelect')
cellSelect = cellSelect.learn;

% define fov and day number
nFOV = size(trueDays,1);
nDays = size(trueDays,2);

% set parameters
binWidth = 5;
cueDil = 2;
% maxBin = 102;
% cueDil = [];
maxBin = [];

% define cue template
rewLoc = {[240 290],[510 560]};

dataActivity = struct();

% loop through fov and day
for ff = 1:nFOV
    for dd = 1:nDays
        % correct for "true" day
        trueInv = find(trueDays(ff,:)==dd);
        if isempty(trueInv)
            continue
        end

        %% Load data and process template

        % change to current directory
        curFolderName = foldersLearning{ff}{trueInv};
        newFolderName = replace(curFolderName,'D:\AD_Project\imagingData\data','/MATLAB Drive/FY2025/SubData');
        newFolderName = replace(newFolderName,'\','/');
        cd(newFolderName)

        % load binned activity and cue template
        d = dir('dfofaveragesmooth_sig*');
        load(d(1).name)
        load('cueAnalysis_sig/tempRL.mat','tempRL')

        %% Loop through all subcategories

        % loop through sexes
        sexIDs = fieldnames(cellSelect)';
        for ss = sexIDs

            % loop through cell types
            cellTypes = fieldnames(cellSelect.(ss{:}))';
            for ct = cellTypes

                % loop through morphology
                morphTypes = fieldnames(cellSelect.(ss{:}).(ct{:}))';
                for mt = morphTypes

                    % loop through genotypes
                    genotypes = fieldnames(cellSelect.(ss{:}).(ct{:}).(mt{:}))';
                    for gg = genotypes
                        % identify current cells
                        curCells = cellSelect.(ss{:}).(ct{:}).(mt{:}).(gg{:}){ff,dd};

                        if isempty(curCells)
                            corrNoLags = [];
                            lagsAbs = [];
                            dfofRdgBkg = [];
                            dfofRdg = [];
                            dfofBkg = [];
                            dfofIO = [];
                            dfofI = [];
                            dfofO = [];
                        else
                            % get current activity data
                            curDfof = mean(dfofaveragesmooth_sig(:,curCells),2,'omitnan');

                            % define evironment variables
                            if trueInv>1
                                rewLocUse = round(rewLoc{2}/5);
                            else
                                rewLocUse = round(rewLoc{1}/5);
                            end

                            dataQuant = distrQuant(curDfof',tempRL,rewLocUse,cueDil,maxBin);

                            corrNoLags = dataQuant.CorrNoLag;
                            lagsAbs = dataQuant.LagToPeakCorr*binWidth;
                            dfofRdgBkg = dataQuant.RdgBkg;
                            dfofRdg = dataQuant.Rdg;
                            dfofBkg = dataQuant.Bkg;

                            I = dataQuant.CueIn;
                            O = dataQuant.CueOut;
                            dfofI = I;
                            dfofO = O;
                            dfofIO = I/O;
                        end
                        
                        % store correlation values
                        dataActivity.corrNoLag.(ct{:}).(ss{:}).(mt{:}).(gg{:}){ff,dd} = corrNoLags;
                        dataActivity.peakLag.(ct{:}).(ss{:}).(mt{:}).(gg{:}){ff,dd} = lagsAbs;
                        dataActivity.dfofRdgBkg.(ct{:}).(ss{:}).(mt{:}).(gg{:}){ff,dd} = dfofRdgBkg;
                        dataActivity.dfofRdg.(ct{:}).(ss{:}).(mt{:}).(gg{:}){ff,dd} = dfofRdg;
                        dataActivity.dfofBkg.(ct{:}).(ss{:}).(mt{:}).(gg{:}){ff,dd} = dfofBkg;
                        dataActivity.dfofIO.(ct{:}).(ss{:}).(mt{:}).(gg{:}){ff,dd} = dfofIO;
                        dataActivity.dfofI.(ct{:}).(ss{:}).(mt{:}).(gg{:}){ff,dd} = dfofI;
                        dataActivity.dfofO.(ct{:}).(ss{:}).(mt{:}).(gg{:}){ff,dd} = dfofO;
                    end
                end
            end
        end

    end
end

cd(p1)
save('data/corrInfoGlobal.mat','dataActivity')


%% Calculate in-cue versus out-cue decoding

% Now in dedicated code: prepareHeatDataNewDecodeIO

% clear; clc;
% p1 = 'D:/AD_Project/imagingData';
% cd(p1)
% 
% % load folder and index data
% load('D:/AD_Project/imagingData/foldersLearning.mat')
% load('data/cellSelect.mat','cellSelect')
% cellSelect = cellSelect.learn;
% 
% % load data
% load('data/decoding/decodingData_grid.mat','decodingData')
% decodeUse.grid = decodingData.meanPerCorrectSpatial;
% 
% load('data/decoding/decodingData_nongrid.mat','decodingData')
% decodeUse.nongrid = decodingData.meanPerCorrectSpatial;
% 
% load('data/decoding/decodingData_all.mat','decodingData')
% decodeUse.allType = decodingData.meanPerCorrectSpatial;
% 
% % define fov and day number
% nFOV = size(trueDays,1);
% nDays = 10;
% 
% % set parameters
% binWidth = 5;
% cueDil = 2;
% maxBin = [];
% cellTypes = {'allType','grid','nongrid'};
% 
% % load cue templates
% cueFolder = 'D:/AnalysisCode/PostAnalysis/Cues/6mEnv2/';
% load([cueFolder 'tempL.mat'])
% load([cueFolder 'tempR.mat'])
% cueTemp = tempL | tempR;
% rewLoc = [510 560]/5;
% 
% decodingActivity = struct();
% 
% % loop through fov and day
% for ff = 1:nFOV
%     for dd = 1:nDays
%         for cc = 1:length(cellTypes)
%             % get current spatial decoding data
%             curDecode = squeeze(decodeUse.(cellTypes{cc})(ff,dd,:));
% 
%             % quantify distriubtion
%             dataQuant = distrQuant(curDecode',cueTemp,rewLoc,cueDil,maxBin);
% 
%             % calculate IO ratios
%             I = dataQuant.CueIn;
%             O = dataQuant.CueOut;
%             decodeIORaw = I/O;
%             decodeIONorm = (I-O)/(I+O);
% 
%             % store correlation values
%             decodingActivity.decodeIORawGlob.(cellTypes{cc}){ff,dd+1} = decodeIORaw;
%             decodingActivity.decodeIONormGlob.(cellTypes{cc}){ff,dd+1} = decodeIONorm;
%             decodingActivity.decodeIGlob.(cellTypes{cc}){ff,dd+1} = I;
%             decodingActivity.decodeOGlob.(cellTypes{cc}){ff,dd+1} = O;
% 
%             % loop through sexes
%             sexIDs = fieldnames(cellSelect)';
%             for ss = sexIDs
% 
%                 % loop through morphology
%                 morphTypes = fieldnames(cellSelect.(ss{:}).(cellTypes{cc}))';
%                 for mt = morphTypes
% 
%                     % loop through genotypes
%                     genotypes = fieldnames(cellSelect.(ss{:}).(cellTypes{cc}).(mt{:}))';
%                     for gg = genotypes
%                         % identify current cells
%                         curCells = cellSelect.(ss{:}).(cellTypes{cc}).(mt{:}).(gg{:}){ff,dd+1};
% 
%                         if isempty(curCells) || ~strcmp(mt{:},'allMorph')
%                             decodingActivity.decodeIORaw.(cellTypes{cc}).(ss{:}).(mt{:}).(gg{:}){ff,dd+1} = NaN;
%                             decodingActivity.decodeIONorm.(cellTypes{cc}).(ss{:}).(mt{:}).(gg{:}){ff,dd+1} = NaN;
%                             decodingActivity.decodeI.(cellTypes{cc}).(ss{:}).(mt{:}).(gg{:}){ff,dd+1} = NaN;
%                             decodingActivity.decodeO.(cellTypes{cc}).(ss{:}).(mt{:}).(gg{:}){ff,dd+1} = NaN;
%                         else
%                             % store correlation values
%                             decodingActivity.decodeIORaw.(cellTypes{cc}).(ss{:}).(mt{:}).(gg{:}){ff,dd+1} = decodeIORaw;
%                             decodingActivity.decodeIONorm.(cellTypes{cc}).(ss{:}).(mt{:}).(gg{:}){ff,dd+1} = decodeIONorm;
%                             decodingActivity.decodeI.(cellTypes{cc}).(ss{:}).(mt{:}).(gg{:}){ff,dd+1} = I;
%                             decodingActivity.decodeO.(cellTypes{cc}).(ss{:}).(mt{:}).(gg{:}){ff,dd+1} = O;
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% cd(p1)
% save('data/decodeInfoGlobal.mat','decodingActivity')

