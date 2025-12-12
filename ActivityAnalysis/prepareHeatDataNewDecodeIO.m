%% Calculate in-cue versus out-cue decoding

clear; clc;
p1 = '/MATLAB Drive/FY2025/imagingData/';
cd(p1)

% load folder and index data
load('foldersLearning.mat')
load('data/cellSelect.mat','cellSelect')
cellSelect = cellSelect.learn;

% load data
useItr = 21;
useDay = 1:11;

load('data/decoding/decodingDataMaster_grid.mat','decodingDataItr')
decodeUse.grid = decodingDataItr(useItr).meanPerCorrectSpatial(:,useDay,:);

load('data/decoding/decodingDataMaster_nongrid.mat','decodingDataItr')
decodeUse.nongrid = decodingDataItr(useItr).meanPerCorrectSpatial(:,useDay,:);

load('data/decoding/decodingDataMaster_all.mat','decodingDataItr')
decodeUse.allType = decodingDataItr(useItr).meanPerCorrectSpatial(:,useDay,:);

% define fov and day number
nFOV = size(trueDays,1);
nDays = length(useDay);

% set parameters
binWidth = 5;
cueDil = 2;
maxBin = [];
cellTypes = {'allType','grid','nongrid'};

% load cue templates
cueFolder = '/MATLAB Drive/FY2025/AnalysisCode/PostAnalysis/Cues/6mEnv2/';
load([cueFolder 'tempL.mat'])
load([cueFolder 'tempR.mat'])
cueTemp = tempL | tempR;
rewLoc = [510 560]/5;

decodingActivity = struct();


%%

% loop through fov and day
for ff = 1:nFOV
    for dd = 1:nDays
        for cc = 1:length(cellTypes)
            % get current spatial decoding data
            curDecode = squeeze(decodeUse.(cellTypes{cc})(ff,dd,:));

            % quantify distriubtion
            dataQuant = distrQuant(curDecode',cueTemp,rewLoc,cueDil,maxBin);

            % calculate IO ratios
            I = dataQuant.CueIn;
            O = dataQuant.CueOut;
            decodeIORaw = I/O;
            decodeIONorm = (I-O)/(I+O);

            % calculate spatialized overall decoding
            All = mean(curDecode,'omitnan');

            % store correlation values
            decodingActivity.decodeIORawGlob.(cellTypes{cc}){ff,useDay(dd)} = decodeIORaw;
            decodingActivity.decodeIONormGlob.(cellTypes{cc}){ff,useDay(dd)} = decodeIONorm;
            decodingActivity.decodeIGlob.(cellTypes{cc}){ff,useDay(dd)} = I;
            decodingActivity.decodeOGlob.(cellTypes{cc}){ff,useDay(dd)} = O;
            decodingActivity.decodeAllGlob.(cellTypes{cc}){ff,useDay(dd)} = All;

            % loop through sexes
            sexIDs = fieldnames(cellSelect)';
            for ss = sexIDs

                % loop through morphology
                morphTypes = fieldnames(cellSelect.(ss{:}).(cellTypes{cc}))';
                for mt = morphTypes

                    % loop through genotypes
                    genotypes = fieldnames(cellSelect.(ss{:}).(cellTypes{cc}).(mt{:}))';
                    for gg = genotypes
                        % identify current cells
                        curCells = cellSelect.(ss{:}).(cellTypes{cc}).(mt{:}).(gg{:}){ff,useDay(dd)};

                        if isempty(curCells) || ~strcmp(mt{:},'allMorph')
                            decodingActivity.decodeIORaw.(cellTypes{cc}).(ss{:}).(mt{:}).(gg{:}){ff,useDay(dd)} = NaN;
                            decodingActivity.decodeIONorm.(cellTypes{cc}).(ss{:}).(mt{:}).(gg{:}){ff,useDay(dd)} = NaN;
                            decodingActivity.decodeI.(cellTypes{cc}).(ss{:}).(mt{:}).(gg{:}){ff,useDay(dd)} = NaN;
                            decodingActivity.decodeO.(cellTypes{cc}).(ss{:}).(mt{:}).(gg{:}){ff,useDay(dd)} = NaN;
                            decodingActivity.decodeAll.(cellTypes{cc}).(ss{:}).(mt{:}).(gg{:}){ff,useDay(dd)} = NaN;
                        else
                            % store correlation values
                            decodingActivity.decodeIORaw.(cellTypes{cc}).(ss{:}).(mt{:}).(gg{:}){ff,useDay(dd)} = decodeIORaw;
                            decodingActivity.decodeIONorm.(cellTypes{cc}).(ss{:}).(mt{:}).(gg{:}){ff,useDay(dd)} = decodeIONorm;
                            decodingActivity.decodeI.(cellTypes{cc}).(ss{:}).(mt{:}).(gg{:}){ff,useDay(dd)} = I;
                            decodingActivity.decodeO.(cellTypes{cc}).(ss{:}).(mt{:}).(gg{:}){ff,useDay(dd)} = O;
                            decodingActivity.decodeAll.(cellTypes{cc}).(ss{:}).(mt{:}).(gg{:}){ff,useDay(dd)} = All;
                        end
                    end
                end
            end
        end
    end
end

cd(p1)
save('data/decodeInfoGlobal_newDecodeIO.mat','decodingActivity')

