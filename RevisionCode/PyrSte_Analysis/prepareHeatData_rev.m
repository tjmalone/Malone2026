%%
% Prepares pre-calculated data for use in timecoursePlot_master.
% Specifically, process inter-day consistency, cue amplitude difference,
% and cue template correlations. Some calculations are based on previously
% run scripts, wihle some are directly calculated.
%


%% Generate cue template correlations

clear; close all; clc

p1 = 'D:\AD_Project\imagingData';
cd(p1)

% define output folder
folderOut = 'D:\AD_Project\Revision\PyrSteAnalysis\';

% load folder and index data
load('foldersLearning.mat')
load([folderOut 'cellSelect.mat'],'cellSelect')
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
        cd(curFolderName)

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

                            % process distribution
                            dataQuant = distrQuant(curDfof',tempRL,rewLocUse,cueDil,maxBin);

                            % calculate anchoring
                            I = dataQuant.CueIn;
                            O = dataQuant.CueOut;
                            dfofI = I;
                            dfofO = O;
                            dfofIO = I/O;
                        end
                        
                        % store correlation values
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
save([folderOut 'corrInfoGlobal.mat'],'dataActivity')

