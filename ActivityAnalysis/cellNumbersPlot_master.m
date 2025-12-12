%% cellNumbersGet
% Get the number of cells that meet various criteria and store data for
% generation of sex/cell type heatmap.
%

%% Normalize cell numbers

clear; close all; clc

p1 = 'D:\AD_Project\imagingData';
cd(p1)

% load active cell number data
load('D:\AD_Project\imagingData\data\cellNums.mat')
cellNumsAll = cellNums.All(:,1);

% load morphology information
load('analysis_ManualSelection\data_PyrSte.mat','dataPS');

% load mutaul index information
load('analysis_ManualSelection\infoMatch.mat','nCells','infoMatch')
useDay = 1;
nFOV = length(infoMatch);

% set data save directory
dataPath = 'data\';

% set save folder name
svFile = [p1 '\Figures\Figures' char(datetime('today','format','yyyyMMdd')) '\cellNumbers'];
if ~isfolder(svFile)
    mkdir(svFile);
    mkdir([svFile '\data']);
end

% define threshold parameters
threshCenter = 158;        % threshold for morphology center
threshRange = 20;            % threshold for morphology range
threshDist = 4;             % threshold for cell matching. Selected based on visualization
threshOvr = 0.4;            % threshold for matching area overlap

% current morphology edges
thLow = threshCenter-threshRange;
thHigh = threshCenter+threshRange;

% calculate morphology classes
morph = zeros(size(dataPS.allArea));
morph(dataPS.allArea<=thLow) = 2;
morph(dataPS.allArea>=thHigh) = 1;

% FOV scale factor (converts to um^2)
sizeScale = 1.47^2/1000^2; % m=m^2/pixel^2

% Load and scale FOV areas
% load FOV areas
load('D:\AD_Project\imagingData\data_FOV\FOVarea.mat')
areasNormActive = FOVareas.Active*sizeScale;
areasNormManual = FOVareas.Manual*sizeScale;

% base code
codeBase = ['C=' num2str(threshCenter) '_R='...
    num2str(threshRange) '_d=' num2str(threshDist) '_o=' num2str(threshOvr)];

% code for plot titles
codePlot = replace(codeBase,'_',', ');

% code for saving
codeSv = ['_' replace(erase(codeBase,'='),'.','-')];


%% Calculate overlap percentages

% initialize data structures for current parameter set
dataIdxs = struct();
dataN = struct();

% calculate man, act, com, ste, pyr, man/act, ste/act, pyr/act,
% ste/com, pyr/com

% loop through FOV
for ff = 1:nFOV
    % current FOV morphology classes
    diamCur = morph(dataPS.fovIdx==ff);

    % store cell indices for base categories
    dataIdxs.manual{ff,1} = 1:nCells.manual(ff);

    dataIdxs.stellate{ff,1} = find(diamCur==1);
    dataIdxs.pyramidal{ff,1} = find(diamCur==2);

    % current match information
    curMatch = infoMatch{ff}{useDay};

    if ~isempty(curMatch)
        % apply distance threshold
        useIdx = curMatch(:,4)<=threshDist & curMatch(:,5)>=threshOvr;

        % match stellate and pyramidal indices
        isStellate = ismember(curMatch(:,1),dataIdxs.stellate{ff});
        isPyramidal = ismember(curMatch(:,1),dataIdxs.pyramidal{ff});

        % store match indices
        dataIdxs.active{ff,1} = curMatch(useIdx & (isStellate | isPyramidal),2);
        dataIdxs.common{ff,1} = curMatch(useIdx & (isStellate | isPyramidal) & curMatch(:,3)==1,2);
        dataIdxs.Act_Man{ff,1} = curMatch(useIdx,2);
        dataIdxs.Act_Ste{ff,1} = curMatch(useIdx & isStellate,2);
        dataIdxs.Act_Pyr{ff,1} = curMatch(useIdx & isPyramidal,2);
        dataIdxs.Com_Man{ff,1} = curMatch(useIdx & curMatch(:,3)==1,2);
        dataIdxs.Com_Ste{ff,1} = curMatch(useIdx & isStellate & curMatch(:,3)==1,2);
        dataIdxs.Com_Pyr{ff,1} = curMatch(useIdx & isPyramidal & curMatch(:,3)==1,2);

    else
        dataIdxs.active{ff,1} = [];
        dataIdxs.common{ff,1} = [];
        dataIdxs.Act_Man{ff,1} = [];
        dataIdxs.Act_Ste{ff,1} = [];
        dataIdxs.Act_Pyr{ff,1} = [];
        dataIdxs.Com_Man{ff,1} = [];
        dataIdxs.Com_Ste{ff,1} = [];
        dataIdxs.Com_Pyr{ff,1} = [];
        dataField.Ste{ff,1} = [];
        dataField.Pyr{ff,1} = [];
    end

    fNames = fieldnames(dataIdxs);
    for mm = 1:length(fNames)
        dataN.(fNames{mm}) = cellfun(@length,dataIdxs.(fNames{mm}));
    end
end

% calculate percentages
perData = struct();
perData.ManAll = dataN.manual;
perData.ManSte = dataN.stellate;
perData.ManPyr = dataN.pyramidal;
perData.ManNormAll = dataN.manual./areasNormManual;
perData.ManNormSte = dataN.stellate./areasNormManual;
perData.ManNormPyr = dataN.pyramidal./areasNormManual;
perData.ActMan_ManAll = dataN.Act_Man./dataN.manual;
perData.ActMan_ManSte = dataN.Act_Ste./dataN.stellate;
perData.ActMan_ManPyr = dataN.Act_Pyr./dataN.pyramidal;
perData.ActManAll_Act = dataN.Act_Man./dataN.active;
perData.ActManSte_Act = dataN.Act_Ste./dataN.active;
perData.ActManPyr_Act = dataN.Act_Pyr./dataN.active;
perData.ActDif = perData.ActManPyr_Act./perData.ActManSte_Act;
perData.ComManAll_Com = dataN.Com_Man./dataN.common;
perData.ComManSte_Com = dataN.Com_Ste./dataN.common;
perData.ComManPyr_Com = dataN.Com_Pyr./dataN.common;
perData.ComDif = perData.ComManPyr_Com./perData.ComManSte_Com;
perData.ComDif(perData.ComDif==inf) = NaN;
perData.ActAll = cellNumsAll;
perData.ActSte = dataN.Act_Ste;
perData.ActPyr = dataN.Act_Pyr;
perData.ActNormAll = cellNumsAll./areasNormActive;
perData.ActNormSte = dataN.Act_Ste./areasNormActive;
perData.ActNormPyr = dataN.Act_Pyr./areasNormActive;

% save data
typeKey = {'# Manually Identified Cells','# Manually Identified Cells (norm. by FOV)',...
    'Percent of Manual Population that is Active','Active Cell Classification',...
    'Common Cell Classification','# Active Cells','# Active Cells (norm. by FOV)'};
save([dataPath 'ManActMatch' codeSv '.mat'],'dataIdxs','dataN','perData','typeKey')


%% Store data in heatmap struct

% load group information
load('groupIDs.mat','groups','groupIDs','sexes')
sexIDs = {'allSex','female','male'};
morphTypes = {'allMorph','ste','pyr'};
nGroups = length(groups);
nSexes = length(sexes);
nMorphs = length(morphTypes);

% initialize data struct
numData = struct();
fields = {'Man','ManNorm','ActMan_Man','ActType_Act','ComType_Com','Act','ActNorm'};
nFields = length(fields);
for ff = 1:nFields
    numData.(fields{ff}) = cell(nSexes,nMorphs,nGroups);
end

% set field call parts
fldStarts = {'Man','ManNorm','ActMan_Man','ActMan','ComMan','Act','ActNorm'};
fldEnds = {'','','','_Act','_Com','',''};
morphTags = {'All','Ste','Pyr'};

% loop through sexes
for ss = 1:nSexes

    % loop through morphologies
    for mm = 1:nMorphs

        % loop through genotypes
        for gg = 1:nGroups

            % current sexes
            curFOV = intersect(groups{gg},sexes{ss});

            for ff = 1:nFields
                curField = [fldStarts{ff} morphTags{mm} fldEnds{ff}];
                numData.(fields{ff}){ss,mm,gg} = perData.(curField)(curFOV);
            end
        end
    end
end


%% Plot heat maps

for ff = 1:nFields
    % define current data
    curMap = numData.(fields{ff});

    % make heat map
    h = heatmapGrid(curMap,['Sex',sexIDs],['Morphology',morphTypes]);
    sgtitle(typeKey{ff},'FontSize',18)

    % save figure
    savefig([svFile '\heatmap_num' fields{ff} '.fig'])
    save([svFile '\data\heatmap_num' fields{ff} '.mat'],'curMap')
end

