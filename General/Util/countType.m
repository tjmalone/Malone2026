function types = countType(cells,folder)
%% countType
%
% Determines the cells types of a given set of cells. Possible cell types
% include grid cells (excludes cue cells), cue cells, speed cells (excludes
% grid and cue cells), and other cells (all cells not of a defined type).
% Output is struct containing indices and percentages of each cell type.
%
% Inputs:
%   (1) cells: cell indices to analyze (default is all cells)
%
%   (2) folder: folder containing RunByRun_sig folder (default is cuurent
%       directory)
%


%% Process inputs

% define base folder
if nargin<2
    folder = pwd;
end

% move to desired folder
p = pwd;
cd(folder)

% identify cells
if nargin==0
    load('allROIs.mat','roi')
    cells = 1:size(roi,3);
end

% number of cells to analyze
cellN = length(cells);


%% Identify cell types

types = struct();

% extract cue cell indices
load('cueAnalysisNew_sig\newScoreShuffleTemplate\Left\cueCells.mat');
cueL = cueCells.cueCellRealIdx;
load('cueAnalysisNew_sig\newScoreShuffleTemplate\Right\cueCells.mat');
cueR = cueCells.cueCellRealIdx;
cueRL = unique([cueL;cueR]);

% extract grid cell indices
load('PValueClassifier_KY2_6_sig\grid.mat');
gridI = grid.indices;
gridI = setdiff(gridI,cueRL);

% extract speed cell indices
load('speed_dfof_sig\speed.mat');
speedI = unique([speed.speedCellPost;speed.speedCellNegt]);
speedI = setdiff(speedI,cueRL);
speedI = setdiff(speedI,gridI);

% determine "other" cells
otherI = cells;
otherI = setdiff(otherI,cueRL);
otherI = setdiff(otherI,gridI);
otherI = setdiff(otherI,speedI);


%% Save cell type indices and percentages

% save all indices
types.cue = intersect(cells,cueRL);
types.grid = intersect(cells,gridI);
types.speed = intersect(cells,speedI);
types.other = otherI;

% calculate percentage
types.cuePer = length(types.cue)/cellN;
types.gridPer = length(types.grid)/cellN;
types.speedPer = length(types.speed)/cellN;
types.otherPer = length(types.other)/cellN;


%%

cd(p)


end












