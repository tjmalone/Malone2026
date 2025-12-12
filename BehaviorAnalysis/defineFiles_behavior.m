%% defineFiles_behavior
% Allows user to define the various sub-categories of virmen sessions for
% behavior analysis. Defines files in 3 stages.
%   1: Direct file paths (files.mat)
%   2: Environment or behavior category for each file (fileTypes.mat)
%   3: More robust set of indices for various categories (fileType.mat)

clear; clc


%% Define file paths by beahvior analysis group

% behavior categories
fields = {'famLearn','novLearn','famRecall','novRecall'};

% manually define file selection
genotype = 'AD';
famLearn = 4;
novLearn = 5:14;
famRecall = 16;
novRecall = 17;

d = dir('*T*.txt');
dd = struct2cell(d);

% generate file struct
files = struct();

for ff = 1:length(fields)
    if ~isempty(eval(fields{ff}))
        files.(fields{ff}) = fullfile(pwd,dd(1,eval(fields{ff})))';
    else
        files.(fields{ff}) = {};
    end
end

files.genotype = genotype;

save('files.mat','files')


%% Define environment and behavior category for each file

clear; clc

% behavior categories
fields = {'famLearn','novLearn','famRecall','novRecall'};

load('files.mat')

d = dir('*T*.txt');
fNames = {d(:).name}';
nNames = length(fNames);

% set environment/recall type by behavior category
envType = zeros(nNames,1);
recallType = zeros(nNames,1);
envs = [1 2 1 2];
recalls = [0 0 1 1];

% generate type variables
for ff = 1:length(fields)
    useField = cellfun(@(x) erase(x,[pwd '\']),files.(fields{ff}),'UniformOutput',0);
    
    curGroup = ismember(fNames,useField);
    
    envType(curGroup) = envs(ff);
    recallType(curGroup) = recalls(ff);
end

save('fileTypes.mat','envType','recallType');


%% Generate robust index set

clear; clc

load('fileTypes.mat');

% behavior categories
typeFields = {'famLearn','novLearn','famRecall','novRecall','famPre'};

% set environment and recall type values
envIdx = [1 2 1 2 1];
recIdx = [0 0 1 1 -1];

% envType = [ones(3,1);envType];
envType(1:3) = 1;

% recallType = [-ones(3,1);recallType];
recallType(1:3) = -1;

% generate type indices
type = arrayfun(@(a,b) find(envType==a & recallType==b),...
    envIdx,recIdx,'UniformOutput',0);
type{5} = [type{1};type{5}];

% generate combined indices
typeLearn = type{2};
typeAll = unique(cat(1,type{:}));
 
% save variables
save('fileType.mat','envType','recallType','typeFields','type','typeLearn','typeAll');

