%%

clear; clc


%% Initialize inputs

genoInDead = {'GCaMP <+/->; Tau P301S<WT/WT>',...
    'GCaMP <-/->; Tau P301S<WT/WT>',...
    'GCaMP <+/->; Tau P301S<*/WT>',...
    'GCaMP <-/->; Tau P301S<*/WT>',...
    'Tau P301S<WT/WT>',...
    'Tau P301S<*/WT>',...
    'GCaMP <+>; Tau P301S<WT/WT>'};

genoInAlive = {'GCaMP <+/->; Tau P301S<WT/WT>',...
    'GCaMP <->; Tau P301S<WT/WT>',...
    'GCaMP <+>; Tau P301S<*/WT>',...
    'GCaMP <->; Tau P301S<*/WT>',...
    'Tau P301S<WT/WT>',...
    'Tau P301S<*/WT>', ...
    'GCaMP <+>; Tau P301S<WT/WT>'};

genoOut = 1:7;

gendIn = {'Male','Female'};
gendOut = 1:2;

CoDpar = {'Hindlimb Paralysis','Slow Movement'};

fileID = {'survivalLog_250417_forFig.xlsx','aliveMice\250417_mice_forFig.xlsx'};

warning('off')

%% Load dead mice

fileCur = fileID{1};

keepCol = [2:3 5:7];
keepRow = [];

opts = detectImportOptions(fileCur);
data = readtable(fileCur,opts);

for ii = 1:size(data,1)
    if ~isempty(data.ID(4))
        keepRow(end+1) = ii;
    end
end

data = data(keepRow,keepCol);
nM = size(data,1);


%% Process dead mice info

mice = struct();
mice.genotype = zeros(nM,1);
mice.gender = zeros(nM,1);
mice.CoD = zeros(nM,1);
mice.age = zeros(nM,1);

for ii = 1:nM
    for jj = 1:length(genoOut)
        if strcmp(data.Genotype{ii},genoInDead{jj})
            mice.genotype(ii) = genoOut(jj);
        end
    end
    
    for jj = 1:length(gendOut)
        if strcmp(data.Gender{ii},gendIn{jj})
            mice.gender(ii) = gendOut(jj);
        end
    end
    
    for jj = 1:length(CoDpar)
        if strcmp(data.Condition{ii},CoDpar{jj})
            mice.CoD(ii) = 1;
            break
        end
    end
    
    mice.age(ii) = round(days(data.DOD(ii)-data.DOB(ii)));
end


%% Load living mice

fileCur = fileID{2};

keepCol = [4 6 8];
keepRow = [];

opts = detectImportOptions(fileCur);
data = readtable(fileCur,opts);

for ii = 1:size(data,1)
    if ~isempty(data.Use{ii}) && round(days(datetime("now")-data.DOB(ii)))>=26
        keepRow(end+1) = ii;
    end
end

data = data(keepRow,keepCol);
nMold = nM;
nM = size(data,1);


%% Process living mice info

for ii = 1:nM
    if strcmp(data.Genotype{ii},'')
        continue
    end
    
    nMold = nMold+1;
    useIdx = nMold;
    
    for jj = 1:length(genoOut)
        if strcmp(data.Genotype{ii},genoInAlive{jj})
            mice.genotype(useIdx) = genoOut(jj);
        end
    end
    
    for jj = 1:length(gendOut)
        if strcmp(data.Sex{ii},gendIn{jj})
            mice.gender(useIdx) = gendOut(jj);
        end
    end
    
    mice.CoD(useIdx) = 0;
    
    mice.age(useIdx) = round(days(datetime("now")-data.DOB(ii)));
end


%% Generate suvival curves

nM = size(mice.age,1);

useVar = 'genotype';
useCol = {1:2,3:4,5,6};
screenVar = 'gender';
screenCol = 1:2;

% useVar = 'gender';
% useCol = {1,2};
% screenVar = 'genotype';
% screenCol = 3:4;

nType = length(useCol);

srvC = cell(nM,nType+1);

row = 1;
for ii = 1:nM
    if ~ismember(mice.(screenVar)(ii),screenCol)
        continue
    end
    
    srvC{row,1} = mice.age(ii)/365*12;
    
    for jj = 1:nType
        if ismember(mice.(useVar)(ii),useCol{jj})
            srvC{row,jj+1} = mice.CoD(ii);
            break
        end
    end
    
    row = row+1;
end

srvC = srvC(1:row-1,:);



%% F1: WT vs. AD

v1 = 2;
v2 = 3;

X1 = cell2mat(srvC(~cellfun('isempty',srvC(:,v1)),[1 v1]));
X1(:,2) = ~X1(:,2);

X2 = cell2mat(srvC(~cellfun('isempty',srvC(:,v2)),[1 v2]));
X2(:,2) = ~X2(:,2);

logrank(X1,X2,0.05,1)


xlim([0 17])










