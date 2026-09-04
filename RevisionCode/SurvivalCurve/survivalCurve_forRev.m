%%

clear; close all; clc

p1 = 'Z:\labMembers\TM\AD_Project\survivalCurve';
cd(p1)


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

fileID = {'survivalLog_260625_forRev.xlsx','aliveMice\260625_mice_forRev.xlsx'};

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

figure

v1 = 2;
v2 = 3;

X1 = cell2mat(srvC(~cellfun('isempty',srvC(:,v1)),[1 v1]));
X1(:,2) = ~X1(:,2);

X2 = cell2mat(srvC(~cellfun('isempty',srvC(:,v2)),[1 v2]));
X2(:,2) = ~X2(:,2);

colors = {[0 0 1],[1 0 0]};
legends = {'WT','PS19'};
logrankTM(X1,X2,colors,legends,0.05,1)

xlim([0 17])


%% Generate suvival curves

nM = size(mice.age,1);

% useVar = 'genotype';
% useCol = {1:2,3:4,5,6};
% screenVar = 'gender';
% screenCol = 1:2;

useVar = 'gender';
useCol = {1,2};
screenVar = 'genotype';
screenCol = 3:4;

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

figure

v1 = 2;
v2 = 3;

X1 = cell2mat(srvC(~cellfun('isempty',srvC(:,v1)),[1 v1]));
X1(:,2) = ~X1(:,2);

X2 = cell2mat(srvC(~cellfun('isempty',srvC(:,v2)),[1 v2]));
X2(:,2) = ~X2(:,2);

colors = {[175 0 0]/255,[255 102 102]/255};
legends = {'Male','Female'};
logrankTM(X1,X2,colors,legends,0.05,1)

xlim([0 17])


%% Age restriction

figure

ageMin = 7.5;
ageMax = 9.5;

% filter early time points
Y1 = X1(X1(:,1)>=ageMin,:);
Y2 = X2(X2(:,1)>=ageMin,:);

% censor late time points
O1 = Y1(:,1)>=ageMax;
O2 = Y2(:,1)>=ageMax;
Y1(O1,1) = ageMax;
Y1(O1,2) = 1;
Y2(O2,1) = ageMax;
Y2(O2,2) = 1;

% set all events to early or late age
% Y1(~O1,1) = ageMin;
% Y2(~O2,1) = ageMax;

% remove late time point (only for testing)
% Y1(O1,:) = [];
% Y2(O2,:) = [];

% perform age-restircted log-range test
logrankTM(Y1,Y2,colors,legends,0.05,1)

xlim([ageMin ageMax])


% Calculate binned survival onset
edges = ageMin:0.5:ageMax;
edgeMean = mean([edges(1:end-1);edges(2:end)],1);
nEdge  = length(edges)-1;

dataX = {X1,X2};
dataOnset = nan(nEdge,2);

nAlive = [];
nPara = [];

for ee = 1:nEdge
    curBin = [edges(ee) edges(ee+1)];
    for xx = 1:2
        % get current data
        curX = dataX{xx};

        % find new paralysis
        curPara = inrange(curX(:,1),curBin) & curX(:,2)==0;

        % find all alive mice
        curAlive = curX(:,1)>curBin(1);

        % find paralysis onset percentage
        dataOnset(ee,xx) = sum(curPara)/sum(curAlive)*100;

        nAlive(ee,xx) = sum(curAlive);
        nPara(ee,xx) = sum(curPara);

    end
end


% plot paralysis onset rate
figure; hold on
for xx = 1:2
    plot(edgeMean,dataOnset(:,xx),'Color',colors{xx})
end

xlim([ageMin ageMax])

[~,p] = ttest(dataOnset(:,1),dataOnset(:,2));

title(['p = ' num2str(p,2)])

ylabel('Probability of developing paralysis (%)')


% add final statistics
nUnits = 'time bins';
testName = 'two-tailed paired Students t-test';
testPair = 1;
testMC = 0;
testLimitP = 0;
testCat = [];
outStats = [testCat ttestEffectSize(dataOnset(:,1),dataOnset(:,2),testName,nUnits,testPair,testMC,testLimitP)];

