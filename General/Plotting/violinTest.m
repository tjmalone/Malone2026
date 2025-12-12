%% Define data

clear; close all; clc

groupNames = {'Weight','Goodness','Cuteness'};
nG = length(groupNames);

barNames = {'Cat','Dog','Rabbit','Badger'};
nB = length(barNames);

% define data and settings
means = [12, 30, 4, 19; 7, 10, 3, 1; 9, 2, 7, 4];
nRange = [20 100];
fctr = 1.5;

% define plot settings
pShow = generatePShow(3,[1 2;1 3;1 4;2 4]);
colors4 = {[1 0 0],[1 0 1],[0 1 1],[0 0 1]};
colors16 = repmat(colors4,3,1);
colors16(2,:) = cellfun(@(x) 0.7*x,colors16(2,:),'UniformOutput',false);
colors16(3,:) = cellfun(@(x) 0.4*x,colors16(3,:),'UniformOutput',false);


% initialize data
data = cell(nG,nB);

% define data
rng(42)
for gg = 1:nG
    for bb = 1:nB
        curN = randi(diff(nRange))+nRange(1);
        curMean = means(gg,bb);
        data{gg,bb} = normrnd(curMean,curMean/fctr,curN,1);
    end
end


%% Validate with bar graph

figure;
tiledlayout(1,3)
sgtitle('Bar Graph')

nexttile(1); hold on
barGroup(barNames,data(1,:),'bar',colors4)

nexttile(2); hold on
barGroup(groupNames,data,'bar',colors4,pShow)
legend(barNames)

nexttile(3); hold on
barGroup(groupNames,data,'bar',colors16,pShow)
legend(barNames)




%% Plot with violins

figure;
tiledlayout(1,3)
sgtitle('Violin Graph')

nexttile(1); hold on
barGroup(barNames,data(1,:),'violin',colors4)

nexttile(2); hold on
h = barGroup(groupNames,data,'violin',colors4,pShow);
legend(h,barNames)

nexttile(3); hold on
h = barGroup(groupNames,data,'violin',colors16,pShow);
legend(h,barNames)
