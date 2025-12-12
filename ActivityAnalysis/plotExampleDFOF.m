%% plotExampleRBR

clear; close all; clc

p1 = 'D:\AD_Project\imagingData';
cd(p1)

binWidth=5;
trackStart = 0;
trackEnd = 600;

% set save folder name
svFile = [p1 '\Figures\Figures' char(datetime('today','format','yyyyMMdd')) '\activityExamples'];
if ~isfolder(svFile)
    mkdir(svFile);
    mkdir([svFile '\data']);
end


%% plot high RBR consistency

cd('D:\AD_Project\imagingData\data\TM230303-2\230420\loc1_new\TSeries-886_1\suite2p')

useIdx = 1;
load('dfof_sig_391_cells.mat')
load('abfFake.mat','abfFake')
load('speedThreshold.mat','speedThreshold')

% calculate moving indices
t = abfFake.t;
speed = diff(abfFake.y);
fastEnough = [false ;speed>speedThreshold];

% get dfof
dfof_sig_mean = mean(dfof_sig,2,'omitnan')*100;
dfof_sig_use = dfof_sig(:,useIdx)*100;

figure
tiledlayout(2,1)
useData = {dfof_sig_use,dfof_sig_mean};
yMax = [400 10];
xRange = [120 220];
yLabs = {'\DeltaF/F (%)','Mean \DeltaF/F (%)'};

for ii= 1:2
    % define current data
    curData = useData{ii};

    % plot data
    nexttile(); hold on
    plot(t,curData)

    % plot stop regions
    inRegion = false;
    for jj = 1:length(fastEnough)
        if ~fastEnough(jj) && ~inRegion
            startIdx = jj;
            inRegion = true;
        elseif fastEnough(jj) && inRegion
            endIdx = jj - 1;
            patch([t(startIdx) t(jj-1) t(jj-1) t(startIdx)],[0 0 yMax(ii) yMax(ii)], ...
                [0.5 0.5 0.5], 'EdgeColor', 'none');
            inRegion = false;
        end
    end

    % set labels
    xlim(xRange)
    ylim([0 yMax(ii)])
    xlabel('Time (s)')
    ylabel(yLabs{ii})
end

% save figure
savefig([svFile '\dfof_example.fig'])

