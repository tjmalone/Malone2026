%%

clear; close all; clc

files = dir('*.txt');
nFile = length(files);

pSlow = zeros(1,nFile);
pLick = zeros(1,nFile);

for ii = 1:nFile
    try
        [pSlow(ii),pLick(ii)] = plotLogData4m_optogenetics_training_final_TM(files(ii).name);
    catch
        pSlow(ii) = NaN;
        
        distBeforeReward=20;
        cmRewardRelated=30;
        numberOfRuns = 0;
        a=logParams_thre05(files(ii).name);
        d=predLicking_allRuns_allLicking_modified(a.y,a.lickIdx,a.rewardIdx,distBeforeReward,cmRewardRelated,numberOfRuns);
        pLick(ii) = d.perPredInAll;
    end
    close all
end

save('behavior_NT.mat','pSlow','pLick')

%%

figure; hold on

xidx = 1:nFile;

plot(xidx,pSlow,'LineWidth',2,'Color',[0,0.9,0]);
plot(xidx,pLick*100,'LineWidth',2,'Color',[0,0.4,0])

set(gca,'fontsize', 18);

hold off

ylim([0 100]);
xlim([0.9 xidx(end)+0.1]);

leg1 = legend('%slow','%lick','%correct',...
    'Orientation', 'horizontal', 'Location', 'northoutside');

xlabel('day');
ylabel('percent of runs')

savefig(figure(1), 'training_NT.fig');







