%%

clear all; close all; clc

p = pwd;
mice = dir('ID_*');

nMice = length(mice);
pSlowing = cell(nMice,1);
pLicking = cell(nMice,1);

for jj = 1:length(mice)
    cd(mice(jj).name)
    
    sessions = dir('*.txt');
    sessions = sessions([1 3:end]);
    
    % 25 35 50 75, 30
    preRewLick=20;
    postRewLick = 30;
    
    % 50 100 150 200, 110
    binSlow = 100;
    postRewSlow = 30;
    trStart = 0;
    trEnd = 600;
    
    nRuns = 0;
    
    nSes = length(sessions);
    pSlowing{jj} = zeros(nSes,1);
    pLicking{jj} = zeros(nSes,1);
    
    for ii = 1:nSes
        curSes = sessions(ii).name;
        %     disp(curSes)
        
        logData = readLog(curSes, 9);
        Pr = logParams(logData);
        
        curSlow = speedChange_Percentile(curSes,binSlow,postRewSlow,...
            nRuns,trStart,trEnd);
        curLick = predLicking_allRuns_allLicking(Pr.y,Pr.lickIdx,...
            Pr.rewardIdx,preRewLick,postRewLick,nRuns,0);
        
        pSlowing{jj}(ii) = curSlow.meanCentile/100;
        pLicking{jj}(ii) = curLick.perPredInAll;
        
    end
    cd(p)
end


%%

days = 0:20;
useCols = [3 4;1 2];
catNames = {'WT','AD'};
colors = {'r','b'};

figure

for ii = 1:2
    if ii==1
        useData = pLicking;
        ttl = ['Lick thresh = ' num2str(preRewLick)];
    else
        useData = pSlowing;
        ttl = ['Slow thresh = ' num2str(binSlow)];
    end
    
    subplot(2,1,ii); hold on
    h = {};
    
    for jj = 1:size(useCols,1)
        for kk = 1:size(useCols,2)
            plot(days,useData{useCols(jj,kk)},[colors{jj} ':'],'LineWidth',1.3)
        end
        curMean = mean(cat(2,useData{useCols(jj,:)}),2);
        h{jj} = plot(days,curMean,colors{jj},'LineWidth',1);
    end
    
    title(ttl)
    legend([h{1}(1),h{2}(1)],catNames,'Location','eastoutside')
end

savefig(['behavior_' num2str(preRewLick) '_' num2str(binSlow)]);
% close


%%

% save('behavior.mat','pSlowing','pLicking')




