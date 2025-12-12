function speedScoreCalculation_TM(~,~)
% modified from original speedScoreCalculation code to vectorize
% speedScoreCalculatin for shuffles
%%
load('abfFake.mat');
m='dfof_*';
d=dir(m);
for n=1:length(d)
    load(d(n).name);
end
nShuffle = 1000;

folderNames = {'speed_dfof','speed_dfof_sig'};
suffix = {'dfof','dfofSig'};

for ff = 1:2
    if ff==1
        dfofCur = dfof;
    elseif ff==2
        dfofCur = dfof_sig;
    end
    
    mkdir(folderNames{ff});
    cd(folderNames{ff});
    clear speed

    [speed.scores,speed.useDfofFilter,speed.useSpeedFilter,speed.lowLim,speed.upLim] = speedScoreVector_TM(dfofCur,abfFake);
    
    figure,plot(speed.scores,'r.')
    xlabel('Cell IDs');
    ylabel('Speed scores');
    title(['Speed scores ' suffix{ff}]);
    
    saveas(gcf,'speedScores.fig');
    close
    
    speed.speedCellPost=[];
    speed.speedCellNegt=[];
    speed.shuffleScores=[];
    for n=1:size(speed.useDfofFilter,2)
        disp(n)
        Fuse=speed.useDfofFilter(:,n);
        
        I=length(Fuse);
        r=randi([ceil(I*0.05),ceil(I*0.95)],nShuffle,1);
        
        FShuffle = zeros(length(Fuse),nShuffle);
        for s=1:nShuffle%because we need to make 1000 random shffules.
            FShuffle(:,s)=[(Fuse(r(s):I));(Fuse(1:(r(s)-1)))];% put s data to the s column of the matrix.
        end
        
        % changed to calculate shuffled scores using the dfof and speed
        % outputs of speedScoreVector without repeating speed calculations
        % or dfof smoothing
        shuffleScores = corr(FShuffle,speed.useSpeedFilter)';
%         [shuffleScores,~,~,~,~]=speedScoreVector(FShuffle,abfFake);

        if speed.scores(n)>=prctile(shuffleScores,99);
            speed.speedCellPost(end+1,1)=n;
        elseif speed.scores(n)<=prctile(shuffleScores,1);
            speed.speedCellNegt(end+1,1)=n;
        end
        speed.shuffleScores(:,n)=shuffleScores;%each column has shuffles scores of one cell
    end
    speed.percentile=99;
    save('speed.mat','speed');
    
    cd ..\
end


