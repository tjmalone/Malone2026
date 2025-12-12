function [licking] = predLicking_allRuns_allLicking(y,lickIdx,rewardIdx,distBeforeReward,cmRewardRelated,numberOfRuns,figOn)

%this function is very similar to predLicking_allRuns.m but use all licks
%as "other licks". "predLicking_allRuns.m" only licks before reward.

%based on observation of many mice, we set that that lick after more than 25cm of the
%reward are considred as new licks (reward unrelated licks).

%This code analyzes the predictive licking of reward on linear track.
%this is different from predLicking.m because that code
%only focus on the runs that there are licking, but this one look at all runs (with our without lick before reward). If ony licked after reward,
%doesn't count.

%input
%(1) y: y position of the animal. The data should be sampled at equal time
%intervals
%(2) rewardIdx:indices of reward delivered. This is the indices in y data
%(3) lickingIdx:indices of reward delivered. This is the indices in y data
%(4) distBeforeReward: the distance before reward that is considered to be predictive licking. The lick is
%calculated within this diatance, called predictive licking. The
%remaining licks between 0 and loc1 (rewardLocation-distBeforeReward) and
%25 cm after reward are other licks.

%(5) numberOfRuns: this is to determine how many runs (counted from the
%beginning) are used for the analysis. The reason for doing this is that
%sometimes the animal's behavior does not look good after a while so only
%use the first couple of runs to evaluate the behavior. If numberOfRuns IS
%0, all the runs will be used. 0.5 means use half. Otherwise, we will just use the indicated
%number runs.

%Output:
% Licking.y:y position
% Licking.rewardIdx:reward indices amoung y
% Licking.yReward:y position of rewards
% Licking.lickIdx:lick indicis
% Licking.loc1:location closer to reward
% Licking.allNormPredLickIndividualRun:individual run normalized predicted lick
% Licking.allNormOtherLickIndividualRun:individual run normalized non predictive lick
% Licking.allNormPredLick:mean of allNormPredLickIndividualRun
% Licking.allNormOtherLick:mean of allNormOtherLickIndividualRun
% licking.lickIdxBeforeReward:lick idx in each run before reward
% licking.lickYBeforeReward:lick y in each run before reward
% Licking.isMorePredLicking:whether there are more licking
% Licking.perPredLick:percentage of runs with licking before reward that showed predictive licking
%
%
C=diff(y);
C=[C;0];

loc1=mean(y(rewardIdx))-distBeforeReward;
% cmRewardRelated=50;%use 30cm after reward as the lick related to reward

A=C<=-max(y)*0.2;%0 and 1. Indices of the teleportation point will be 1; in the analysis,
%in some previous resonant galvo data, there are multiple points during
%telepotation, for example, two points in teleportation and their speed
%difference was -300 and -700, so using this -max(y)*0.15 will mostly cover
%all these intermediate points. Therefore, contiguous points that are 1 need to be
%identified

N=contiguous(A,[1]);
NN=N{1,2};
NN=[[0 0];NN];
NN=[NN;[length(y)+1 length(y)+1]];

if numberOfRuns==0;%use all runs
    useNN=NN;
elseif numberOfRuns==0.5;%use half runs
    useNN=NN(1:ceil(size(NN,1)/2),:);
else
    if size(NN,1)<numberOfRuns+1;
        licking.y=y;%y position
        licking.rewardIdx=rewardIdx;%reward indices amoung y
        licking.yReward=y(rewardIdx);%y position of rewards
        licking.lickIdx=lickIdx;%lick indicis
        licking.distBeforeReward=distBeforeReward;%location closer to reward
        licking.cmRewardRelated=cmRewardRelated;
        licking.totalNumberOfRuns=size(NN,1)-1;%real total number of runs
        licking.useNumberOfRuns=numberOfRuns;%used number of runs
        licking.allNormPredLickIndividualRun=NaN;%individual run normalized predicted lick
        licking.allNormOtherLickIndividualRun=NaN;%individual run normalized non predictive lick
        licking.allNormPredLick=NaN;%mean of allNormPredLickIndividualRun
        licking.allNormOtherLick=NaN;%mean of allNormOtherLickIndividualRun
        licking.lickIdxBeforeReward=NaN;%lick idx in each run before reward
        licking.lickYBeforeReward=NaN;%lick y in each run before reward
        licking.lickIdxAfterReward=NaN;%lick idx in each run After reward
        licking.lickYAfterReward=NaN;%lick y in each run After reward
        
        %%%%%%%%these are important%%%%%%%%%%%%
        licking.lickIdxPred=NaN;%indices of predictive licking
        licking.lickYPred=NaN; %y posistions of predictive licking
        licking.lickIdxOther=NaN;%indices of non-predictive licking
        licking.lickYOther=NaN;%y positions of non-redictive licking
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        licking.isMorePredLicking=NaN;%whether there are more licking
        licking.perPredLick=NaN;%percentage of runs with licking before reward that showed predictive licking
        
        licking.NPredLickAll=NaN;% total number of pred lick
        licking.NOtherLickAll=NaN;% total number of other lick
        
        %
        % P=allNormPredLickIndividualRun(~isnan(allNormPredLickIndividualRun));%since L1 for each run is similar but slightly different, just add each run: will not be off too much
        % P=sum(P);
        %
        % O=allNormOtherLickIndividualRun(~isnan(allNormOtherLickIndividualRun));%since L1 for each run is similar but slightly different, just add each run: will not be off too much
        % O=sum(O);
        
        licking.NormNPredLickAll=NaN;% Norm total number of pred lick
        licking.NormNOtherLickAll=NaN;% Norm total number of other lick
        licking.perPredInAll= NaN;%predictive licking among all licking
        licking.perPredInAllNorm= NaN;%predictive licking among all licking
        
        return
    else
        useNN=NN(1:numberOfRuns+1,:);%use whatever number of runs determined by the user
    end
end


%find the indices within each run
runIdx=[];

NPredLickAll=0;%number of pred licks: all runs
NOtherLickAll=0;%number of other lick: all runs
for m=1:size(useNN,1)-1;
    runIdx(m,1)=useNN(m,2)+1;%made change from 1 to 2
    runIdx(m,2)=useNN(m+1,1)-1;
end

allNormPredLick=[];%one number, mean of all pred lick / length of pred lick track (between reward and loc1) for individual runs (mean of allNormPredLickIndividualRun);
allNormOtherLick=[];%one number, mean of none pred lick / length of pred lick track (between 0 and loc1)for individual runs (mean of allNormOtherLickIndividualRun);

allNormPredLickIndividualRun=[];%one number, all pred lick / length of pred lick track (between reward and loc1)
allNormOtherLickIndividualRun=[];%one number, all none pred lick / length of pred lick track (between 0 and loc1)

isMorePredLicking=[];
lickIdxBeforeReward={};
lickYBeforeReward={};
lickIdxAfterReward={};
lickYAfterReward={};

lickIdxPred={};
lickYPred={};
lickIdxOther={};
lickYOther={};


for m=1:size(runIdx,1);
    rIdx=rewardIdx(rewardIdx<=runIdx(m,2)&rewardIdx>=runIdx(m,1));
    lickIdxThisRun=lickIdx(lickIdx<=runIdx(m,2)&lickIdx>=runIdx(m,1));%lick indices this run
    lickYThisRun=y(lickIdxThisRun);%y position of these licks
    
    
    if ~isempty(rIdx);
        rY=y(rIdx);
    else
        rY=mean(y(rewardIdx)); %if there is no reward in the run, use the mean reward location
    end
    L1=rY-loc1;
    %      L2=loc1;
    
    endOfTrack=y(runIdx(m,2));
    if rY+cmRewardRelated>=endOfTrack; %use 25cm after reward as the lick related to reward
        L2=loc1; %there is no space to lick after reward, so just use before reward
    else
        L2=loc1+(endOfTrack-(rY+cmRewardRelated));
    end
    
    if ~isempty(rIdx);
        lickIdxBefore=lickIdx(lickIdx<=rIdx&lickIdx>=runIdx(m,1));%only focus on licks before reward
    else
        lickIdxBefore=lickIdxThisRun(lickYThisRun<=rY&lickIdxThisRun>=runIdx(m,1));%only focus on licks before reward
    end
    
    lickIdxAfter=lickIdxThisRun(lickYThisRun>(rY+cmRewardRelated)&lickIdxThisRun<=runIdx(m,2));
    allLicksThisRun=[lickIdxBefore;lickIdxAfter];
    
    
    if ~isempty(lickIdxBefore);
        lickYBefore=y(lickIdxBefore);
        
        lickIdxBeforeReward{m}=lickIdxBefore;
        lickYBeforeReward{m}=lickYBefore;
        lickIdxAfterReward{m}=lickIdxAfter;
        lickYAfterReward{m}=y(lickIdxAfter);
        
        NPredLick=length(find(lickYBefore>=loc1));%here it's ok to use location because we already know that these licks occured before reward delivery time
        %  NOtherLick=length(find(lickYBefore<loc1));
        NOtherLick=length(find(lickYBefore<loc1))+length(lickIdxAfter);
        
        NPredLickAll=NPredLickAll+ NPredLick;
        NOtherLickAll=NOtherLickAll+ NOtherLick;
        
        
        lickIdxPred{m}=lickIdxBefore(lickYBefore>=loc1)';
        lickYPred{m}=lickYBefore(lickYBefore>=loc1)';
        lickIdxOther{m}=[lickIdxBefore(lickYBefore<loc1);lickIdxAfter]';
        lickYOther{m}=[lickYBefore(lickYBefore<loc1);y(lickIdxAfter)]';
        
        NormPredLick=NPredLick/L1;%number of licks normalize by distance
        NormNOtherLick=NOtherLick/L2;%number of licks normalize by distance
        
        allNormPredLickIndividualRun(m)= NormPredLick;
        allNormOtherLickIndividualRun(m)= NormNOtherLick;
        
        if NormPredLick>NormNOtherLick;
            isMorePredLicking(m,1)=1;
        else
            isMorePredLicking(m,1)=0;
        end
        
    else %if lickIdxBefore=[];
        lickIdxBeforeReward{m}=[];
        lickYBeforeReward{m}=[];
        lickIdxAfterReward{m}=lickIdxAfter;
        lickYAfterReward{m}=y(lickIdxAfter);
        
        lickIdxPred{m}=[];
        lickYPred{m}=[];
        lickIdxOther{m}=[lickIdxAfter]';
        lickYOther{m}=[y(lickIdxAfter)]';
        
        
        NPredLick=0;
        NOtherLick=length(lickIdxAfter);
        NPredLickAll=NPredLickAll+ NPredLick;
        NOtherLickAll=NOtherLickAll+ NOtherLick;
        
        NormPredLick=NPredLick/L1;%number of licks normalize by distance
        NormNOtherLick=NOtherLick/L2;%number of licks normalize by distance
        
        
        allNormPredLickIndividualRun(m)= NormPredLick;
        allNormOtherLickIndividualRun(m)= NormNOtherLick;
        
        if NormPredLick>NormNOtherLick;
            isMorePredLicking(m,1)=1;
        else
            isMorePredLicking(m,1)=0;
        end
        
    end
    %  else
    %      lickIdxBeforeReward{m}=nan;
    %          lickYBeforeReward{m}=nan;
    %       allNormPredLickIndividualRun(m)= nan;
    %  allNormOtherLickIndividualRun(m)= nan;
    %      isMorePredLicking(m,1)=nan;
    %  end
end
a=allNormPredLickIndividualRun(~isnan(allNormPredLickIndividualRun));
b=allNormOtherLickIndividualRun(~isnan(allNormOtherLickIndividualRun));
allNormPredLick=mean(a);
allNormOtherLick=mean(b);
c=isMorePredLicking(~isnan(isMorePredLicking));
% perPredLick=length(find(c==1))/length(c);
perPredLick=length(find(c==1))/(size(useNN,1)-1);

licking.y=y;%y position
licking.rewardIdx=rewardIdx;%reward indices amoung y
licking.yReward=y(rewardIdx);%y position of rewards
licking.lickIdx=lickIdx;%lick indicis
licking.distBeforeReward=distBeforeReward;%location closer to reward
licking.cmRewardRelated=cmRewardRelated;
licking.totalNumberOfRuns=size(NN,1)-1;%real total number of runs
licking.useNumberOfRuns=size(useNN,1)-1;%used number of runs
licking.allNormPredLickIndividualRun=allNormPredLickIndividualRun;%individual run normalized predicted lick
licking.allNormOtherLickIndividualRun=allNormOtherLickIndividualRun;%individual run normalized non predictive lick
licking.allNormPredLick=allNormPredLick;%mean of allNormPredLickIndividualRun
licking.allNormOtherLick=allNormOtherLick;%mean of allNormOtherLickIndividualRun
licking.lickIdxBeforeReward=lickIdxBeforeReward;%lick idx in each run before reward
licking.lickYBeforeReward=lickYBeforeReward;%lick y in each run before reward
licking.lickIdxAfterReward=lickIdxAfterReward;%lick idx in each run After reward
licking.lickYAfterReward=lickYAfterReward;%lick y in each run After reward

%%%%%%%%these are important%%%%%%%%%%%%
licking.lickIdxPred=lickIdxPred;%indices of predictive licking
licking.lickYPred=lickYPred; %y posistions of predictive licking
licking.lickIdxOther=lickIdxOther;%indices of non-predictive licking
licking.lickYOther=lickYOther;%y positions of non-redictive licking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

licking.isMorePredLicking=isMorePredLicking;%whether there are more licking
licking.perPredLick=perPredLick;%percentage of runs with licking before reward that showed predictive licking

licking.NPredLickAll=NPredLickAll;% total number of pred lick
licking.NOtherLickAll=NOtherLickAll;% total number of other lick


P=allNormPredLickIndividualRun(~isnan(allNormPredLickIndividualRun));%since L1 for each run is similar but slightly different, just add each run: will not be off too much
P=sum(P);

O=allNormOtherLickIndividualRun(~isnan(allNormOtherLickIndividualRun));%since L1 for each run is similar but slightly different, just add each run: will not be off too much
O=sum(O);

licking.NormNPredLickAll=P;% Norm total number of pred lick
licking.NormNOtherLickAll=O;% Norm total number of other lick
licking.perPredInAll= NPredLickAll/(NPredLickAll+NOtherLickAll);%predictive licking among all licking
licking.perPredInAllNorm= P/(P+O);%predictive licking among all licking

if figOn
    figure,
    plot(y,'k');
    hold on
    plot(rewardIdx,y(rewardIdx),'g.','MarkerSize',10)
    %plot loc1
    hold on
    line([1 length(y)],[loc1 loc1],'Color','b');
    hold on
    line([1 length(y)],[mean(y(rewardIdx))+cmRewardRelated  mean(y(rewardIdx))+cmRewardRelated],'Color','b');
    %plot predictive licking
    % for m=1:size(runIdx,1);
    %     if ~isnan(lickYBeforeReward{m});
    %         a=find(lickYBeforeReward{m}>=loc1);
    %         b=find(lickYBeforeReward{m}<loc1);
    %         hold on
    %         plot(lickIdxBeforeReward{m}(a),lickYBeforeReward{m}(a),'r.','MarkerSize',10)
    %         hold on
    %         plot(lickIdxBeforeReward{m}(b),lickYBeforeReward{m}(b),'b.','MarkerSize',10)
    %         hold on
    %         plot(lickIdxAfterReward{m},lickYAfterReward{m},'b.','MarkerSize',10)
    %     end
    % end
    
    plot(cell2mat(lickIdxPred),cell2mat(lickYPred),'r.','MarkerSize',10);
    hold on
    plot(cell2mat(lickIdxOther),cell2mat(lickYOther),'b.','MarkerSize',10);
    %plot the track with more predictive licking
    i=find(isMorePredLicking==1);
    if max(i)>length(rewardIdx);
        i=i(1:end-1);
    end
    predRewardIdx=rewardIdx(i);
    yPredReward=licking.yReward(i);
    hold on
    plot(predRewardIdx,yPredReward,'mo','MarkerSize',10)
    
    
    a=allNormPredLick-allNormOtherLick;
    if a>0
        title(['all runs more pred licking',num2str(perPredLick*100),'% runs pred']);
    else
        title(['all runs no more pred licking',num2str(perPredLick*100),'% runs pred']);
    end
end

end
