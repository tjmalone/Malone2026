function [licking] = predLickingTwo_allRuns_allLicking(y,lickIdx,rewardIdx,distBeforeReward,cmRewardRelated,numberOfRuns,plotFig)

%this code is different from "predLickingOneOfTwo_allRuns.m'. This code
% considerS BOTH rewards and use the mean number of predictive licking as
%the reward predicive licking. It also consider all other runs on the track
%as "other track". This code is more specific to the two reward 10m track.
%rewardLoc1=535cm, rewardLoc2=35cm.


%This code analyzes the predictive licking of reward on linear track.
%this is different from predLickingOneOfTwo.m because that code
%only focus on the runs that there are licking, but this one look at all runs (with our without lick before reward). If ony licked after reward,
%doesn't count.

%input
%(1) y: y position of the animal. The data should be sampled at equal time
%intervals
%(2) rewardIdx:indices of reward delivered. This is the indices in y data
%(3) lickingIdx:indices of reward delivered. This is the indices in y data
%(4) loc1: the location that is closer to the reward. The lick is
%calculated between loc1 and reward, called predictive licking. The
%remaining licks between 0 and loc1 are other licks.
%if one reward is too close to the beginning of the track (within "distBeforeReward"), all distances
%before reward will be considered as predictive licking range.

%(5) numberOfRuns: this is to determine how many runs (counted from the
%beginning) are used for the analysis. The reason for doing this is that
%sometimes the animal's behavior does not look good after a while so only
%use the first couple of runs to evaluate the behavior. If numberOfRuns IS
%0, all the runs will be used. 0.5 means use half. Otherwise, we will just use the indicated
%number runs.
%(6) rewardLoc1: the location of one of the reward. This is the reward that
%will be analyzed
%(7) rewardLoc2: the location of the other reward. This one will not be
%analyzed

%Output:
% (1) Licking.y:y position
% (2) Licking.rewardIdx:reward indices amoung y
% (3) Licking.yReward:y position of rewards
% (4) Licking.lickIdx:lick indicis
% (5) Licking.loc1:location closer to reward
% (6) Licking.allNormPredLickIndividualRun:individual run normalized predicted lick
% (7) Licking.allNormOtherLickIndividualRun:individual run normalized non predictive lick
% (8) Licking.allNormPredLick:mean of allNormPredLickIndividualRun
% (9) Licking.allNormOtherLick:mean of allNormOtherLickIndividualRun
% (10) licking.lickIdxBeforeReward:lick idx in each run before reward
% (11) licking.lickYBeforeReward:lick y in each run before reward
% (12) Licking.isMorePredLicking:whether there are more licking
% (13) Licking.perPredLick:percentage of runs with licking before reward that showed predictive licking
% (14) licking.NPredLickAll: total number of pred lick
% (15) licking.NOtherLickAll: total number of other lick
% (16) licking.NormNPredLickAll: Norm total number of pred lick: normalized by
% distance
% (17) licking.NormNOtherLickAll: Norm total number of other lick: normalized by
% distance
% (18) licking.perPredInAll: predictive licking among all
% licking=licking.NPredLickAll/(licking.NPredLickAll+licking.NOtherLickAll)
% (19) licking.perPredInAllNorm: normalized predictive licking among all
% licking=licking.NormNPredLickAll/(licking.NormNPredLickAll+licking.NormNOtherLickAll)
%


%rewardLoc1 is more distal
%rewardLoc2 is more proximal
yReward=y(rewardIdx);
%
% %there must be two rewards, so figure out the location
% % [p,x]=ksdensity(yReward);
% % pp=findpeaks(p);
% % rewardLoc1=x(find(p==pp(2)));
% % rewardLoc2=x(find(p==pp(1)));
pp=round(yReward);
ppp=unique(pp);
NNumber=[];
for n=1:length(ppp);
    NNumber(n)=length(find(pp==ppp(n)));
end

[~,i]=sort(NNumber);

rewardLocs=ppp(i(end-1:end));


rewardLoc1=max(rewardLocs);
rewardLoc2=min(rewardLocs);


%first got the indices in each run
C=diff(y);
C=[0;C];

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
        % licking.loc1=loc1;%location closer to reward
        licking.totalNumberOfRuns=size(NN,1)-1;%real total number of runs
        licking.useNumberOfRuns=NaN;%used number of runs
        licking.allNormPredLickIndividualRun=NaN;%individual run normalized predicted lick
        licking.allNormOtherLickIndividualRun=NaN;%individual run normalized non predictive lick
        licking.allNormPredLick=NaN;%mean of allNormPredLickIndividualRun
        licking.allNormOtherLick=NaN;%mean of allNormOtherLickIndividualRun
        licking.lickIdxBeforeReward=NaN;%lick idx in each run before reward
        licking.lickYBeforeReward=NaN;%lick y in each run before reward
        licking.isMorePredLicking=NaN;%whether there are more licking
        licking.perPredLick=NaN;%percentage of runs with licking before reward that showed predictive licking
        
        licking.NPredLickAll=NaN;% total number of pred lick
        licking.NOtherLickAll=NaN;% total number of other lick
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

for m=1:size(useNN,1)-1;
    runIdx(m,1)=useNN(m,2)+1;%made change from 1 to 2
    runIdx(m,2)=useNN(m+1,1)-1;
end



% cmRewardRelated=50;%use 30cm after reward as the lick related to reward

%since it is too complicated to use real reward location as reward
%location at run by run basis, the following critiria all use the same
%number: rewardLoc1 and rewardLoc2
% if distBeforeReward>=rewardLoc2;%in our tast the reward 2 is very close to the biginning so it may have an issue. The reward 1 shouldn;'t matter
%these areas are related to predictive licking
locBeforeR1=max([rewardLoc1-distBeforeReward 0]);%this is the location before reward 1. Beyound this location is predictive licking
locBeforeR2=max([rewardLoc2-distBeforeReward 0]);

%do not include the previous run
%     locBeforeR2TrackEnd=max(y)-(disBeforeReward-rewardLoc2);%at the end of the track where the licking could occur early
%     locAfterR1=rewardLoc1+cmRewardRelated;
%     locAfterR1=rewardLoc2+cmRewardRelated;
% else %distBeforeReward<rewardLoc2;%in our tast the reward 2 is very close to the biginning so it may have an issue. The reward 1 shouldn;'t matter
%these areas are related to predictive licking
%      locBeforeR1=max([rewardLoc1-distBeforeReward 0]);%this is the location before reward 1. Beyound this location is predictive licking
%     locBeforeR2=rewardLoc2-distBeforeReward;
%do not include the previous run
%     locBeforeR2TrackEnd=nan;
%        locAfterR1=rewardLoc1+cmRewardRelated;
%     locAfterR1=rewardLoc2+cmRewardRelated;
% end

lickIdxPredR1={};
lickYPredR1={};
lickIdxPredR2={};
lickYPredR2={};
lickIdxOther={};
lickYOther={};
lickR1={};
lickR2={};
for m=1:size(runIdx,1);
    rIdx=rewardIdx(rewardIdx<=runIdx(m,2)&rewardIdx>=runIdx(m,1));
    lickIdxThisRun=lickIdx(lickIdx<=runIdx(m,2)&lickIdx>=runIdx(m,1));%lick indices this run
    lickYThisRun=y(lickIdxThisRun);%y position of these licks
    rIdxAll{m}=rIdx';
    
    if isempty(rIdx);%if there is no reward
        locBeforeR1Pred=rewardLoc1;
        locBeforeR2Pred=rewardLoc2;
        locAfterR1=rewardLoc1+cmRewardRelated;
        locAfterR2=rewardLoc2+cmRewardRelated;
        
        lickIdxPredR1{m}=lickIdxThisRun(lickYThisRun>=locBeforeR1 & lickYThisRun<locBeforeR1Pred)';
        lickYPredR1{m}=lickYThisRun(lickYThisRun>=locBeforeR1 & lickYThisRun<locBeforeR1Pred)';
        
        lickIdxPredR2{m}=lickIdxThisRun(lickYThisRun>=locBeforeR2 & lickYThisRun<locBeforeR2Pred)';
        lickYPredR2{m}=lickYThisRun(lickYThisRun>=locBeforeR2 & lickYThisRun<locBeforeR2Pred)';
        
        lickR1{m}=[];%reward 1 related lick
        lickR2{m}=[];%reward 2 related lick
        
        otherLickIdx=setdiff(lickIdxThisRun,[lickIdxPredR1{m} lickIdxPredR2{m} lickR1{m} lickR2{m}]);
        otherLickY=y(otherLickIdx);
        if size(otherLickIdx,1)~=1;
            lickIdxOther{m}=otherLickIdx';
        else
            lickIdxOther{m}=otherLickIdx;
        end
        
        if size(otherLickY,1)~=1;
            lickYOther{m}=otherLickY';
        else
            lickYOther{m}=otherLickY;
        end
        
    elseif length(rIdx)==2;%having both reward
        rY=y(rIdx);
        lickIdxPredR1{m}=lickIdxThisRun(lickYThisRun>=locBeforeR1 & lickIdxThisRun<max(rIdx))';
        lickYPredR1{m}=lickYThisRun(lickYThisRun>=locBeforeR1 & lickIdxThisRun<max(rIdx))';
        
        lickIdxPredR2{m}=lickIdxThisRun(lickYThisRun>=locBeforeR2 & lickIdxThisRun< min(rIdx))';
        lickYPredR2{m}=lickYThisRun(lickYThisRun>=locBeforeR2 & lickIdxThisRun< min(rIdx))';
        
        lickR1{m}=lickIdxThisRun(lickYThisRun<=(max(rY)+cmRewardRelated) & lickIdxThisRun> max(rIdx))';
        lickR2{m}=lickIdxThisRun(lickYThisRun<=(min(rY)+cmRewardRelated) & lickIdxThisRun> min(rIdx))';%reward 2 related lick
        
        otherLickIdx=setdiff(lickIdxThisRun,[lickIdxPredR1{m} lickIdxPredR2{m} lickR1{m} lickR2{m}]);
        otherLickY=y(otherLickIdx);
        if size(otherLickIdx,1)~=1;
            lickIdxOther{m}=otherLickIdx';
        else
            lickIdxOther{m}=otherLickIdx;
        end
        
        if size(otherLickY,1)~=1;
            lickYOther{m}=otherLickY';
        else
            lickYOther{m}=otherLickY;
        end
        
    else %now length(rIdx)==1; only one reward
        rY=y(rIdx);
        %now see which reward it is
        if abs(rY-rewardLoc1)<2 %this is reward 1
            locBeforeR2Pred=rewardLoc2;
            locAfterR2=rewardLoc2+cmRewardRelated;
            
            lickIdxPredR1{m}=lickIdxThisRun(lickYThisRun>=locBeforeR1 & lickIdxThisRun<rIdx)';
            lickYPredR1{m}=lickYThisRun(lickYThisRun>=locBeforeR1 & lickIdxThisRun<rIdx)';
            lickIdxPredR2{m}=lickIdxThisRun(lickYThisRun>=locBeforeR2 & lickYThisRun<locBeforeR2Pred)';
            lickYPredR2{m}=lickYThisRun(lickYThisRun>=locBeforeR2 & lickYThisRun<locBeforeR2Pred)';
            
            lickR1{m}=lickIdxThisRun(lickYThisRun<=(rY+cmRewardRelated) & lickIdxThisRun> rIdx)';
            lickR2{m}=[];%reward 2 related lick
            otherLickIdx=setdiff(lickIdxThisRun,[lickIdxPredR1{m} lickIdxPredR2{m} lickR1{m} lickR2{m}]);
            otherLickY=y(otherLickIdx);
            if size(otherLickIdx,1)~=1;
                lickIdxOther{m}=otherLickIdx';
            else
                lickIdxOther{m}=otherLickIdx;
            end
            
            if size(otherLickY,1)~=1;
                lickYOther{m}=otherLickY';
            else
                lickYOther{m}=otherLickY;
            end
            
        else %this is reward 2
            
            locBeforeR1Pred=rewardLoc1;
            locAfterR1=rewardLoc1+cmRewardRelated;
            
            lickIdxPredR1{m}=lickIdxThisRun(lickYThisRun>=locBeforeR1 & lickYThisRun<locBeforeR1Pred)';
            lickYPredR1{m}=lickYThisRun(lickYThisRun>=locBeforeR1 & lickYThisRun<locBeforeR1Pred)';
            
            lickIdxPredR2{m}=lickIdxThisRun(lickYThisRun>=locBeforeR2 & lickIdxThisRun< rIdx)';
            lickYPredR2{m}=lickYThisRun(lickYThisRun>=locBeforeR2 & lickIdxThisRun< rIdx)';
            
            lickR1{m}=[];%reward 1 related lick
            lickR2{m}=lickIdxThisRun(lickYThisRun<=(rY+cmRewardRelated) & lickIdxThisRun> rIdx)';%reward 2 related lick
            otherLickIdx=setdiff(lickIdxThisRun,[lickIdxPredR1{m} lickIdxPredR2{m} lickR1{m} lickR2{m}]);
            otherLickY=y(otherLickIdx);
            if size(otherLickIdx,1)~=1;
                lickIdxOther{m}=otherLickIdx';
            else
                lickIdxOther{m}=otherLickIdx;
            end
            
            if size(otherLickY,1)~=1;
                lickYOther{m}=otherLickY';
            else
                lickYOther{m}=otherLickY;
            end
        end
    end
    
    
end


%test whether the above numbers are correct:
% figure,
% plot(y,'k')
% hold on
% plot(cell2mat(lickIdxPredR1),cell2mat(lickYPredR1),'r.','MarkerSize',10);
% hold on
% plot(cell2mat(lickIdxPredR2),cell2mat(lickYPredR2),'m.','MarkerSize',10);
% hold on
% plot(cell2mat(lickIdxOther),cell2mat(lickYOther),'b.','MarkerSize',10);
% hold on
% plot(a.rewardIdx,y(a.rewardIdx),'g.','MarkerSize',10);
% hold on
% plot(cell2mat(lickR1),y(cell2mat(lickR1)),'y.','MarkerSize',10);
% hold on
% plot(cell2mat(lickR2),y(cell2mat(lickR2)),'y.','MarkerSize',10);


allNormPredLick=[];%one number, mean of all pred lick / length of pred lick track (between reward and loc1) for individual runs (mean of allNormPredLickIndividualRun);
allNormOtherLick=[];%one number, mean of none pred lick / length of pred lick track (between 0 and loc1)for individual runs (mean of allNormOtherLickIndividualRun);

allNormPredLickIndividualRun=[];%one number, all pred lick / length of pred lick track (between reward and loc1)
allNormOtherLickIndividualRun=[];%one number, all none pred lick / length of pred lick track (between 0 and loc1)

isMorePredLicking=[];

%go through data and get these numbers

for m=1:length(lickIdxPredR1);
    predLength=2*distBeforeReward;
    otherLength=max(y)-min(y)-predLength-2*cmRewardRelated;
    allNormPredLickIndividualRun(m)=(length(lickIdxPredR1{m})+length(lickIdxPredR2{m}))/predLength;
    allNormOtherLickIndividualRun(m)=length(lickIdxOther{m})/otherLength;
    
    if  allNormPredLickIndividualRun(m)>allNormOtherLickIndividualRun(m);
        isMorePredLicking(m,1)=1;
    else
        isMorePredLicking(m,1)=0;
    end
    
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

%%%%%%%%these are important%%%%%%%%%%%%
licking.lickIdxPredR1=lickIdxPredR1;%indices of predictive licking for reward 1
licking.lickYPredR1=lickYPredR1; %y posistions of predictive licking for reward 1
licking.lickIdxPredR2=lickIdxPredR2;%indices of predictive licking for reward 1
licking.lickYPredR2=lickYPredR2; %y posistions of predictive licking for reward 1
licking.lickR1=lickR1;%indices of post licking for reward 1
licking.lickR2=lickR2; %y posistions of post licking for reward 1
licking.lickIdxOther=lickIdxOther;%indices of non-predictive licking
licking.lickYOther=lickYOther;%y positions of non-redictive licking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

licking.isMorePredLicking=isMorePredLicking;%whether there are more licking
licking.perPredLick=perPredLick;%percentage of runs with licking before reward that showed predictive licking

licking.NPredLickAll=length(cell2mat(lickIdxPredR1))+length(cell2mat(lickIdxPredR2));% total number of pred lick: use mean of the two rewards
licking.NOtherLickAll=length(cell2mat(lickIdxOther));% total number of other lick


P=allNormPredLickIndividualRun(~isnan(allNormPredLickIndividualRun));%since L1 for each run is similar but slightly different, just add each run: will not be off too much
P=sum(P);

O=allNormOtherLickIndividualRun(~isnan(allNormOtherLickIndividualRun));%since L1 for each run is similar but slightly different, just add each run: will not be off too much
O=sum(O);

licking.NormNPredLickAll=P;% Norm total number of pred lick
licking.NormNOtherLickAll=O;% Norm total number of other lick
licking.perPredInAll= licking.NPredLickAll/(licking.NPredLickAll+licking.NOtherLickAll);%predictive licking among all licking
licking.perPredInAllNorm= P/(P+O);%predictive licking among all licking


if plotFig
    figure,
    plot(y,'k');
    hold on
    plot(rewardIdx,y(rewardIdx),'g.','MarkerSize',10)
    %plot loc1
    hold on
    line([1 length(y)],[rewardLoc1+cmRewardRelated rewardLoc1+cmRewardRelated],'Color','b');
    hold on
    line([1 length(y)],[rewardLoc2+cmRewardRelated  rewardLoc2+cmRewardRelated],'Color','b');
    hold on
    line([1 length(y)],[locBeforeR1 locBeforeR1],'Color','b');
    hold on
    line([1 length(y)],[locBeforeR1  locBeforeR1],'Color','b');
    %plot predictive licking and other licking
    hold on
    plot(cell2mat(lickIdxPredR1),cell2mat(lickYPredR1),'r.','MarkerSize',10);
    hold on
    plot(cell2mat(lickIdxPredR2),cell2mat(lickYPredR2),'m.','MarkerSize',10);
    hold on
    plot(cell2mat(lickIdxOther),cell2mat(lickYOther),'b.','MarkerSize',10);
    
    
    
    %plot the track with more predictive licking
    i=find(isMorePredLicking==1);
    rIdxAllPred=rIdxAll(i);
    predRewardIdx=cell2mat(rIdxAllPred);
    yPredReward=y(predRewardIdx);
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


