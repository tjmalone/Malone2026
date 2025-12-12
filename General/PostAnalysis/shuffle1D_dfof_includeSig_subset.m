function[fields]=shuffle1D_dfof_includeSig_subset(cellNumber,imageStartNumber,imageEndNumber,abfStartNumber,abfEndNumber,trackStart,trackEnd,binWidth,speedThreshold,Pvalue1,Pvalue2,Pvalue3,dfof,dfof_sig,abf)
%this code included dfof_sig just to use the criteria in Dombeck paper that
%at least 20% runs in field have sig trans

%shuffled time for this one is 1000 times
%this code will eventrally give 1-P value
%this analysis is down based on dfof generated from the PCAICA code (manually selection of baseline) where the results of
%dfof fluorescence  are stored in the dfof
% keep only the images that were saved during clampex recording
% Pvalue1 is the upper threshold for the high number to be in field
% Pvalue2 is the upper threshold for the adjacent bin to be in field
% Pvalue3 is the lower threshold for the bin to be out of field
% in this code, Faverage is averaged activity based on segTime. not real
% fluorescence
% in and out fields were highlighted. In fields are red, out fields are
% gray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%general processing of the data: extract the data to be analyzed, and speed
%threshold the data



I=imageEndNumber-imageStartNumber+1;
Fuseallcell = dfof(imageStartNumber:imageEndNumber,:);
Fuse=Fuseallcell(:,cellNumber);

%some dfof are all NaNs

if isempty(find(~isnan(Fuse)))
    fields.dfofUse=Fuseallcell(:,cellNumber);
    dfofaverage=zeros((trackEnd-trackStart)/binWidth,1);
    fields.dfofaveragesmooth=dfofaverage;
    fields.dfofaveragesmoothFields=zeros(size(dfofaverage));
    fields.Pvalues=zeros((trackEnd-trackStart)/binWidth,1);
    fields.inFieldBins=NaN;
    fields.outFieldBins=NaN;
    fields.transitions=NaN;
    fields.percentageBins=NaN;
    fields.fieldStartEnds=NaN;
    fields.fieldCenters=NaN;
    fields.fieldSpacings=NaN;
    fields.fieldWidths=NaN;
    fields.meanInField=NaN;
    fields.meanOutField=NaN;
    fields.ratioInOut=NaN;
    return
end


if isempty(find(Fuse~=0))
    fields.dfofUse=Fuseallcell(:,cellNumber);
    dfofaverage=zeros((trackEnd-trackStart)/binWidth,1);
    fields.dfofaveragesmooth=dfofaverage;
    fields.dfofaveragesmoothFields=zeros(size(dfofaverage));
    fields.Pvalues=zeros((trackEnd-trackStart)/binWidth,1);
    fields.inFieldBins=NaN;
    fields.outFieldBins=NaN;
    fields.transitions=NaN;
    fields.percentageBins=NaN;
    fields.fieldStartEnds=NaN;
    fields.fieldCenters=NaN;
    fields.fieldSpacings=NaN;
    fields.fieldWidths=NaN;
    fields.meanInField=NaN;
    fields.meanOutField=NaN;
    fields.ratioInOut=NaN;
    return
end

% FuseSigallcell=dfof_sig(imageStartNumber:imageEndNumber,:);
% FuseSig=FuseSigallcell(:,cellNumber);
% %for the 18 track, add 900 to the track y position
%
if trackEnd==1800
    abf.y=abf.y+900;
end

% IDENTIFY WHEN RUNNING SPEED WAS FAST ENOUGH

% image indices in the abf file that correspond to relevant images
abfIms = abfStartNumber:abfEndNumber;
% speed at each index of abfIms
speed=diff(abf.y(abf.imageIndex(abfIms)));
% elements of abfIms during fast running
fastEnough=[false ;speed>speedThreshold];
% indices of the abf file during fast running
abfFastEnough = abf.imageIndex(abfIms(fastEnough));
trackLength=trackEnd-trackStart;
track=[trackStart:binWidth:trackEnd];
[y,nbin]=histc(abf.y(abf.imageIndex(abfStartNumber:abfEndNumber)),track);
%note: the nbin here could have 0 number if the mouse is in negative
%position

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%start making figures: two figures: (1) original averaged firing rate (2)
%in and out fields identified by P

%(1) original averaged firing rate

figure(1);clf;

% %heat map
% subplot(211)
% plotPointHeatMap(abf.y(abfFastEnough),...
%     1:length(abfFastEnough),...abf.t(abfFastEnough),...
%     Fuse(fastEnough,1)','plotSpec',gca);
% set(gca,'color',[0 0 0]);%set the background to be black
% %set(gca,'xtick',[],'ytick',[]);%remove labeling of all axis
% axis tight;
% %remove all lables on axis
% axis off
% %xlabel('track');
% %ylabel('time');
% title(['cell',num2str(cellNumber)]);%add cell number to the title
% %get rid of extra spaces
% %tightfig;
%
% %to maximize the figure
% %maxfig(f,1);

%plot run by run using the bined dfof
subplot(411)
[dfofM,dfofA,C]=plotdfofrunbyrun(cellNumber,imageStartNumber,imageEndNumber,abfStartNumber,abfEndNumber,binWidth,trackStart,trackEnd,speedThreshold,dfof,abf);
title(['cell',num2str(cellNumber),'   run by run']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%(2) in and out fields identified by P
%make first figure to be the figure for original firing rate

%get the bin data for original data set
for binNumber=1:trackLength/binWidth
    %make a parameter called Fcaverage
    dfofaverage(binNumber) = mean(Fuseallcell(nbin==binNumber&fastEnough,cellNumber));
end

binCenters=(binWidth/2:binWidth:trackLength);
gaussianWindow=gausswin(3*5/binWidth,1);
[dfofaverage]=gaussianSmoothWithNan(dfofaverage,gaussianWindow);
%plot Fcaverage and binCenters
% plot(binCenters,smoothFaverage,'black','LineWidth',2);
% Fmax=max(smoothFaverage);
% axis([0 trackLength 0 Fmax*1.1]);

%plot dfofaverage
subplot(412);
plot(binCenters,dfofaverage,'black','LineWidth',2);
Fmax=max(dfofaverage);
if dfofaverage==0
    axis([0 trackLength -1 1]);
elseif min(dfofaverage)==Fmax*1.1
    axis([0 trackLength -1 1]);
else
    axis([0 trackLength min(dfofaverage) Fmax*1.1]);
end
title('mean dfof');

if dfofaverage==0
    fields.dfofUse=Fuseallcell(:,cellNumber);
    fields.dfofaveragesmooth=dfofaverage;
    fields.dfofaveragesmoothFields=zeros(size(dfofaverage));
    fields.Pvalues=zeros((trackEnd-trackStart)/binWidth,1);
    fields.inFieldBins=NaN;
    fields.outFieldBins=NaN;
    fields.transitions=NaN;
    fields.percentageBins=NaN;
    fields.fieldStartEnds=NaN;
    fields.fieldCenters=NaN;
    fields.fieldSpacings=NaN;
    fields.fieldWidths=NaN;
    fields.meanInField=NaN;
    fields.meanOutField=NaN;
    fields.ratioInOut=NaN;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%then make the in and out field plot under the original firing rate on the
%same plot: based on comparison with shuffle data
%create shuffle data under certain bins
%create shuffle data from 1000 random shuffles
% generate a set of random numbers between 5%-95% of the data set
clear F
clear Fbin

r=randi([ceil(I*0.05),ceil(I*0.95)],1000,1);
%create a matrix that has 1000 columns and row number equals to the total
%number of frames
F=zeros(I,1000);
Fbin=zeros(trackLength/binWidth,1000);
for s=1:1000%because we need to make 1000 random shffules.
    %F(:,s)=[(Fuse(r(s):I));(Fuse(1:(r(s)-1)))];% put s data to the s column of the matrix.
    F(:,s)=circshift(Fuse,r(s));
end

for binNumber=1:trackLength/binWidth;
    %make a parameter called Fcaverage
    Fbin(binNumber,:) = mean(F(nbin==binNumber&fastEnough,:),1);
    %this Fbin should have all 100 shuffle data (in columns) and under this
    %bin conditions

    %smooth the shuffle

end

%     gaussianWindow=gausswin(3,1);
%make the gaussianWindow to be adaptable to the binWidth, for 5cm,
%gaussianwin(3,1), for 2.5cm binWidth, gaussianwin(6,1), for 1cm,
%gaussianwin(15,1);
gaussianWindow=gausswin(3*5/binWidth,1);
for s=1:1000
    [Fbin(:,s)]=gaussianSmoothWithNan(Fbin(:,s),gaussianWindow);
end
%make the original fluorescent bined data as one column data
Fori=dfofaverage';
%now is to compare each number in Fori with Fbin in the same row
%first create a matrix the row number is the bin numbers and column number
%is the shuffle time: 1000
E=zeros(trackLength/binWidth,1000);
for v=1:trackLength/binWidth;
    %compare and get a logic array that each shuffle data in the bin if it
    %is>the one in the original data, the number is 1. so the E is a
    %matrix containing 1 and 0.
    E(v,:)=Fbin(v,:)>=Fori(v,:);
    %then the P value is the percentage of bin data that is larger than the one
    %in the original data, means the percentage of 1 in each row. Then 1-p
    %value is 1-this percentage number
    P1=1-sum(E,2)/1000; %'2' means the second dimension in the table E.
    %so here P1 is the 1-p value
end

% %smooth P value: this is changed step for the P value in 20150327_2
% gaussianWindow=gausswin(3*5/binWidth,1);
% [P1]=gaussianSmoothWithNan(P1,gaussianWindow);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get in and out field bins based on P value
PInField=[P1>=Pvalue1];
POutField=[P1<=Pvalue3];

%%%%%%%%%%%%%%%%%%%%%%%
%test wether there is any PInField
if isempty(find(PInField,1));
    fields.dfofUse=Fuseallcell(:,cellNumber);
    fields.dfofaveragesmooth=dfofaverage;
    fields.dfofaveragesmoothFields=zeros(size(dfofaverage));
    fields.Pvalues=P1;
    fields.inFieldBins=NaN;
    fields.outFieldBins=NaN;
    fields.transitions=NaN;
    fields.percentageBins=NaN;
    fields.fieldStartEnds=NaN;
    fields.fieldCenters=NaN;
    fields.fieldSpacings=NaN;
    fields.fieldWidths=NaN;
    fields.meanInField=NaN;
    fields.meanOutField=NaN;
    fields.ratioInOut=NaN;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot in and out fields and binCenters
% figure,plot(binCenters,PInField,'r',binCenters,POutField,'b');
% ylim([0,1.2]);

%next also include the fields that next to these ones with higher P value
%and the P values are just bigger than 0.7 (or other number, defined by
%Pvalue2

%first find all the 1 in PInField
[rowIn,colIn]=find(PInField);
%make a new matrix for extended fields
PInFieldNew=PInField;
for p=1:length(rowIn);
    if rowIn(p)==1;
        PInFieldNew(rowIn(p))=1;
    elseif P1(rowIn(p)-1)>=Pvalue2;
        PInFieldNew(rowIn(p)-1)=1;
    else PInFieldNew(rowIn(p)-1)=0;
    end
    if rowIn(p)==trackLength/binWidth;
        PInFieldNew(rowIn(p))=1;
    elseif P1(rowIn(p)+1)>=Pvalue2;
        PInFieldNew(rowIn(p)+1)=1;
    else PInFieldNew(rowIn(p)+1)=0;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%next need to see wether there is really continuous 1s in PInFieldNew,
%because the numbers in the current PInFieldNew (even the previous
%PInField), could be individual 1s that are not continous (like 101010001,
%not 11100011). If do not do this first, the 'contiguous' function below (runsN=contiguous(PInFieldNew,[1]);)will show error.
%If do this test, there will be no erros from contigous but will return to
%empty output. This empty output will be used in the final classifier
%'isGridCell'
[rowNew,colNew]=find(PInFieldNew);
D=diff(rowNew);
%if there are continuous numbers, their must be 1 in D. the intersect below
%will tell whether there is any number in D that contains 1
if (isempty(intersect(1,D)));
    fields.dfofUse=Fuseallcell(:,cellNumber);
    fields.dfofaveragesmooth=dfofaverage;
    fields.dfofaveragesmoothFields=zeros(size(dfofaverage));
    fields.Pvalues=P1;
    fields.inFieldBins=NaN;
    fields.outFieldBins=NaN;
    fields.transitions=NaN;
    fields.percentageBins=NaN;
    fields.fieldStartEnds=NaN;
    fields.fieldCenters=NaN;
    fields.fieldSpacings=NaN;
    fields.fieldWidths=NaN;
    fields.meanInField=NaN;
    fields.meanOutField=NaN;
    fields.ratioInOut=NaN;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this codes below was from previous ones. this code is wrong because it started from finding continuous 1s, which ignores individual 1.This is not right. Defining whether individual field has enough bins should be done after expanding the field.

%%next also include the fields that next to these ones with higher P value
%and the P values are just bigger than 0.7 (or other number, defined by
%Pvalue2)
% runs=contiguous(PInField,[1]);
% %this contiguous function tells the boundary of raws that are continuous for one number;
% allInField=runs{1,2};%get the 2-column data in the runs. Each row is 1 block of continuous 1, the first column is the smaller boundary, the second column is the bigger boundary
% [m,n] = size(allInField);%this tells the number of rows of the allInField, basically m tells how many continous blocks that all contains 1.
% fieldStart=allInField(:,1);
% fieldEnd=allInField(:,2);
% PInFieldNew=PInField;
% for a=1:m;
%     if fieldStart(a)==1;
%         PInFieldNew(fieldStart(a))=1;
%     elseif P1(fieldStart(a)-1)>=Pvalue2;
%         PInFieldNew(fieldStart(a)-1)=1;
%     else PInFieldNew(fieldStart(a)-1)=0;
%     end
%     if fieldEnd(a)==trackLength/binWidth;
%         PInFieldNew(fieldEnd(a))=1;
%     elseif P1(fieldEnd(a)+1)>=Pvalue2;
%         PInFieldNew(fieldEnd(a)+1)=1;
%     else PInFieldNew(fieldEnd(a)+1)=0;
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%make an empty row vector to further identify in field features:(1) three
%or more continuous bins; one exception is: except at the beginning and end of the track
%where two adjacent bins were sufficient. This exception can be taken care
%by the function contiguous, which already finds continous 1s (more than
%1). So all the field here should be more than 1 bins
%first recheck the contiguous 1 numbers in the new PInField (PInFieldNew)
runsN=contiguous(PInFieldNew,[1]);

%%%%%%%%%%%%%%%%%%%
% %if there is no continuous 1 in PInFieldNew, this code will show error like
% %this:
% Attempted to access indexVect(0); index must be a positive
% integer or logical.
%
% Error in contiguous (line 51)
%     shiftVect = [indexVect(2:end);indexVect(end)];
%
% 51      shiftVect = [indexVect(2:end);indexVect(end)];

%but based on the current code, there shouldn't be errors, because whether
%there are continous 1s has been tested above already (see above)
%%%%%%%%%%%%%%%%%%%%

allInFieldN=runsN{1,2};
[mN,nN] = size(allInFieldN);
fieldStartN=allInFieldN(:,1);
fieldEndN=allInFieldN(:,2);
for aN=1:mN;
    if fieldStartN(aN)==1;%accept 2 bin field under 5cm and 4bin field under 2.5cm bins if it is at the start of the track
        if fieldEndN(aN)-fieldStartN(aN)>=2*5/binWidth-1;
            PInFieldNew(fieldStartN(aN):fieldEndN(aN))=1;
        else PInFieldNew(fieldStartN(aN):fieldEndN(aN))=0;
        end
    elseif fieldEndN(aN)==length(PInFieldNew);%accept 2 bin field if it is at the end of the track
        if fieldEndN(aN)-fieldStartN(aN)>=2*5/binWidth-1;
            PInFieldNew(fieldStartN(aN):fieldEndN(aN))=1;
        else PInFieldNew(fieldStartN(aN):fieldEndN(aN))=0;
        end

        %below makes sure that for others, there will be at least 3 bins under 5cm bin, at least 6 bins under 2.5cm bin (maks this number to be adaptable to binWidth by equation).
        %     elseif fieldEndN(aN)-fieldStartN(aN)>2;
    elseif fieldEndN(aN)-fieldStartN(aN)>=3*5/binWidth-1;
        PInFieldNew(fieldStartN(aN):fieldEndN(aN))=1;
    else
        PInFieldNew(fieldStartN(aN):fieldEndN(aN))=0;
    end
end

%check whehter PInFieldNew is empty in this step
if isempty(find(PInFieldNew,1));
    fields.dfofUse=Fuseallcell(:,cellNumber);
    fields.dfofaveragesmooth=dfofaverage;
    fields.dfofaveragesmoothFields=zeros(size(dfofaverage));
    fields.Pvalues=P1;
    fields.inFieldBins=NaN;
    fields.outFieldBins=NaN;
    fields.transitions=NaN;
    fields.percentageBins=NaN;
    fields.fieldStartEnds=NaN;
    fields.fieldCenters=NaN;
    fields.fieldSpacings=NaN;
    fields.fieldWidths=NaN;
    fields.meanInField=NaN;
    fields.meanOutField=NaN;
    fields.ratioInOut=NaN;
    return
end

%(1) if want to plot P values and Faverage in the same plot, uncomment below:
% %scale PInFieldNew and POutField to fit the same plot of Faverage
% PInFieldNewS=zeros(size(PInFieldNew,1),1);
% POutFieldS=zeros(size(POutField,1),1);
% for n=1:size(PInFieldNew,1);
%     if PInFieldNew(n)==0;
%         PInFieldNewS(n)=min(Faverage);
%     else
%         PInFieldNewS(n)=max(Faverage);
%     end
% end
%
% for n=1:size(POutField,1);
%     if POutField(n)==0;
%         POutFieldS(n)=min(Faverage);
%     else
%         POutFieldS(n)=max(Faverage);
%     end
% end
%
% plot(binCenters,Faverage,'black','LineWidth',2);
% Fmax=max(Faverage);
% axis([0 trackLength min(Faverage) Fmax*1.1]);
% hold on
% plot(binCenters,PInFieldNewS,'r',binCenters,POutFieldS,'b');

% (2)If want to plot patch colors under in and out field places
%find the continuous 1 in new PInFieldNew again (this PInFieldNew has been
%updated in for loop from line 150.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%next need to test whether in the 10% if runs there is at leat one binned dfof in in-field bins are above
%average. to do this, use the dfofM from the function: plotdfofrunbyfun

%first renew the field information again
clear runsN
clear allInFieldN
clear fieldStartN
clear fieldEndN
runsN=contiguous(PInFieldNew,[1]);
allInFieldN=runsN{1,2};
[mN,nN] = size(allInFieldN);
fieldStartN=allInFieldN(:,1);
fieldEndN=allInFieldN(:,2);

%%%% this is the old version for dfof
% %mean of dfof in all runs: should calculate just the numbers in dfofaverage
% %that are not NaN
% %this is old threshold.
% % threshold=mean(mean(dfofaverage(~isnan(dfofaverage))))+std(dfofaverage(~isnan(dfofaverage)));
%
% %this is new threshold: mean dfof of all out field bins
% threshold=mean(dfofaverage(find(POutField,1)));
%
% for n=1:mN,%for each field
%     for m=1:size(dfofM,1);%for each run (it may include the runs that unfinished, so for each field, the number of runs is different), find whether there is at least one bin (one field includes many bins) that its mean dfof higher than meandfofaverage of out field bins
%         dfofRunByRun=find(dfofM(m,allInFieldN(n,1):allInFieldN(n,2))>=threshold,1);
%         dfofRunByRunBig(m,n)=(~isempty(dfofRunByRun));%so this is 0 or 1 numbers. Next step is to see whether there are move than 20% runs like that.
% %this line below is to check whether there is >10% of runs that averaged dfof in
% %the field (mean in all the bins within that field) is above mean dfofaverage. if we do this, a lot of cell may not
% %pass because they do not have a lot of this kind of runs
% % dfofRunByRunBig(m,n)=mean(dfofM(m,allInFieldN(n,1):allInFieldN(n,2)))>=meandfofaverage;
%     end
% end


%new version using dfof_sig

%so this dfofM has nan when speed doesn't reach the threshold
[dfofM_sig]=getRunByRunActivityNoSmooth(cellNumber,cellNumber,imageStartNumber,imageEndNumber,abfStartNumber,abfEndNumber,binWidth,trackStart,trackEnd,speedThreshold,dfof_sig,abf);

dfofM_sig=dfofM_sig{1};
dfofRunByRunBig=[];
for n=1:mN,%for each field
    for m=1:size(dfofM_sig,1);%for each run (it may include the runs that unfinished, so for each field, the number of runs is different), find whether there is at least one bin (one field includes many bins) having at least one siginificant transient
        dfofRunByRun=find(dfofM_sig(m,allInFieldN(n,1):allInFieldN(n,2))>0,1);%when there is a significant transient, it will be bigger than 0
        dfofRunByRunBig(m,n)=(~isempty(dfofRunByRun));%so this is 0 or 1 numbers. Next step is to see whether there are move than 20% runs like that.
        %this line below is to check whether there is >20% of runs that averaged dfof in
        %the field (mean in all the bins within that field) is above mean dfofaverage. if we do this, a lot of cell may not
        %pass because they do not have a lot of this kind of runs
        % dfofRunByRunBig(m,n)=mean(dfofM(m,allInFieldN(n,1):allInFieldN(n,2)))>=meandfofaverage;
    end
end

%Now, the next step is to calcuate the numner of runs for each bin,
%this is to find the telepore point where the speed is a big negative
[row,col]=find(speed<-100);
startIndicesP=[1 (row+1)']';
endIndicesP=[row' abfEndNumber]';

%because sometimes the fast aquisation then average of every three frames to get the y will cause that there are two points in the middle during teleportation, so these numbers are larger than real run number, solve it like this
startIndices=startIndicesP(diff(startIndicesP)>1);
startIndices=[startIndices' startIndicesP(end)]';
%do the problem happens to the endIndices, solve it like this
endIndices=endIndicesP(diff([1 endIndicesP'])>1);


binRuns=zeros(size(endIndices,1),trackLength/binWidth);
for r=1:size(endIndices,1);
    %get all the binNumber for this run
    nThisRun=nbin(startIndices(r,1):endIndices(r,1));%find all the n within this run
    nThisRun=nThisRun(find(nThisRun));%here remove all the 0
    %get the smallerst bin number
    minnThisRun=min(nThisRun);
    %get the largest bin number
    maxnThisRun=max(nThisRun);
    %all 1s in the bins that occupied
    binRuns(r,nThisRun)=1;
    %this is correct. this one has many 0, when there is nan (below
    %threshold) or beginning or end of track, it will show 0, otherwise
    %this is 1

    %below is not correct
    %     binRuns(r,minnThisRun:maxnThisRun)=1;
    %this binRuns will only reveal that in one run, what range of bins were
    %covered (basically to find at the beginning and end of each run  what bins were missed due to the incompletion of the field), if there are nan's in the middle, it will not see. but this
    %can be revealed in the next step: see right below
end

%now check whether for each field, whether there is at least 10% of runs where there is at least one bin the dfof is higher than meandfofaverage
%in the original cristina paper, this percentile is 20% but those are
%spikes. GCaMP6f willd definitly miss traisnents


for n=1:mN;
    if isempty(dfofRunByRunBig);
        NRunspass=0;
    else
        NRunspass=length(find(dfofRunByRunBig(:,n)));%get number of runs that having at least one bin bigger than threshold
    end
    NTotalRuns=0;
    binRunsThisField=binRuns(:,allInFieldN(n,1):allInFieldN(n,2));
    for m=1:size(binRunsThisField,1);
        %if all numbers are 0, discard this run
        if isempty(find(binRunsThisField(m,:)));
            NTotalRuns=NTotalRuns;
        else
            NTotalRuns=NTotalRuns+1;
        end
    end
    %     NTotalRuns=min(length(find(binRuns(:,allInFieldN(n,1)))),length(find(binRuns(:,allInFieldN(n,2)))));
    %

    perRuns=NRunspass/NTotalRuns;
    % changed from 0.1 to 0.2 on 3/9/23
    if perRuns>=0.2;
        PInFieldNew(allInFieldN(n,1):allInFieldN(n,2),1)=1;
    else
        PInFieldNew(allInFieldN(n,1):allInFieldN(n,2),1)=0;
    end
end

%check whehter PInFieldNew is empty in this step
if isempty(find(PInFieldNew,1));
    fields.dfofUse=Fuseallcell(:,cellNumber);
    fields.dfofaveragesmooth=dfofaverage;
    fields.dfofaveragesmoothFields=zeros(size(dfofaverage));
    fields.Pvalues=P1;
    fields.inFieldBins=NaN;
    fields.outFieldBins=NaN;
    fields.transitions=NaN;
    fields.percentageBins=NaN;
    fields.fieldStartEnds=NaN;
    fields.fieldCenters=NaN;
    fields.fieldSpacings=NaN;
    fields.fieldWidths=NaN;
    fields.meanInField=NaN;
    fields.meanOutField=NaN;
    fields.ratioInOut=NaN;
    return
end


%now end of testing in fields

%%%%%%%%%%%%%%%%%%%%%%

%indeces of in fields
inFieldBins=find(PInFieldNew);

%convert dfofaverage (smoothed by gaussian window as did above) to
%the dfofaveragesmoothFields, which only in fields keep its own numbers,
%out fields and unidentified fields are all 0.
dfofaveragesmoothFields=zeros(1,length(dfofaverage));
for n=1:length(inFieldBins);
    dfofaveragesmoothFields(1,inFieldBins(n,1))=dfofaverage(1,inFieldBins(n,1));
end

%if in the previous step, the PInFieldNew is not empty, now there must be
%contiguous numbers in PInFieldNew:
clear runsN
clear allInFieldN
clear mN
clear nN

runsN=contiguous(PInFieldNew,[1]);
allInFieldN=runsN{1,2};
[mN,nN] = size(allInFieldN);
%first column in allInField is field start, second column in allInField is
%field end
% fieldStartN=allInFieldN(:,1);
% fieldEndN=allInFieldN(:,2);

%calculate field Centers
fieldCenters=zeros(mN,1);
%calculate field width
fieldWidths=zeros(mN,1);
for n=1:mN;
    %this is field start,this is the start point of the bins
    fieldStartEnd{n}(1,1)=(allInFieldN(n,1)-1)*binWidth;
    %this is field end, this is the end point of the bins
    fieldStartEnd{n}(2,1)=allInFieldN(n,2)*binWidth;
    fieldCenters(n,1)=mean(fieldStartEnd{n});
    fieldWidths(n,1)=(allInFieldN(n,2)-allInFieldN(n,1))*binWidth;
end

%calculate the ratio of mean in field firing ratio
%firing rate
%but before doing this, should normalize dfofaverage so that it doesn't
%contain negative numbers. otherwise, mOutField can be negative, just zero
%all the negative numberss

dfofaverageNoZero=dfofaverage;
% dfofaverageNoZero(dfofaverageNoZero<0)=0;
dfofaverageNoZero=dfofaverageNoZero-min(dfofaverageNoZero);
mInField=mean(dfofaverageNoZero(PInFieldNew==1));

%number of bins for in field
I=length(find(PInFieldNew));

subplot(414)
plot(binCenters,dfofaveragesmoothFields,'black','LineWidth',2);
FFieldmax=max(dfofaveragesmoothFields);
axis([0 trackLength 0 FFieldmax*1.1]);
title('Fields Only');

%plot in field in patch colors
%make x y numbers for in-field patches
xInField=zeros(1,mN*4);
yInField=zeros(1,mN*4);
for n=1:mN;
    xInField(1,(1+(n-1)*4))=(allInFieldN(n,1)-1)*binWidth;
    xInField(1,(2+(n-1)*4))=(allInFieldN(n,1)-1)*binWidth;
    xInField(1,(3+(n-1)*4))=allInFieldN(n,2)*binWidth;
    xInField(1,(4+(n-1)*4))=allInFieldN(n,2)*binWidth;
    yInField(1,(1+(n-1)*4))=0;
    yInField(1,(2+(n-1)*4))=max(P1);
    yInField(1,(3+(n-1)*4))=max(P1);
    yInField(1,(4+(n-1)*4))=0;
end

%plot P value and patches
subplot(413)

patch(xInField,yInField,[ 1 0 0]);
hold on

%plot this one at last is to bring these curves to the front
plot(binCenters,P1,'black','LineWidth',2);
axis([0 trackLength 0 max(P1)*1.1]);
title('P value and in- and out-fields');
%plot threshold lines:
plot(binCenters,Pvalue1,'g','LineWidth',1);
plot(binCenters,Pvalue2,'g','LineWidth',1);
plot(binCenters,Pvalue3,'g','LineWidth',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% now test out-of-field

%if there is no out fields
if isempty(find(POutField,1));
    fields.dfofUse=Fuseallcell(:,cellNumber);
    fields.dfofaveragesmooth=dfofaverage;
    fields.dfofaveragesmoothFields=dfofaveragesmoothFields;
    fields.Pvalues=P1;
    fields.inFieldBins=inFieldBins;
    fields.outFieldBins=NaN;
    fields.transitions=NaN;
    fields.percentageBins=I/(trackLength/binWidth);
    fields.fieldStartEnds=fieldStartEnd;
    fields.fieldCenters=fieldCenters;
    fields.fieldSpacings=diff(fieldCenters);
    %if there is only one field, so there is no field spacing, set it to 0
    if length(fieldCenters)>1;
        fields.fieldSpacings=diff(fieldCenters);
    else
        fields.fieldSpacings=0;
    end
    fields.fieldWidths=fieldWidths;
    fields.meanInField=mInField;
    fields.meanOutField=NaN;
    fields.ratioInOut=NaN;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if there are out-of-field, test this criteria: the an out-of-field should have
% more than one bins
% like did previously for in-fields, to avoid the function 'contiguous'
% generates error, first test whether there are any continous 1s in the
% out-of-field vector

[rowOut,colOut]=find(POutField);
Dout=diff(rowOut);
%if there are continuous numbers, their must be 1 in D. the intersect below
%will tell whether there is any number in D that contains 1
if (isempty(intersect(1,Dout)));
    fields.dfofUse=Fuseallcell(:,cellNumber);
    fields.dfofaveragesmooth=dfofaverage;
    fields.dfofaveragesmoothFields=dfofaveragesmoothFields;
    fields.Pvalues=P1;
    fields.inFieldBins=inFieldBins;
    fields.outFieldBins=NaN;
    fields.transitions=NaN;
    fields.percentageBins=I/(trackLength/binWidth);
    fields.fieldStartEnds=fieldStartEnd;
    fields.fieldCenters=fieldCenters;
    fields.fieldSpacings=diff(fieldCenters);
    %if there is only one field, so there is no field spacing, set it to 0
    if length(fieldCenters)>1;
        fields.fieldSpacings=diff(fieldCenters);
    else
        fields.fieldSpacings=0;
    end
    fields.fieldWidths=fieldWidths;
    fields.meanInField=mInField;
    fields.meanOutField=NaN;
    fields.ratioInOut=NaN;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%next, if the previous step passed, check for a given out-field, there are
%more than one bins

runsOut=contiguous(POutField,[1]);

allOutField=runsOut{1,2};
[mO,nO] = size(allOutField);
fieldStartO=allOutField(:,1);
fieldEndO=allOutField(:,2);
for aO=1:mO;
    %below makes sure that for others, there will be at least 2 bins under 5cm bin, at least 4 bins under 2.5cm bin (maks this number to be adaptable to binWidth by equation).
    if fieldEndO(aO)-fieldStartO(aO)>=2*5/binWidth-1;
        POutField(fieldStartO(aO):fieldEndO(aO))=1;
    else
        POutField(fieldStartO(aO):fieldEndO(aO))=0;
    end
end

%%%%%%%%%%%

%next test whether now the POutField is empty
%if there is no out fields
if isempty(find(POutField,1));
    fields.dfofUse=Fuseallcell(:,cellNumber);
    fields.dfofaveragesmooth=dfofaverage;
    fields.dfofaveragesmoothFields=dfofaveragesmoothFields;
    fields.Pvalues=P1;
    fields.inFieldBins=inFieldBins;
    fields.outFieldBins=NaN;
    fields.transitions=NaN;
    fields.percentageBins=I/(trackLength/binWidth);
    fields.fieldStartEnds=fieldStartEnd;
    fields.fieldCenters=fieldCenters;
    fields.fieldSpacings=diff(fieldCenters);
    %if there is only one field, so there is no field spacing, set it to 0
    if length(fieldCenters)>1;
        fields.fieldSpacings=diff(fieldCenters);
    else
        fields.fieldSpacings=0;
    end
    fields.fieldWidths=fieldWidths;
    fields.meanInField=mInField;
    fields.meanOutField=NaN;
    fields.ratioInOut=NaN;
    return
end
%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%now all the in and out field bins were defined. calculate all the necessary
%perameters for in-fields




%calculate how many bins were assigned to in and out fields
%number of bins for in field
I=length(find(PInFieldNew));
%number of bins for out field
O=length(find(POutField));
%total number of in and out field bins
totalBins=I+O;

%calculate the ratio of mean in field firing ratio and mean out field
%firing rate

mOutField=mean(dfofaverageNoZero(POutField==1));
r=mInField/mOutField;

%calculate number of transitions

%indeces of out fields
outFieldBins=find(POutField);
%calculate number of transitions
[NTrans]=transitions(inFieldBins,outFieldBins);




%make x y numbers for in-field patches
%find continuous 1 in POutField
runsM=contiguous(POutField,[1]);
allOutField=runsM{1,2};
[mM,nM] = size(allOutField);
outfieldStartM=allOutField(:,1);
outfieldEndM=allOutField(:,2);

%plot out field in patch colors
%set patch parameters
xOutField=zeros(1,mM*4);
yOutField=zeros(1,mM*4);
for n=1:mM;
    xOutField(1,(1+(n-1)*4))=(outfieldStartM(n)-1)*binWidth;
    xOutField(1,(2+(n-1)*4))=(outfieldStartM(n)-1)*binWidth;
    xOutField(1,(3+(n-1)*4))=outfieldEndM(n)*binWidth;
    xOutField(1,(4+(n-1)*4))=outfieldEndM(n)*binWidth;
    yOutField(1,(1+(n-1)*4))=0;
    yOutField(1,(2+(n-1)*4))=max(P1);
    yOutField(1,(3+(n-1)*4))=max(P1);
    yOutField(1,(4+(n-1)*4))=0;
end


%this figure below is to replace the previous subplot 313 if the code can
%run up to here. without doing this, the outField batches will be on top of
%the P value so that the P value in the out field are not visible

%plot P value and patches
subplot(413)
patch(xInField,yInField,[ 1 0 0]);
hold on
patch(xOutField,yOutField,[0.3 0.3 0.3]);
hold on
%plot this one at last is to bring these curves to the front
plot(binCenters,P1,'black','LineWidth',2);
axis([0 trackLength 0 max(P1)*1.1]);
title('P value and in- and out-fields');
%plot threshold lines:
hold on
line([binCenters(1) binCenters(end)],[Pvalue1 Pvalue1],'Color','g','LineWidth',1);
hold on
line([binCenters(1) binCenters(end)],[Pvalue2 Pvalue2],'Color','g','LineWidth',1);
hold on
line([binCenters(1) binCenters(end)],[Pvalue3 Pvalue3],'Color','g','LineWidth',1);

%store the used dfof
fields.dfofUse=Fuseallcell(:,cellNumber);
%bined dfof. This smoothed dfofaverage
fields.dfofaveragesmooth=dfofaverage;%row vector
%field only dfof (smoothed)
fields.dfofaveragesmoothFields=dfofaveragesmoothFields;%row vector
%p value of all bins
fields.Pvalues=P1;
%indeces of in fields
fields.inFieldBins=inFieldBins;
%indeces of out fields
fields.outFieldBins=outFieldBins;
%number of transitions
fields.transitions=NTrans;
%percentage of bins assigned as either in and out fields
fields.percentageBins=totalBins/(trackLength/binWidth);
%field start and end locations (in cm). cell structure, each cell is one
%field, then it is one column vector, first raw is field start, second
%raw is field end.
fields.fieldStartEnds=fieldStartEnd;
%field centers (in cm)
fields.fieldCenters=fieldCenters;

%field spacings (in cm)
if length(fieldCenters)>1;
    fields.fieldSpacings=diff(fieldCenters);
else
    fields.fieldSpacings=0;
end

%field widths (in cm)
fields.fieldWidths=fieldWidths;
%mean binned dfof of all in fields
fields.meanInField=mInField;
%mean binned dfof of all out fields
fields.meanOutField=mOutField;
%ratio of mean binned dfof of in and out fields
fields.ratioInOut=r;

end

