function dfofMCorrelationFull(trackEnd,binWidth)

name='dfof_*';
fn= dir(name);
load(fn(1).name);
load(fn(2).name);

load('abfFake.mat');

N=size(dfof,2);%N is the number of cells that needs to get dfof
%load parameters for the slice fit
nCellStart=1;
nCellEnd=N;
imageStartNumber=1;
imageEndNumber=length(abfFake.t);
abfStartNumber=imageStartNumber;
abfEndNumber=imageEndNumber;

[speedThreshold]= speedThreshold1D( 1,length(abfFake.t),abfFake,0);
close

trackStart=0;


%%


name='dfofaveragesmooth_*';
fn= dir(name);
for n=1:length(fn);
    load(fn(n).name);
end

disp('interp averaged activity')
%generate interpolated dfofaveragesmooth and dfofaveragesmooth_sig
dfofaveragesmooth_interp=[];
dfofaveragesmooth_sig_interp=[];

for n=1:size(dfofaveragesmooth,2);
    disp(n)
    A=dfofaveragesmooth(:,n);
    if length(find(~isnan(A)))>1;%at least two numbers are not nan
        dfofaveragesmooth_interp(:,n)=naninterp(A);
    else
        A(isnan(A))=0;
        dfofaveragesmooth_interp(:,n)=A;
    end

    A=dfofaveragesmooth_sig(:,n);
    if length(find(~isnan(A)))>1;%at least two numbers are not nan
        dfofaveragesmooth_sig_interp(:,n)=naninterp(A);
    else
        %if there are less than 2 non-nan numbers, keep nan
        dfofaveragesmooth_sig_interp(:,n)=A;
    end

end

save('dfofaveragesmooth_interp.mat','dfofaveragesmooth_interp');
save('dfofaveragesmooth_sig_interp.mat','dfofaveragesmooth_sig_interp');

%%

disp('correlation RunByRun dfof')

clear dfofM
%get run by run
mkdir('RunByRun_dfof');
cd('RunByRun_dfof');
[dfofM]=getRunByRunActivity(nCellStart,nCellEnd,imageStartNumber,imageEndNumber,abfStartNumber,abfEndNumber,binWidth,trackStart,trackEnd,speedThreshold,dfof,abfFake);
save('dfofM.mat','dfofM');

r=[];%generally now many R,how many C
c=[];
%find a general size of dfofM
if length(dfofM)>=5;
    for n=1:5;
        if ~isempty(dfofM{n});
            r=size(dfofM{n},1);
            c=size(dfofM{n},2);
        end
    end
else
    for n=1:length(dfofM);
        if ~isempty(dfofM{n});
            r=size(dfofM{n},1);
            c=size(dfofM{n},2);
        end
    end
end

dfofMInterpM={};%remove first and last runs
dfofEnds = cellfun(@(x) x(end,end-1:end),dfofM,'UniformOutput',0);
dfofEmptys = cellfun(@isempty,dfofM);
dfofEnds = cat(2,dfofEnds{dfofEmptys==0});
nanEnd = all(isnan(dfofEnds));
for n=1:length(dfofM);
    dfofMInterpM{n}=[];
    thisdfofM=dfofM{n};
    if isempty(thisdfofM);
        if nanEnd
            dfofMInterpM{n}=NaN(r-2,c);
        else
            dfofMInterpM{n}=NaN(r-1,c);
        end
    else
        for m=2:size(thisdfofM,1);
            % include last lap if one of last two bins is non-nan (TM 2/22/24)
            if m==size(thisdfofM,1) && nanEnd
                continue
            end
            if length(find(~isnan(thisdfofM(m,:))))>1;%at least two numbers are not nan
                dfofMInterpM{n}(m-1,:)=naninterp(thisdfofM(m,:));
                %if there are less than 2 non-nan numbers, keep nan
            else
                dfofMInterpM{n}(m-1,:)=thisdfofM(m,:);
            end
        end
    end
end

save('dfofMInterpM.mat','dfofMInterpM');

corrInfo=[];
%first compute the correlation of each run to the mean dfofaveragesmooth

corrInfo.noNaN=[];%indicate runs with no nans. If a run has no nan, it will be 1, if it has nan, it will be 0.
for n=1:length(dfofMInterpM);
    for m=1:size(dfofMInterpM{n},1);
        if ~isempty(find(isnan(dfofMInterpM{n}(m,:))));
            corrInfo.noNaN(m,n)=0;
        else
            corrInfo.noNaN(m,n)=1;
        end
    end
end

%first compute the correlation of each run to the mean dfofaveragesmooth
corrInfo.toMean=[];%first compute the correlation of each run to the mean dfofaveragesmooth
corrInfo.meantoMean=[]; %mean of to mean
for n=1:length(dfofMInterpM);
    for m=1:size(dfofMInterpM{n},1);
        if corrInfo.noNaN(m,n)==1
            corrInfo.toMean(m,n)=corr(dfofMInterpM{n}(m,:)',dfofaveragesmooth_interp(:,n));
        else
            corrInfo.toMean(m,n)=NaN;
        end
    end
    A=corrInfo.toMean(:,n);
    corrInfo.meantoMean(1,n)=mean(A(~isnan(A)));
end

%second compute the adjescent run correlation
corrInfo.toNext=[];
corrInfo.meantoNext=[];
for n=1:length(dfofMInterpM);
    for m=1:size(dfofMInterpM{n},1)-1;
        A=dfofMInterpM{n}(m,:);
        B=dfofMInterpM{n}(m+1,:);
        if corrInfo.noNaN(m,n)*corrInfo.noNaN(m+1,n)==1;
            corrInfo.toNext(m,n)=corr(A',B');
        else
            corrInfo.toNext(m,n)=NaN;
        end
    end
    C=corrInfo.toNext(:,n);
    corrInfo.meantoNext(1,n)=mean(C(~isnan(C)));
end

%3rd compute the corr of this run and all the others after it
corrInfo.toOthers=[];
corrInfo.meantoOthers=[];
for n=1:length(dfofMInterpM);
    H=[];
    for m=1:size(dfofMInterpM{n},1)-1;
        for i=m+1:size(dfofMInterpM{n},1);
            A=dfofMInterpM{n}(m,:);
            B=dfofMInterpM{n}(i,:);
            if corrInfo.noNaN(m,n)*corrInfo.noNaN(i,n)==1;
                H(end+1,1)=corr(A',B');
            else
                H(end+1,1)=NaN;
            end
        end
    end
    corrInfo.toOthers(:,n)=H;
    corrInfo.meantoOthers(1,n)=mean(H(~isnan(H)));
end

save('corrInfo.mat','corrInfo');
cd ../


%%
disp('correlation RunByRun sig')
clear dfofM_sig
%get run by run
mkdir('RunByRun_sig');
cd('RunByRun_sig');
[dfofM_sig]=getRunByRunActivity(nCellStart,nCellEnd,imageStartNumber,imageEndNumber,abfStartNumber,abfEndNumber,binWidth,trackStart,trackEnd,speedThreshold,dfof_sig,abfFake);
save('dfofM_sig.mat','dfofM_sig');


dfofMInterpM_sig={};%remove first and last runs

dfofMInterpM_sig={};%remove first and last runs
for n=1:length(dfofM_sig);
    dfofMInterpM_sig{n}=[];
    thisdfofM=dfofM_sig{n};
    if isempty(thisdfofM);
        if nanEnd
            dfofMInterpM_sig{n}=NaN(r-2,c);
        else
            dfofMInterpM_sig{n}=NaN(r-1,c);
        end
    else
        for m=2:size(thisdfofM,1);
            % include last lap if one of last two bins is non-nan (TM 2/22/24)
            if m==size(thisdfofM,1) && nanEnd
                continue
            end
            if length(find(~isnan(thisdfofM(m,:))))>1;%at least two numbers are not nan
                dfofMInterpM_sig{n}(m-1,:)=naninterp(thisdfofM(m,:));
                %if there are less than 2 non-nan numbers, keep nan
            else
                dfofMInterpM_sig{n}(m-1,:)=thisdfofM(m,:);
            end
        end
    end
end

save('dfofMInterpM_sig.mat','dfofMInterpM_sig');

corrInfo=[];
%first compute the correlation of each run to the mean dfofaveragesmooth

corrInfo.noNaN=[];%indicate runs with no nans. If a run has no nan, it will be 1, if it has nan, it will be 0.
for n=1:length(dfofMInterpM_sig);
    for m=1:size(dfofMInterpM_sig{n},1);
        if ~isempty(find(isnan(dfofMInterpM_sig{n}(m,:))));
            corrInfo.noNaN(m,n)=0;
        else
            corrInfo.noNaN(m,n)=1;
        end
    end
end

%first compute the correlation of each run to the mean dfofaveragesmooth
corrInfo.toMean=[];%first compute the correlation of each run to the mean dfofaveragesmooth
corrInfo.meantoMean=[]; %mean of to mean
for n=1:length(dfofMInterpM_sig);
    for m=1:size(dfofMInterpM_sig{n},1);
        if corrInfo.noNaN(m,n)==1
            corrInfo.toMean(m,n)=corr(dfofMInterpM_sig{n}(m,:)',dfofaveragesmooth_sig_interp(:,n));
        else
            corrInfo.toMean(m,n)=NaN;
        end
    end
    A=corrInfo.toMean(:,n);
    corrInfo.meantoMean(1,n)=mean(A(~isnan(A)));
end

%second compute the adjescent run correlation
corrInfo.toNext=[];
corrInfo.meantoNext=[];
for n=1:length(dfofMInterpM_sig);
    for m=1:size(dfofMInterpM_sig{n},1)-1;
        A=dfofMInterpM_sig{n}(m,:);
        B=dfofMInterpM_sig{n}(m+1,:);
        if corrInfo.noNaN(m,n)*corrInfo.noNaN(m+1,n)==1;
            corrInfo.toNext(m,n)=corr(A',B');
        else
            corrInfo.toNext(m,n)=NaN;
        end
    end
    C=corrInfo.toNext(:,n);
    corrInfo.meantoNext(1,n)=mean(C(~isnan(C)));
end

%3rd compute the corr of this run and all the others after it
corrInfo.toOthers=[];
corrInfo.meantoOthers=[];
for n=1:length(dfofMInterpM_sig);
    H=[];
    for m=1:size(dfofMInterpM_sig{n},1)-1;
        for i=m+1:size(dfofMInterpM_sig{n},1);
            A=dfofMInterpM_sig{n}(m,:);
            B=dfofMInterpM_sig{n}(i,:);
            if corrInfo.noNaN(m,n)*corrInfo.noNaN(i,n)==1;
                H(end+1,1)=corr(A',B');
            else
                H(end+1,1)=NaN;
            end
        end
    end
    corrInfo.toOthers(:,n)=H;
    corrInfo.meantoOthers(1,n)=mean(H(~isnan(H)));
end

save('corrInfo.mat','corrInfo');
cd ../
