function dfofMCorrelation_clean(trackEnd,binWidth)

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

[speedThreshold]= speedThreshold1D( 1,length(abfFake.t),abfFake );
close

trackStart=0;

binName = num2str(binWidth);
binName = replace(binName,'.','-');


%%

disp('correlation RunByRun dfof')

clear dfofM
%get run by run
mkdir('RunByRun_dfof');
cd('RunByRun_dfof');
[dfofM]=getRunByRunActivity(nCellStart,nCellEnd,imageStartNumber,imageEndNumber,abfStartNumber,abfEndNumber,binWidth,trackStart,trackEnd,speedThreshold,dfof,abfFake);
save(['dfofM_' binName '.mat'],'dfofM');

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

save(['dfofMInterpM_' binName '.mat'],'dfofMInterpM');

cd ..

%%
disp('correlation RunByRun sig')
clear dfofM_sig
%get run by run
mkdir('RunByRun_sig');
cd('RunByRun_sig');
[dfofM_sig]=getRunByRunActivity(nCellStart,nCellEnd,imageStartNumber,imageEndNumber,abfStartNumber,abfEndNumber,binWidth,trackStart,trackEnd,speedThreshold,dfof_sig,abfFake);
save(['dfofM_sig_' binName '.mat'],'dfofM_sig');


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

save(['dfofMInterpM_sig_' binName '.mat'],'dfofMInterpM_sig');


cd ../
