function [allCellsCorrected]=correctPValue_clean(allCells,binWidth,trackLength)

%this code is to correct p value because sometimes P value identify
%multiple fields as fewer fields because of the slow decay of calcium.
%This code uses KY method to identify field based on the peak
%height and if in p value there are fused fields, they will be separated
%based on the peak height method.

%input:
%allCells, grid, nongrid,params: this is the output of p value;
%binWidth: unit is cm. this is binWidth for P value, normally it is 2.5cm
%for current p value. for previous old p value it was 5cm
%trackLength: give in cm

%these are things do not change
allCellsCorrected.dfofUse=allCells.dfofUse;
allCellsCorrected.dfofaveragesmooth=allCells.dfofaveragesmooth;
allCellsCorrected.dfofaveragesmoothFields=allCells.dfofaveragesmoothFields;%this will not change because no infield bins will be removed and no infield bins will be added
allCellsCorrected.Pvalues=allCells.Pvalues;
allCellsCorrected.inFieldBins=allCells.inFieldBins;
allCellsCorrected.outFieldBins=allCells.outFieldBins;
allCellsCorrected.transitions=allCells.transitions; %this will not change number of transitions because it doesn't add out fields betwen fields
allCellsCorrected.percentageBins=allCells.percentageBins;
allCellsCorrected.meanInField=allCells.meanInField;
allCellsCorrected.meanOutField=allCells.meanOutField;
allCellsCorrected.ratioInOut=allCells.ratioInOut;
allCellsCorrected.indices=allCells.indices;

%add binCenters, and grid, non-grid indices to it
allCellsCorrected.binCenters=[binWidth*0.5:binWidth:trackLength-0.5*binWidth];

%these are things may change
allCellsCorrected.fieldIndices={};%indices of all fields
allCellsCorrected.fieldActivities={};%activity of each field
allCellsCorrected.fieldStartEnds={};
allCellsCorrected.fieldCenters={};%this is center of mask
allCellsCorrected.fieldCenters_Peak={};%this is center of the peak acticity in each field
allCellsCorrected.fieldSpacings={};%this is determined by COM of field
allCellsCorrected.fieldSpacings_Peak={};%this is field spacings determined using peak locations
allCellsCorrected.fieldWidths={};
allCellsCorrected.minSpacing=[];%this is determined by COM of field
allCellsCorrected.minSpacing_Peak=[];%this is field spacings determined using peak locations
allCellsCorrected.maxWidth=[];
allCellsCorrected.KYFieldInfo={};

%%

for nCell=1:size(allCells.dfofaveragesmooth,2);
    if allCells.Pvalues(:,nCell)==0;
        
        allCellsCorrected.fieldStartEnds{nCell}=allCells.fieldStartEnds{nCell};
        allCellsCorrected.fieldCenters{nCell}=allCells.fieldCenters{nCell};%this is center of mask
        if ~isempty(allCells.fieldSpacings);
            allCellsCorrected.fieldSpacings{nCell}=allCells.fieldSpacings{nCell};%this is determined by COM of field
        end
        allCellsCorrected.fieldWidths{nCell}=allCells.fieldWidths{nCell};
        allCellsCorrected.minSpacing(nCell)=allCells.minSpacing(nCell);%this is determined by COM of field
        allCellsCorrected.maxWidth(nCell)=allCells.maxWidth(nCell);
        allCellsCorrected.fieldIndices{nCell}=NaN;
        allCellsCorrected.fieldActivities{nCell}={};
        allCellsCorrected.fieldCenters_Peak{nCell}=NaN;
        allCellsCorrected.fieldSpacings_Peak{nCell}=NaN;
        allCellsCorrected.minSpacing_Peak(nCell)=NaN;
        
    elseif allCells.Pvalues(:,nCell)==1;
        allCellsCorrected.fieldStartEnds{nCell}=allCells.fieldStartEnds{nCell};
        allCellsCorrected.fieldCenters{nCell}=allCells.fieldCenters{nCell};%this is center of mask
        if ~isempty(allCells.fieldSpacings);
            allCellsCorrected.fieldSpacings{nCell}=allCells.fieldSpacings{nCell};%this is determined by COM of field
        end
        allCellsCorrected.fieldWidths{nCell}=allCells.fieldWidths{nCell};
        allCellsCorrected.minSpacing(nCell)=allCells.minSpacing(nCell);%this is determined by COM of field
        allCellsCorrected.maxWidth(nCell)=allCells.maxWidth(nCell);
        allCellsCorrected.fieldIndices{nCell}=NaN;
        allCellsCorrected.fieldActivities{nCell}={};
        allCellsCorrected.fieldCenters_Peak{nCell}=NaN;
        allCellsCorrected.fieldSpacings_Peak{nCell}=NaN;
        allCellsCorrected.minSpacing_Peak(nCell)=NaN;
    else
        %first get fieldInfo using KY method
        [allCellsCorrected.KYFieldInfo{nCell}] = extractFrateRandom_GetFieldInfo(allCells.dfofaveragesmooth(:,nCell),binWidth);
        %this is the idx of all individual fields
        fieldInfo=allCellsCorrected.KYFieldInfo{nCell};
        b=fieldInfo.wholeFieldIdx;
        
        %next get the idx of all P value fields
        nBin=trackLength/binWidth;
        a=zeros(nBin,1);
        
        if isnan(allCells.inFieldBins{nCell});
            allCellsCorrected.fieldStartEnds{nCell}=allCells.fieldStartEnds{nCell};
            allCellsCorrected.fieldCenters{nCell}=allCells.fieldCenters{nCell};%this is center of mask
            if ~isempty(allCells.fieldSpacings);
                allCellsCorrected.fieldSpacings{nCell}=allCells.fieldSpacings{nCell};%this is determined by COM of field
            end
            allCellsCorrected.fieldWidths{nCell}=allCells.fieldWidths{nCell};
            allCellsCorrected.minSpacing(nCell)=allCells.minSpacing(nCell);%this is determined by COM of field
            allCellsCorrected.maxWidth(nCell)=allCells.maxWidth(nCell);
            allCellsCorrected.fieldIndices{nCell}=NaN;
            allCellsCorrected.fieldActivities{nCell}={};
            allCellsCorrected.fieldCenters_Peak{nCell}=NaN;
            allCellsCorrected.fieldSpacings_Peak{nCell}=NaN;
            allCellsCorrected.minSpacing_Peak(nCell)=NaN;
            
        else
            
            
            a(allCells.inFieldBins{nCell})=1;
            p=contiguous(a,[1]);
            aa=p{1,2};
            c={};
            for n=1:size(aa,1);
                c{n}=aa(n,1):1:aa(n,2);
            end
            
            %check which field has fused multiple fields
            k=zeros(length(c),1);
            koverlap={};
            for n=1:length(c);
                t=c{n};
                koverlap{n}=[];
                for m=1:length(b);
                    i=intersect(c{n},b{m});
                    if isempty(i);
                        k(n);
                        koverlap{n};
                    else
                        k(n)=k(n)+1;
                        koverlap{n}(end+1)=m;
                    end
                end
            end
            
            if isempty(find(k>1));
                allCellsCorrected.fieldStartEnds{nCell}=allCells.fieldStartEnds{nCell};
                allCellsCorrected.fieldCenters{nCell}=allCells.fieldCenters{nCell};%this is center of mask
                if ~isempty(allCells.fieldSpacings);
                    allCellsCorrected.fieldSpacings{nCell}=allCells.fieldSpacings{nCell};%this is determined by COM of field
                end
                allCellsCorrected.fieldWidths{nCell}=allCells.fieldWidths{nCell};
                allCellsCorrected.minSpacing(nCell)=allCells.minSpacing(nCell);%this is determined by COM of field
                allCellsCorrected.maxWidth(nCell)=allCells.maxWidth(nCell);
                allCellsCorrected.fieldIndices{nCell}=c;
                allCellsCorrected.fieldActivities{nCell}={};
                for nn=1:length(c);
                    allCellsCorrected.fieldActivities{nCell}{nn}=allCells.dfofaveragesmooth(c{nn},nCell);
                    [~,hh]=max(allCellsCorrected.fieldActivities{nCell}{nn});
                    allCellsCorrected.fieldCenters_Peak{nCell}(nn)=c{nn}(hh)*binWidth-binWidth*0.5;
                end
                
                if length(c)>1;
                    allCellsCorrected.fieldSpacings_Peak{nCell}=diff(allCellsCorrected.fieldCenters_Peak{nCell});
                else
                    allCellsCorrected.fieldSpacings_Peak{nCell}=0;
                end
                allCellsCorrected.minSpacing_Peak(nCell)=min(allCellsCorrected.fieldSpacings_Peak{nCell});
                
            else
                ck=c(find(k>1));
                ckoverlap=koverlap(find(k>1));
                newFieldIdx={};
                for n=1:length(ck);
                    oldIdx=ck{n};
                    newFieldEnd=[];
                    for m=1:length(ckoverlap{n});
                        thisNewIdx=b{ckoverlap{n}(m)};
                        if m==1;
                            newFieldIdx{end+1}=min(oldIdx):1:max(thisNewIdx);
                            newFieldEnd(m,1)=max(thisNewIdx);
                        elseif m==length(ckoverlap{n});
                            newFieldIdx{end+1}=min(thisNewIdx):1:max(oldIdx);
                            newFieldEnd(m,1)=max(oldIdx);
                        else
                            newFieldIdx{end+1}=intersect(oldIdx,thisNewIdx);
                            newFieldEnd(m,1)=max(intersect(oldIdx,thisNewIdx));
                        end
                    end
                    
                    %if there are (just in case) indices that are not assigned to any field
                    %(they can only be the indices in the middle), just put it to its close
                    %field
                    
                    newFieldIdxOfThisOldField=cell2mat(newFieldIdx(end-length(ckoverlap{n})+1:end));
                    
                    df=setdiff(oldIdx,newFieldIdxOfThisOldField);
                    if ~isempty(df);
                        for h=1:length(df);
                            v=abs(newFieldEnd-df(h));
                            [~,ii]=min(v);
                            %put it to the field that it is more close to
                            newFieldIdx{end-length(ckoverlap{n})+ii}=unique(sort([newFieldIdx{end-length(ckoverlap{n})+ii} df(h)]));
                        end
                    end
                    
                end
                
                
                allNewFieldIdx=[c(find(k<=1)) newFieldIdx];
                r=[];
                for n=1:length(allNewFieldIdx);
                    r(n)=allNewFieldIdx{n}(1);
                end
                
                [~,i]=sort(r);
                allNewFieldIdx=allNewFieldIdx(i);
                
                fieldStartEnds=[];
                fieldWidths=[];
                fieldActivities={};
                %this is based on center of mask
                fieldCenters=[];
                fieldCenters_Peak=[];
                
                for n=1:length(allNewFieldIdx);
                    fieldStartEnds(n,1)=(allNewFieldIdx{n}(1)-1)*binWidth;
                    fieldStartEnds(n,2)=allNewFieldIdx{n}(end)*binWidth;
                    fieldWidths(n,1)=fieldStartEnds(n,2)-fieldStartEnds(n,1);
                    fieldCenters(n)=mean(fieldStartEnds(n,:));
                    fieldActivities{n}=allCells.dfofaveragesmooth(allNewFieldIdx{n},nCell);
                    [~,hh]=max(fieldActivities{n});
                    fieldCenters_Peak(n)=allNewFieldIdx{n}(hh)*binWidth-binWidth*0.5;
                end
                
                allCellsCorrected.fieldStartEnds{nCell}=fieldStartEnds;
                allCellsCorrected.fieldCenters{nCell}=fieldCenters;%this is center of mask
                allCellsCorrected.fieldSpacings{nCell}=diff(fieldCenters);%this is determined by COM of field
                allCellsCorrected.fieldWidths{nCell}=fieldWidths;
                allCellsCorrected.minSpacing(nCell)=min(diff(fieldCenters));%this is determined by COM of field
                allCellsCorrected.maxWidth(nCell)=max(fieldWidths);
                allCellsCorrected.fieldIndices{nCell}=allNewFieldIdx;
                allCellsCorrected.fieldActivities{nCell}=fieldActivities;
                allCellsCorrected.fieldCenters_Peak{nCell}=fieldCenters_Peak;
                allCellsCorrected.fieldSpacings_Peak{nCell}=diff(fieldCenters_Peak);
                allCellsCorrected.minSpacing_Peak(nCell)=min(diff(fieldCenters_Peak));
                
                
            end
            
        end
        
    end
end

save('allCellsCorrected.mat','allCellsCorrected');

end

