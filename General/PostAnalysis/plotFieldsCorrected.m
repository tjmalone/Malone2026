clear; close all; clc

load('abf.mat','abf')
load('gridAnalysis_sig\allCells.mat')
load('gridAnalysis_sig\allCellsCorrected.mat')
load('gridAnalysis_sig\params.mat')

d = dir('dfof_sig*');
load(d(1).name)
dfof = dfof_sig;

binWidth = 5;
trackLength = 600;

useCell = 96;

savedir = 'gridAnalysis_sig\allCorrected';
mkdir(savedir);

%compute use run number
abfStartNumber=1;
abfEndNumber=length(abf.y);
abfIms = abfStartNumber:abfEndNumber;
% speed at each index of abfIms
speed=diff(abf.y(abf.imageIndex(abfIms)));

%this is to find the telepore point where the speed is a big negative
[row,col]=find(speed<-100);
endIndicesP=[row' abfEndNumber]';
startIndicesP=[1 (row+1)']';

startIndices=startIndicesP(diff(startIndicesP)>1);
startIndices=[startIndices' startIndicesP(end)]';
%do the problem happens to the endIndices, solve it like this
endIndices=endIndicesP(diff([1 endIndicesP'])>1);

for nCell = useCell

    filename1 = sprintf('%s_%d.fig', 'cell', nCell);
    filename2 = sprintf('%s_%d.tif', 'cell', nCell);

    if allCells.Pvalues(:,nCell)==0;
        figure,
        saveas(gcf,fullfile(savedir,filename1));
    else

        useRunNumbers=params{nCell}.useRunNumbers;
        useIdx=[];
        for iuse=1:length(useRunNumbers);
            useIdx=[useIdx,startIndices(useRunNumbers(iuse)):1:endIndices(useRunNumbers(iuse))];
        end


        abfFakeNew.t=abf.t(useIdx);
        abfFakeNew.y=abf.y(useIdx);
        abfFakeNew.imageIndex=1:1:length(useIdx);
        usedfof=dfof(useIdx,:);
        abfEndNumber=length(useIdx);
        imageEndNumber=abfEndNumber;
        abfStartNumber=1;
        imageStartNumber=1;


        figure,
        subplot(411);
        plotdfofrunbyrun(nCell,imageStartNumber,imageEndNumber,abfStartNumber,abfEndNumber,...
            params{nCell}.binWidth,params{nCell}.trackStart,params{nCell}.trackEnd,...
            params{nCell}.speedThreshold,usedfof,abfFakeNew,1);
        title(['cell',num2str(nCell),'   run by run']);
        subplot(412);
        plot(allCellsCorrected.binCenters,allCellsCorrected.dfofaveragesmooth(:,nCell),'black','LineWidth',2);
        Fmax=max(allCellsCorrected.dfofaveragesmooth(:,nCell));

        if ~isempty(find(allCells.Pvalues(:,nCell)~=1));

            axis([0 trackLength min(allCellsCorrected.dfofaveragesmooth(:,nCell)) Fmax*1.1]);
        end
        title('mean dfof');
        if isnan(allCellsCorrected.inFieldBins{nCell});
            saveas(gcf,fullfile(savedir,filename1));

            close
            continue
        end


        xInField=zeros(1,length(allCellsCorrected.fieldIndices{nCell})*4);
        yInField=zeros(1,length(allCellsCorrected.fieldIndices{nCell})*4);
        for n=1:length(allCellsCorrected.fieldIndices{nCell});
            xInField(1,(1+(n-1)*4))=(allCellsCorrected.fieldIndices{nCell}{n}(1)-1)*binWidth;
            xInField(1,(2+(n-1)*4))=(allCellsCorrected.fieldIndices{nCell}{n}(1)-1)*binWidth;
            xInField(1,(3+(n-1)*4))=allCellsCorrected.fieldIndices{nCell}{n}(end)*binWidth;
            xInField(1,(4+(n-1)*4))=allCellsCorrected.fieldIndices{nCell}{n}(end)*binWidth;
            yInField(1,(1+(n-1)*4))=0;
            yInField(1,(2+(n-1)*4))=max(allCellsCorrected.Pvalues(:,nCell));
            yInField(1,(3+(n-1)*4))=max(allCellsCorrected.Pvalues(:,nCell));
            yInField(1,(4+(n-1)*4))=0;
        end
        subplot(413)
        patch(xInField,yInField,[ 1 0 0]);
        hold on
        plot(allCellsCorrected.binCenters,allCellsCorrected.Pvalues(:,nCell),'black','LineWidth',2);
        axis([0 trackLength 0 max(allCellsCorrected.Pvalues(:,nCell))*1.1]);
        title('P value and in- and out-fields');
        %plot threshold lines:
        hold on
        line([allCellsCorrected.binCenters(1) allCellsCorrected.binCenters(end)],[params{nCell}.Pvalue1 params{nCell}.Pvalue1],'Color','g','LineWidth',1);
        hold on
        line([allCellsCorrected.binCenters(1) allCellsCorrected.binCenters(end)],[params{nCell}.Pvalue2 params{nCell}.Pvalue2],'Color','g','LineWidth',1);
        hold on
        line([allCellsCorrected.binCenters(1) allCellsCorrected.binCenters(end)],[params{nCell}.Pvalue3 params{nCell}.Pvalue3],'Color','g','LineWidth',1);

        subplot(414)
        plot(allCellsCorrected.binCenters,allCellsCorrected.dfofaveragesmoothFields(:,nCell),'black','LineWidth',2);
        FFieldmax=max(allCellsCorrected.dfofaveragesmoothFields(:,nCell));
        axis([0 trackLength 0 FFieldmax*1.1]);
        title('Fields Only');

        if isnan(allCellsCorrected.outFieldBins{nCell});
            saveas(gcf,fullfile(savedir,filename1));

            close
            continue
        end

        out=zeros(trackLength/binWidth,1);
        out(allCellsCorrected.outFieldBins{nCell})=1;
        runsM=contiguous(out,[1]);
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
            yOutField(1,(2+(n-1)*4))=max(allCellsCorrected.Pvalues(:,nCell));
            yOutField(1,(3+(n-1)*4))=max(allCellsCorrected.Pvalues(:,nCell));
            yOutField(1,(4+(n-1)*4))=0;
        end

        %plot P value and patches
        subplot(413)
        patch(xInField,yInField,[ 1 0 0]);
        hold on
        patch(xOutField,yOutField,[0.3 0.3 0.3]);
        hold on
        %plot this one at last is to bring these curves to the front
        plot(allCellsCorrected.binCenters,allCellsCorrected.Pvalues(:,nCell),'black','LineWidth',2);
        axis([0 trackLength 0 max(allCellsCorrected.Pvalues(:,nCell))*1.1]);
        title('P value and in- and out-fields');
        %plot threshold lines:
        hold on
        line([allCellsCorrected.binCenters(1) allCellsCorrected.binCenters(end)],[params{nCell}.Pvalue1 params{nCell}.Pvalue1],'Color','g','LineWidth',1);
        hold on
        line([allCellsCorrected.binCenters(1) allCellsCorrected.binCenters(end)],[params{nCell}.Pvalue2 params{nCell}.Pvalue2],'Color','g','LineWidth',1);
        hold on
        line([allCellsCorrected.binCenters(1) allCellsCorrected.binCenters(end)],[params{nCell}.Pvalue3 params{nCell}.Pvalue3],'Color','g','LineWidth',1);

        %save
        saveas(gcf,fullfile(savedir,filename1));
        close
    end

end
