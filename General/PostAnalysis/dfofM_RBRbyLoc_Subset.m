function corrInfoLocation = dfofM_RBRbyLoc_Subset(dfofM,dfofSmooth)
%% dfofM_RBRbyLoc
% calculates run-by-run consistency on interpolated activity as a
% function of track location

% number of bins in rolling average
rollN = 5;

corrInfoLocation=[];
corrInfoLocation.toMean={};
corrInfoLocation.toMeanMean=[];
corrInfoLocation.toMeanMeanMean=[];
corrInfoLocation.toMeanMeanSEM=[];

corrInfoLocation.toOthers={};
corrInfoLocation.toOthersMean=[];
corrInfoLocation.toOthersMeanMean=[];
corrInfoLocation.toOthersMeanSEM=[];

corrInfoLocation.toNext={};
corrInfoLocation.toNextMean=[];
corrInfoLocation.toNextMeanMean=[];
corrInfoLocation.toNextMeanSEM=[];

% calculate RBR dirstribution by cell
for n=1:length(dfofM)
    corrInfoLocation.toMean{n} = [];
    corrInfoLocation.toOthers{n} = [];
    corrInfoLocation.toNext{n} = [];
    corrInfoLocation.EMDtoOthers{n} = [];

    A = dfofM{n};

    % current mean activity
    meanCur = dfofSmooth(:,n);

    % perform RBR correlation for rolling average
    for m=1:size(A,2)-(rollN-1)
        curRange = m:m+rollN-1;

        % mean activity correlation
        corrMean = corr(A(:,curRange)',meanCur(curRange));
        corrInfoLocation.toMean{n}(:,m) = corrMean;

        % lap to lap correlations
        c = corr(A(:,curRange)',A(:,curRange)');

        % other run correlation
        corrOther  = c(tril(true(size(c)),-1));
        corrInfoLocation.toOthers{n}(:,m) = corrOther;

        % adjacent run correlation
        corrNext = diag(c,-1);
        corrInfoLocation.toNext{n}(:,m) = corrNext;

        % lap to lap EMD
        e = calculateEMD(A(:,curRange)',A(:,curRange)');

        % other run EMD
        emdOther  = e(tril(true(size(e)),-1));
        corrInfoLocation.EMDtoOthers{n}(:,m) = emdOther;
    end

    % calculate mean RBR distribution by cell
    corrInfoLocation.toMeanMean(n,:) = mean(corrInfoLocation.toMean{n},1,'omitnan');
    corrInfoLocation.toOthersMean(n,:) = mean(corrInfoLocation.toOthers{n},1,'omitnan');
    corrInfoLocation.toNextMean(n,:) = mean(corrInfoLocation.toNext{n},1,'omitnan');
    corrInfoLocation.EMDtoOthersMean(n,:) = mean(corrInfoLocation.EMDtoOthers{n},1,'omitnan');

end

% calculate mean RBR distribution
corrInfoLocation.toMeanMeanMean = mean(corrInfoLocation.toMeanMean,1,'omitnan');
corrInfoLocation.toMeanMeanSEM = nansem(corrInfoLocation.toMeanMean,1);

corrInfoLocation.toOthersMeanMean = mean(corrInfoLocation.toOthersMean,1,'omitnan');
corrInfoLocation.toOthersMeanSEM = nansem(corrInfoLocation.toOthersMean,1);

corrInfoLocation.toNextMeanMean = mean(corrInfoLocation.toNextMean,1,'omitnan');
corrInfoLocation.toNextMeanSEM = nansem(corrInfoLocation.toNextMean,1);

corrInfoLocation.EMDtoOthersMeanMean = mean(corrInfoLocation.EMDtoOthersMean,1,'omitnan');
corrInfoLocation.EMDtoOthersMeanSEM = nansem(corrInfoLocation.EMDtoOthersMean,1);

end

