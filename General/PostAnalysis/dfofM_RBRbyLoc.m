function dfofM_RBRbyLoc(p)
%% dfofM_RBRbyLoc
% calculates run-by-run consistency on interpolated activity as a
% function of track location
    
p1 = pwd;

if nargin>0
    cd(p)
end

%%

sfxs = {'_dfof','_sig'};
sfxsSec = {'','_sig'};
% number of bins in rolling average
rollN = 5;

for ii = 1:2
    % load mean activity
    load(['dfofaveragesmooth' sfxsSec{ii} '_interp.mat'])

    cd(['RunByRun' sfxs{ii}])

    % load run-by-run activity
    load(['dfofMInterpM' sfxsSec{ii} '.mat']);
    load('corrInfo.mat');

    if ii==1
        dfofM = dfofMInterpM;
        dfofSmooth = dfofaveragesmooth_interp;
    else
        dfofM = dfofMInterpM_sig;
        dfofSmooth = dfofaveragesmooth_sig_interp;
    end

    useRun = [];
    jj = 1;
    while isempty(useRun) & jj<=size(corrInfo.noNaN,2)
        useRun = find(corrInfo.noNaN(:,jj));
        jj = jj+1;
    end

    % use at most the first 20 runs
    if length(useRun)>20
        useRun=useRun(1:20);
    end

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

        A = dfofM{n}(useRun,:);

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
        end

        % calculate mean RBR distribution by cell
        corrInfoLocation.toMeanMean(n,:) = mean(corrInfoLocation.toMean{n},1,'omitnan');
        corrInfoLocation.toOthersMean(n,:) = mean(corrInfoLocation.toOthers{n},1,'omitnan');
        corrInfoLocation.toNextMean(n,:) = mean(corrInfoLocation.toNext{n},1,'omitnan');

    end

    % calculate mean RBR distribution
    corrInfoLocation.toMeanMeanMean = mean(corrInfoLocation.toMeanMean,1,'omitnan');
    corrInfoLocation.toMeanMeanSEM = nansem(corrInfoLocation.toMeanMean,1);

    corrInfoLocation.toOthersMeanMean = mean(corrInfoLocation.toOthersMean,1,'omitnan');
    corrInfoLocation.toOthersMeanSEM = nansem(corrInfoLocation.toOthersMean,1);

    corrInfoLocation.toNextMeanMean = mean(corrInfoLocation.toNextMean,1,'omitnan');
    corrInfoLocation.toNextMeanSEM = nansem(corrInfoLocation.toNextMean,1);

    save('corrInfoLocation.mat','corrInfoLocation');

    cd ..
end

% return to start directory
cd(p1)
