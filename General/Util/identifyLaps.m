function [runIdx,nRuns] = identifyLaps(y,startLap,stopLap,type)
%% identifyLaps
% Identifies the lap indices of a Virmen or abf file. Outputs a matrix
% containing the start and stop indices (columns) for the desried lap range
% (rows). Type should equal 'virmen' (default) or 'abf'. Inputs 1-3 are
% currently only used for 'virmen' analysis.
%

if nargin<2 || isempty(startLap)
    startLap = 1;
end

if nargin<3 || isempty(stopLap)
    stopLap = 0;
end

if nargin<4 || isempty(type)
    type = 'virmen';
end

switch type
    case 'virmen'

        % find velocity
        C = diff(y);
        C = [C;0];

        % Checks for full or partial teleportations points
        A=C<=-max(y)*0.2;

        % find contiguous teleportation points
        N=contiguous(A,1);
        NN=N{1,2};
        NN=[[0 0]; NN; [length(y)+1 length(y)+1]];

        if stopLap==0 % use all runs
            useNN = NN(startLap:end,:);
        else % set laps (last lap is always excluded)
            useNN = NN(startLap:min(stopLap+1,size(NN,1)-1),:);
        end

        nRuns = size(useNN,1)-1;

        runIdx = zeros(nRuns,2);
        for m=1:nRuns
            runIdx(m,1)=useNN(m,2)+1;
            runIdx(m,2)=useNN(m+1,1)-1;
        end

        if runIdx(end,1)==runIdx(end,2)
            runIdx(end,:) = [];
            nRuns = nRuns-1;
        end

    case 'abf'
        % load behavior sync
        load('md.mat','md')
        laps = md(:,4)*100;
        y = md(:,3)*100;

        % find lap ends
        N = contiguous(round(laps));
        NN = cat(1,N{:,2});

        nRuns = size(NN,1);
        for ii = 1:nRuns
            C = diff(y(NN(ii,1):NN(ii,2)));

            % remove teleportation at start of lap
            if C(1)<=-max(y)*0.2
                idxFirst = find(C>-max(y)*0.2,1,'first');
                NN(ii,1) = NN(ii,1) + idxFirst - 1;
            end

            % remove teleportation at end of lap
            if C(end)<=-max(y)*0.2
                idxLast = find(C>-max(y)*0.2,1,'last');
                NN(ii,2) = NN(ii,2)+idxLast-diff(NN(ii,:));
            end
        end

        % save output
        runIdx = NN;

    case 'abf_noLoad'
        % identify teleportation points
        speed = diff(y)';
        idxTele = [find(speed<-100) length(y)];
        A = diff([1 idxTele]);

        % identify start and end indices
        idxStart = [1 idxTele(A==1)+1];
        idxEnd = idxTele(A>1);

        % correct for mismatched index numbers
        iS=length(idxStart);
        iE=length(idxEnd);
        if iS~=iE
            iSE = min([iS iE]);
            idxStart = idxStart(1:iSE);
            idxEnd = idxEnd(1:iSE);
        end

        % store outputs
        runIdx = [idxStart' idxEnd'];
        nRuns = size(runIdx,1);
end