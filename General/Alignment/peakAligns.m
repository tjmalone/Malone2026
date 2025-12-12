function [alignsFinal,delRowsAll] = peakAligns(alignData,rsltFold,keepDupGroup)
%% peakAligns.m
% Finds the best alignments among a full set of alignments (including
% duplicates) based on alignment confidence as calculated in alignmentPost.
% If a duplicaiton groups contains multiple indpendent alignments, these
% groups will be identified by iteratively removing alignments matching the
% peak alignments.
%
% Inputs: 
%   alignData - combined data matrix contains the following:
%       dupGroup - defines cells in common duplication groups
%       allCommon - defines all possible alignments
%       allConf - confidence levels for the corresponding alignment
%   keepDupGroup - whether to store orignal duplication group in
%       alignsFinal. Default (=0) converts all dupGroup to 0. Storing
%       duplication groups can cause issues with saveAligns.
%
% Output:
%   alignsFinal - Final alignments
%
       

%% Process inputs

% separate input into component matrices
dupGroup = alignData(:,1);
allCommon = alignData(:,2:end-1);
allConf = alignData(:,end);

% load results data
if nargin<2 || isempty(rsltFold)
    rsltFold = 'results';
end
load([rsltFold '\' rsltFold '.mat'],'results')

if nargin<3 || isempty(keepDupGroup)
    keepDupGroup = 0;
end


%% Find peak aligns

% final row indices to be removed
delRowsAll = [];
extraRows = [];
extraConf = [];

for ii = 1:max(dupGroup)
    % all indices in current duplication group
    curIdx = find(dupGroup==ii);
    
    % skip empty groups
    if isempty(curIdx); continue; end
    
    % set current rows and confidence levels
    curRows = allCommon(curIdx,:);
    curConf = allConf(curIdx,:);
    
    % initial while loop
    k = 1;
    mx = [];
    refRow = [];
    delRows = zeros(size(curConf));
    
    %% Iteratively identify indpendent alignments
    while true
        % idenify peak align
        [~,mx(k)] = max(curConf);
        refRow(k,:) = curRows(mx(k),:);
        
        % set deletion rows with matching cell values
        rmRows = curRows==refRow(k,:);
        rmRows = rmRows & curRows~=0;
        delRows(any(rmRows,2)) = 1;
        delRows(mx) = 0;
        
        % update confidence array with current align group
        curConf(mx) = 0;
        curConf(delRows==1) = 0;
        
        % exit loop when all aligns match a peak align
        if all(curConf==0)
            break
        end
        
        k = k+1;
    end

    % update global deletion rows
    delRowsAll = cat(1,delRowsAll,curIdx(delRows==1));

    % define all use indices
    useRows = curRows(delRows==0,:);

    % screen for lost indices
    lostRows = curRows;
    for jj = 1:size(lostRows,1)
        for kk = 1:size(useRows,1)
            lostRows(jj,useRows(kk,:)==lostRows(jj,:)) = 0;
        end
    end

    while true
        % isolate lost indices and confidence
        lostRows(sum(lostRows~=0,2)==0,:) = [];
        lostRows = unique(lostRows,'rows');

        % break when no lost rows
        if isempty(lostRows); break; end

        % find confidnce
        lostConf = findConfidence(lostRows,results);

        % find peak lost row
        [curConf,mx] = max(lostConf);
        curRow = lostRows(mx,:);

        % add to extra row
        extraRows(end+1,:) = curRow;
        extraConf(end+1,:) = curConf;

        % remove matched values
        for jj = 1:size(lostRows,1)
            lostRows(jj,curRow==lostRows(jj,:)) = 0;
        end
    end
    
end


%% Save output

alignsFinal = alignData;
alignsFinal(delRowsAll,:) = [];
alignsFinal(end+1:end+size(extraRows,1),2:end) = [extraRows,extraConf];

if keepDupGroup==0
    alignsFinal(:,1) = 0;
end

end

