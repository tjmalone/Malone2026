function alignData = uniqueAligns(useIdx,rsltFold,fMaster)
%% uniqueAligns
% Generates a preliminary alignment array for a  given set of alignment
% results. 
%
% Input: Inputs can be left empty or skipped to use default behavior.
%   useIdx - Alignment set numbers to include in analysis. Cannot exceed
%       the max number of FOV in results set. Example: 1 = all cells
%       aligned to only 1 FOV. 11 = all cells aligned to 11 FOV. Default is
%       all alignment set numbers
%   svID - name alignment set to analyze. Should be the full folder/file
%       name generated in alignmentPost. Default is 'results'
%   fMaster - folder to perform analysis in. Default is current folder
%
% Outputs:
%   alignData - alignment info array.
%       dupGroup: First column is duplication group. All rows with the same
%           duplication group are duplicates or are linked to mutual
%           duplictaes.
%       allCommon: Middle columns are cell IDs for all possible alignments
%       allConf - Last column is confidence levels for the corresponding
%           alignment
%


%% Process inputs

% set primary folder
if nargin<3 || isempty(fMaster)
    fMaster = pwd;
end
orFolder = pwd;
cd(fMaster)

% load results data
if nargin<2 || isempty(rsltFold)
    rsltFold = 'results';
end
load([rsltFold '\' rsltFold '.mat'],'results')

% define alignment set numbers to use
if nargin<1 || isempty(useIdx)
    useIdx = 1:length(results.days);
end


%% Generate combined data matrix

% initialize combined common cell/confidence matrices
allCommon = [];
allConf = [];

% concatenate common cell/confidence matrices for all use indices
for ii = 1:length(useIdx)
    curCommon = results.(['commonCells_' num2str(useIdx(ii))]);
    curConf = results.(['CCconfidence_' num2str(useIdx(ii))]);
    
    allCommon = cat(1,allCommon,curCommon);
    allConf = cat(1,allConf,curConf);
end

% FOV number number
nFOV = size(allCommon,2);


%% Find duplicates

% initialize duplication group array and ID
dupGroup = zeros(size(allConf));
cG = 1;

% cycle through each FOV as reference
for ii = 1:nFOV
   % identify all unique cell values in reference FOV
   [unqCur] = unique(allCommon(:,ii),'first');
    
   % cycle through unique cell values
   for jj = 1:length(unqCur)
       % skip all 0 alignments
       if unqCur(jj)==0; continue; end
       
       % find all duplicates of current cell value
       dups = find(allCommon(:,ii)==unqCur(jj));
       
       % skip truly unique cell values
       if length(dups)==1; continue; end
       
       % identify previous duplicate groups of current duplicates
       curDups = dupGroup(dups);
       
       % check if any duplicate was part of previous duplication group
       if any(curDups~=0)
           % if yes, get minimum group number (excluding 0)
           curDupsI = curDups;
           curDupsI(curDupsI == 0) = inf;
           minDup = min(curDupsI);
           
           % reset the duplication group number for all linked groups to
           % new minimum
           for kk = 1:length(curDups)
               if curDups(kk)==0; continue; end
               
               idxDup = dupGroup==curDups(kk);
               dupGroup(idxDup) = minDup;
           end           
               
       else
           % if no, set duplication group number to new value
           minDup = cG;
           cG = cG +1;
       end
       
       % update duplication group number for current duplications
       dupGroup(dups) = minDup;
   end

end

% concatenate duplication group, alignments, and confidence
alignData = cat(2,dupGroup,allCommon,allConf);

% return to original folder
cd(orFolder)

end

