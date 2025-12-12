function alignsFinal = saveAligns(alignData,svFolder)
%% saveAligns.m
% Takes combined alignmnet data matrix and saves selected aligns in new
% matrix. Selected cells are all cells with non-zero value in column 1.
% This can be manually pasted into matlab matrix or calculated using
% peakAligns.
%
% Inputs:
%   alignData - combined data matrix contains the following:
%       dupGroup - defines cells in common duplication groups
%       allCommon - defines all possible alignments
%       allConf - confidence levels for the corresponding alignment
%   svFolder - folder to save final alignments to. Should be the same
%       folder that alignment resilts are saved in. Default is current
%       directory.
%
% Output:
%   alignsFinal - Final alignments. Does not include duplication group or
%       alignment confidence. Is also saved in svFolder
%
       

% prune unselected  cells and extraneous columns
useIdx = alignData(:,1)==0;
alignsFinal = alignData(useIdx,2:end-1);
alignsFinal = sortrows(alignsFinal,1);

% set output folder
if nargin<2 || isempty(svFolder)
    svFolder = '';
end

% save output in desired folder
save([pwd '\' svFolder '\alignsFinal.mat'],'alignsFinal')

end

