function idxLog = idx2common(curCellSel,alignsLearning)
%% idx2common
% Converts true indices of common cells to a logical index  of only the set
% of commonn cells.

nFOV = size(curCellSel,1);
nDays = size(curCellSel,2);

idxLog = cell(nFOV,nDays);

for ff = 1:nFOV
    curSel = curCellSel{ff,1};
    curAlign = alignsLearning{ff}(:,1);

    idxCur = ismember(curAlign,curSel);

    for dd = 1:nDays
idxLog{ff,dd} = idxCur;

    end
end

end

