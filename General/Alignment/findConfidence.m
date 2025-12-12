function ccConf = findConfidence(curCommon,results)

% calculate alignment confidence
nRows = size(curCommon,1);
nIdx = size(curCommon,2);
ccConf = zeros(nRows,1);

for n = 1:nRows
    grCommon = zeros(nIdx,nIdx);

    % cycle through each FOV as base
    for ii = 1:nIdx
        % get common cell value for base FOV
        CCi = curCommon(n,ii);

        % skip empty FOV
        if CCi==0; continue; end

        % get alignment index for base common cell
        cIdx = find(results.days{ii}.allRefOthersSortedInOrder(:,ii)==CCi,1);

        % cycle through each FOV as compare
        for jj = 1:nIdx
            % get common cell value for compare FOV
            CCj = curCommon(n,jj);

            % only check alignment if compare is 0 or if ii<jj so that
            % each comparison is only run once (0-0 is always skipped,
            % 0-# is only run when # is reference)
            if CCj==0 || ii<jj
                pr = sort([ii jj]);
                % get alignment for test pair (true alignment)
                testPair = results.days{ii}.allRefOthersSortedInOrder(cIdx,pr);
                % compare test pair to proposed alignment
                if isequal(testPair,curCommon(n,pr))
                    grCommon(ii,jj) = 1;
                else
                    grCommon(ii,jj) = -1;
                end
            end
        end
    end

    % calculate percent of proposed pairwise alignments that are true
    % alignments (0-0 alignments are excluded from calculation)
    ccConf(n) = sum(grCommon==1)/sum(grCommon~=0);
end

