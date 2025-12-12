function C = combCells(A,B)
%% combCells
% Combines two non-overlaping cell arrays of the same size into a single
% cell array
%

% combine disjoint inputs
C = cell(size(A));
for ii = 1:numel(A)
    if isempty(A{ii})
        C{ii} = B{ii};
    elseif isempty(B{ii})
        C{ii} = A{ii};
    else
        error('Cell array entries must not overlap')
    end
end

end

