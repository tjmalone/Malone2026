function p = prctileInv(A,P)
%% prctileInv
% Inverse of built-in MATLAB function prctile. Return the percentages of
% percentile elements P within the input data A.
%
% Inputs:
%       A - Input data. Must be 1D numerical array.
%       P - Pecentile values at which to calculate percentages. Can be
%       scalar or 1D array.
%
% Outputs:
%       p - percentages of percentile elements within A.
%

% sort input data
sortedA = sort(A);

% calculate unique values
[C,~,idx] = unique(sortedA,'stable');

% calculate input data percentages
listPercentages = (100/length(A))*(0.5:1:length(A)-0.5);

% calculate mean percentage of unique data values
val = accumarray(idx,listPercentages,[],@mean);

if length(val)>1
    % calculate percentages at input percentiles
    p = interp1(C,val,P);

    % find percentiles outside range
    Pmin = P<C(1);
    Pmax = P>C(end);

    % set percentages for percentiles outside range to 0 or 100
    p(Pmin) = 0;
    p(Pmax) = 100;
else
    % set p if every position is equal
    p = val*ones(size(P));
end

end
