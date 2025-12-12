function tileInv = invPercentile(probeData,refData)
%% invPercentile
% calculates the inverse percentile, i.e. at what percetile of the
% reference data, the probeData falls.
%
% Inputs:
%   probeData: data must be assessed within reference population. Must by a
%       1-by-n array.
%   refData: data to use as percentile reference. Must be an m-by-n matrix.
%

refSort = sort(refData,1);

curLess = sum(refSort<probeData,1);
curEqual = sum(refSort==probeData,1);

tileInv = (curLess+curEqual/2)/size(refData,1)*100;
if any(tileInv>100)
    error('Invalid Percentile')
end

end

