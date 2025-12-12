function pShow = generatePShow(nGroups,cmps)
% Automatically generate bar combinations (pShow) for use with barGroup.All
% inputs are required.
%
% Inputs:
%   nGroups - number of bar groups 
%   cmps - comparisons desired within each group
%

nCmp = size(cmps,1);
pShow = zeros(nCmp*nGroups,2);

idx = 1;
for ii = 1:nGroups
    curP = (cmps-1)*nGroups+ii;
    pShow(idx:idx+nCmp-1,:) = curP;
    idx = idx+nCmp;
end

end