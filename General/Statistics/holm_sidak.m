function adjusted_pvals = holm_sidak(pvals)
%Holm-Sidak procedure for multiple Student's t tests.
%This file is applicable for equal or unequal sample sizes
% Modified 03/14/2024 - Taylor Malone, Yi Gu Lab
%

% reshape p-values as vector
k = numel(pvals);
pval_vector = reshape(pvals,k,1);

% sort p-value vector
[pval_vector,I] = sort(pval_vector);
[~,Iinv] = sort(I);

%Sidak alpha corrected values
pAdj = nan(k,1);
pAdj(1) = 1 - (1-pval_vector(1))^k;
for jj = 2:k
    pAdj(jj) = max(pAdj(jj-1),1 - (1-pval_vector(jj))^(k-jj+1));
end

adjusted_pval_vector = pAdj(Iinv);

adjusted_pvals = reshape(adjusted_pval_vector, size(pvals));

end