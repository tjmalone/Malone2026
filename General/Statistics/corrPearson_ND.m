function coef = corrPearson_ND(x,y,dim)
%% corrPearson_ND
% Implements a simplified version of Matlab's internal corr function.
% Calculates the Pearson correlation for an N-D array along a specific
% dimension. Currently cannot use alternate correlation methods, caluclate
% p-values or handle NaNs.

% set default dimension
if nargin<3 || isempty(dim)
    dim = 1;
end

% remove means
xz = bsxfun(@minus, x, mean(x,dim));
yz = bsxfun(@minus, y, mean(y,dim));

% standard Pearson correlation coefficient formula
x2 = xz .^ 2;
y2 = yz .^ 2;
ab = xz .* yz;
coef = squeeze(sum(ab, dim) ./ sqrt(sum(x2, dim) .* sum(y2, dim)));

end

