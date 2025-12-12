function R2 = findR2(A,B)
%% Rsquare
% Calculates R-square value between two vectors or matrices. If inputs are
% matrics, each column should be one sample. Any rows with NaN in either
% input will be fully ignored.

useIdx = any(~isnan(A),2) & any(~isnan(B),2);
useA = A(useIdx,:);
useB = B(useIdx,:);

Re = sum((useA-useB).^2,1);

Rt = sum((useA-mean(useA)).^2,1);

R2 = 1-Re./Rt;

end