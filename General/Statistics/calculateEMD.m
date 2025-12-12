function emd = calculateEMD(P, Q)
% calculate the 1D earth movers distance between two vectors. When P and Q
% are matrices, will calculate the emd for every pair of columns

% get number of columns
colsP = size(P,2);
colsQ = size(Q,2);

% normalize distributions
useColsP = ~all(P==0,1);
useColsQ = ~all(Q==0,1);
P(:,useColsP) = P(:,useColsP)./sum(P(:,useColsP),1);
Q(:,useColsQ) = Q(:,useColsQ)./sum(Q(:,useColsQ),1);

% Compute the cumulative distribution functions (CDFs)
cdfP = cumsum(P);
cdfQ = cumsum(Q);

% Calculate the EMD as the sum of absolute differences between the CDFs
emd = zeros(colsP,colsQ);
for pp = 1:colsP
    for qq = 1:colsQ
        emd(pp,qq) = sum(abs(cdfP(:,pp) - cdfQ(:,qq)));
    end
end

end

