function [W, chi2, pval] = KendallW(X)
% Compute the Kendall's coefficient of concordance of the matrix X.
% E.g. W = KendallW(RankMatrix)
%
% Input:
%           X must be a N-by-K matrix, N is the number of
%           "candidate" and the K is the number of "judge"
% Outputs:
%           W = Kendall's coefficient of concordance
%
% Created by Lijie Huang, 2010/6/5
% Modified by Craig Atencio, 2/1/16
%==========================================================================


% Method 3
[nObs,nVars] = size(X);
RankMatrix = tiedrank(X);

bigT = 0;
for j = 1:size(RankMatrix,2)
    [freq,prop] = hist(RankMatrix(:,j));
    freq = freq(freq>1);
    bigT = bigT + sum(freq.^3 - freq);
end

ranksum = sum(RankMatrix,2);
S1 = sum(ranksum.^2,1);
S2 = (sum(ranksum))^2;
S = S1 - S2 / nObs;
temp = nObs^3 - nObs;
W = 12*S/(nVars^2*(nObs^3-nObs) - nVars*bigT);

chi2 = nVars*(nObs-1)*W;

pval = 1 - chi2cdf(chi2,nObs-1);

return;











