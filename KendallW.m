function [w, chi2, pval] = kendallw(x)
% kendallw  Compute the Kendall's coefficient of concordance, W, of the matrix X.
% 
%     [w, chi2, pval] = kendallw(x)
% 
%     x : m x n matrix. m = number of observations, n = number of
%     variables.
% 
%     w : Kendall Coefficient of Concordance statistic. w is adjusted for
%     ties in ranks for each column of x.
% 
%     chi2 : chi-squared statistic for w. chi2 = n * (m-1) * w
% 
%     pval : significance value associated with chi2. The test is
%     one-tailed, since it only measures positive associativity.
% 
%     Craig Atencio, 2/1/16


% Method 3
[nObs,nVars] = size(x);
rankmat = tiedrank(x);

% sum ranks across columns, then take overall mean
rowSum = sum(rankmat,2);
rowSumMean = mean(rowSum);

% deviance of row sums from overall mean
S = sum( (rowSum - rowSumMean).^2 );

% calculate adjustment factor for ties
bigT = 0;
for j = 1:size(rankmat,2)
    [freq,prop] = hist(rankmat(:,j));
    freq = freq(freq>1);
    bigT = bigT + sum(freq.^3 - freq);
end

% adjusted w statistic, chi-squared statistic, and pval for H0: w = 0
w = 12*S / (nVars^2*(nObs^3-nObs) - nVars*bigT);
chi2 = nVars*(nObs-1)*w;
pval = 1 - chi2cdf(chi2,nObs-1);


return;











