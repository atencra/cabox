function [sta_sig, siglevel, rand_dist] = ca_sig_sta_from_stim_obs_resp(sta, locator, stim, nreps, pval)
% rand_test_sta_locator Significance test for STA calculated using a locator
% 
%    [sta_sig, siglevel, rand_dist] = ...
%       rand_test_sta_locator(sta, locator, stimulus, nreps, pval)
%
%    computes a randomization test on STA to determine the significant pixels
%    in STA. The test is performed by shifting the spike train, in LOCATOR, by
%    an amount equal to 
% 
%       shiftsize = round( length(locator)/(nreps+1) );
% 
%    and then calculating a random STA. NREPS determines the shift size of the
%    spike train, and also the size of random distribution. PVAL is the
%    significance level of the test.
% 
%    STA_SIG is the significant part of STA at the level PVAL. SIGLEVEL is 
%    the level that was used to determine signifcant parts of STA, and
%    RAND_DIST is the distribution of randomized STA values. It is a vector 
%    of length = NREPS*length(STA(:))
% 
%    Inputs
%    ----------------------------------------------------------
%    sta : spike-triggered average previously estimated
%    locator : spike train of 1's and 0's
%    stim : matrix of the stimulus. Each row is one trial.
%    nreps : number randomization to perform
%    pval : significance level
% 
%    Outputs
%    ----------------------------------------------------------
%    sta_sig : significant parts of sta
%    siglevel : value used to threshold sta
%    rand_dist : distribution used to obtain siglevel
% 
%    Craig Atencio
%    2/1/13

size(sta)

shiftsize = round( length(locator)/(nreps+1) );

rand_dist = [];

for i = 1:nreps
    shift = i * shiftsize;
    loc_rand = circshift(locator, [0 shift]);
    sta_rand = loc_rand * stim; % one sta per row
    rand_dist = [rand_dist; abs( sta_rand(:) ) ];
end

rand_dist = sort(rand_dist);


sig_index = ceil(pval*length(rand_dist));
siglevel = rand_dist(end-sig_index);

sta_abs = abs(sta);
sta_sig = zeros(size(sta));
sta_sig(sta_abs > siglevel) = sta(sta_abs > siglevel);

return;






