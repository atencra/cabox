function [xprior, xposterior] = ca_sta_stimulus_projection(sta, locator, stimulus)
% ca_sta_stimulus_projection - projection for training and test data sets
% 
%     Calculates all projections onto filters over the entire stimulus
%     duration. Each filter was calculated from 1 of the 4 training sets.
% 
%     Input arguments:
% 
%     sta : receptive field matrix.
% 
%     locator : a vector that describes whether a spike occurred during each
%     of the stimulus trials. Values are >= 0.
% 
%     stimulus : the entire stimulus that was played to the neuron. It is
%     an nf x ntrials size matrix. nf represents the number for frequencies
%     in the stimulus. ntrials is the total number of time trials.
% 
%     Output arguments:
% 
%     xprior : all projections of the stimulus onto the sta
% 
%     xposterior : projections corresponding to a spike in locator



ntrials = length(locator);


% Project all trials onto filter to get the prior distribution.
% This distribution will include the training set projections and the
% test set projections. All trials will have a projection value
% associated with them.


% Find the prior projection values for all filters at the same time:
xprior = zeros(ntrials, 1);


[nr, nc] = size(sta); % # frequencies, # time bins

for i = nc:ntrials

   xprior(i) =  sum( sum ( stimulus( :, i-nc+1:i ) .* sta ) ); % inner product

end % (for i)


xposterior = xprior( locator > 0 ); % includes only posterior


return;



