function [xprior, xposterior] = ca_get_sta_stim_projection(v, spktrain, stim)
% ca_get_sta_stim_projection Projection from STA and stimulus matrix
% 
%     [xprior, xposterior] = ca_get_sta_stim_projection(v, spktrain, stim)
% 
%     v : filter, a vector. If sta was originally a matrix
%     it will usually be obtained from sta = sta(:)
% 
%     spktrain : binned spike times. Each element contains the number of 
%     spikes in the time bin.
% 
%     stim : stimulus matrix with one trial per row. Will have the size:
%     Ntrials X Ndims, where Ndims = #freqs * #time bins in the STA.
% 
%     xprior : distribution of all projection values of the stimulus onto the
%     STA.
% 
%     xposterior : distribution of the projection values corresponding to spike
%     times, i.e. where spktrain(i) > 0



% Check and process input arguments
[nrows, ncols] = size(v);
if ( nrows > 1 && ncols > 1 )
    error('STA should be a vector');
end

if ( nrows < ncols )
    v = v(:);
end

[nrows, ncols] = size(stim);
if ( ncols > nrows )
    error('Reshape stim to have one observation per row.');
end


[nrows, ncols] = size(spktrain);
if ( nrows > 1 && ncols > 1 )
    error('spktrain should be a vector');
end

if ( nrows < ncols )
    spktrain = spktrain(:);
end


xraw = stim * v; % all projection values
xprior = xraw ./ std(xraw); % scale to get units of standard deviation
xposterior = xprior( spktrain > 0 ); % includes only posterior





% fprintf('\nRunning train_test_projection ...\n\n');
% 
% 
% ntrials = length(locator);
% 
% % Prior distribution of projection values for test set spikes the 
% % onto training set STA
% 
% 
% % Project all trials onto filter to get the prior distribution.
% % This distribution will include the training set projections and the
% % test set projections. All trials will have a projection value
% % associated with them.
% 
% 
% % Find the prior projection values for all filters at the same time:
% xprior = zeros(ntrials, 1);
% 
% 
% [nr, nc] = size(sta); % # frequencies, # time bins
% 
% for i = nc:ntrials
% 
%    xprior(i) =  sum( sum ( stimulus( :, i-nc+1:i ) .* sta ) ); % inner product
% 
%    if ( ~mod(i, 100000) )
%       fprintf('i = %.0f\n', i);
%    end
% end % (for i)
% 
% 
% xposterior = xprior( locator > 0 ); % includes only posterior
% 
% % [length(xtrain) length(xtest) length(xtrain_locator)
% % length(xtest_locator)]


return;



