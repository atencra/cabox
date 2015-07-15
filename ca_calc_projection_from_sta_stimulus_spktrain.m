function [xprior, xposterior] = ca_calc_projection_from_sta_stimulus_spktrain(v, stimulus, spktrain)
% ca_calc_projection_from_sta_stimulus_spktrain Projection from STA and stimulus matrix
% 
%     [xprior, xposterior] = ca_calc_projection_from_sta_stimulus_spktrain(v, stimulus, spktrain)
% 
%     v : filter, a matrix. Number of rows in v will be number of rows
%     in stimulus.
% 
%     stimulus : stimulus matrix with one frequency per row, and one time
%     bin per column. The stimulus and the filter must have the same number
%     of rows, or frequencies.
%
%     spktrain : binned spike times. Each element contains the number of 
%     spikes in the time bin. length(spktrain) = #cols in stimulus.
%  
%     xprior : distribution of all projection values of the stimulus onto the
%     filter.
% 
%     xposterior : distribution of the projection values corresponding to spike
%     times, i.e. where spktrain(i) > 0



% Check and process input arguments
[nrows, ncols] = size(v);
if ( nrows == 1 || ncols == 1 )
    error('STA should be a matrix');
end


if ( size(v,1) ~= size(stimulus,1) )
    error('Filter and stimulus must have same number of rows.');
end

[nrows, ncols] = size(spktrain);
if ( size(spktrain,1) > 1 && size(spktrain,2) > 1 )
    error('spktrain should be a vector');
end


% Calculate stimulus projection onto filter


[nfreqs, nlags] = size(v); % # frequencies, # time bins

ntrials = length(spktrain);

xraw = zeros(length(nlags:ntrials), 1);


for i = nlags:ntrials
   xraw(i) =  sum( sum ( stimulus( :, i-nlags+1:i ) .* v ) ); % inner product
end % (for i)


xprior = xraw ./ std(xraw); % scale to get units of standard deviation
xposterior = xprior( spktrain > 0 ); % includes only posterior


return;



