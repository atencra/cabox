function ifraction = ca_subset_info_from_data_fraction(xprior, xposterior, fraction, nsamples, min_spikes)
% ifraction = ca_subset_info_from_data_fraction(xprior, xposterior, fraction, nsamples, min_spikes)
% 
%     xprior : prior distribution of projection values
% 
%     xposterior : distribution for spikes
% 
%     fraction : data percentages to estimate info. Example: [90 92.5 95 97.5 100]
% 
%     min_spikes : minimum number of spikes for subsets of data
% 
%     nsamples : number of times to resample data when the number of spikes is
%     not equal to the minimum number.
% 
%     ifraction : jxk cell array. j indicates sample number, k indicates fraction
%     number.


if ( length(xposterior) == min_spikes )
    for j = 1:1
        xposterior_sample = xposterior;
        for k = 1:length(fraction)
            [ifrac] = info_from_data_fraction( xprior, xposterior_sample, fraction(k) );
            ifraction{j}{k} = ifrac;
        end % (for j)
    end % (for j)

else 
    for j = 1:nsamples
        num_sample_spikes = min([length(xposterior) min_spikes]);
        xposterior_sample = randsample(xposterior, num_sample_spikes);
        for k = 1:length(fraction)
            [ifrac] = info_from_data_fraction( xprior, xposterior_sample, fraction(k) );
            ifraction{j}{k} = ifrac;
        end % (for j)
    end % (for j)
end

return;





