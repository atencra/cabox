function [subsetinfo] = ca_calc_sta_info_from_subsets(exp_site_cadata)
% ca_get_cell_assembly_sta_info  Info from cell assemblies
%
% [subsetinfo] = ca_calc_sta_info_from_subsets(exp_site_cadata)
% ------------------------------------------------------------------------
%
%



library('cell_assembly');
library('infobox');


fraction = [90 95 97.5 99 100];
nsamples = 10;
spike_count_threshold = 200; % only calc info if each subset has this many spikes

cadata = exp_site_cadata.cadata;
spktrain = cadata.spktrain;
nlags = cadata.nlags;
nf = cadata.nf;
nlags = cadata.nlags;
dft = exp_site_cadata.df;
dff = 5;
rn = exp_site_cadata.stim; % stim as a string
rn = str2num(rn(3:end)); % get RN number

[stimstr] = ca_get_ripple_noise_stimulus([], rn, dft, dff );

spksubset = exp_site_cadata.spksubset;

for i = 1:length(spksubset)

    subsetinfo(i).exp = exp_site_cadata.exp;
    subsetinfo(i).site = exp_site_cadata.site;
    subsetinfo(i).stim = exp_site_cadata.stim;
    subsetinfo(i).df = exp_site_cadata.df;
    subsetinfo(i).fs = exp_site_cadata.fs;
    subsetinfo(i).probetype = exp_site_cadata.probetype;
    subsetinfo(i).depth = exp_site_cadata.depth;
    subsetinfo(i).neuron = spksubset(i).neuron;
    subsetinfo(i).assembly = spksubset(i).assembly;

    % Get binned spike trains for each subset
    locator_all = spksubset(i).spktrain_all;
    locator_wo_CA = spksubset(i).spktrain_wo_CA;
    locator_w_CA = spksubset(i).spktrain_w_CA;


    fprintf('CA SUBSET INFO: %.0f of %.0f\n', i, length(spksubset));

    if ( sum(locator_all) > spike_count_threshold && ...
         sum(locator_wo_CA) > spike_count_threshold && ...
         sum(locator_w_CA) > spike_count_threshold )

        % Estimate STAs for each subset
        sta_all = ca_sta_from_locator_stimulus(locator_all, stimstr.stimulus, nlags);
        sta_wo_CA = ca_sta_from_locator_stimulus(locator_wo_CA, stimstr.stimulus, nlags);
        sta_w_CA = ca_sta_from_locator_stimulus(locator_w_CA, stimstr.stimulus, nlags);

         % Plot STAs to make sure everything is working okay.
         % Uncomment below to check a neuron.
%         figure;
%         subplot(1,3,1); imagesc(sta_all);
%         subplot(1,3,2); imagesc(sta_wo_CA);
%         subplot(1,3,3); imagesc(sta_w_CA);
%         pause;


        % Get all the projection values; prior means without regard to a
        % spike; posterior means for a spike
        [xprior_all, xposterior_all] = ca_sta_stimulus_projection(sta_all, locator_all, stimstr.stimulus);
        [xprior_wo_CA, xposterior_wo_CA] = ca_sta_stimulus_projection(sta_wo_CA, locator_wo_CA, stimstr.stimulus);
        [xprior_w_CA, xposterior_w_CA] = ca_sta_stimulus_projection(sta_w_CA, locator_w_CA, stimstr.stimulus);



        % Determine minimum number of samples/spikes. Datasets will be 
        % sampled based on this number.
        min_spikes = min([length(xposterior_all) length(xposterior_wo_CA) length(xposterior_w_CA)]);


      % Calculate info for different data fractions
      % -----------------------------------------------------------
        all_ifraction = ca_subset_info_from_data_fraction(xprior_all, xposterior_all, fraction, nsamples, min_spikes);
        wo_CA_ifraction = ca_subset_info_from_data_fraction(xprior_wo_CA, xposterior_wo_CA, fraction, nsamples, min_spikes);
        w_CA_ifraction = ca_subset_info_from_data_fraction(xprior_w_CA, xposterior_w_CA, fraction, nsamples, min_spikes);



        % Get mean/std of information for the data fractions

        [all_info_frac_mn, all_info_frac_std] = ...
            info_from_data_fraction_mean_std(fraction, all_ifraction);

        [wo_CA_info_frac_mn, wo_CA_info_frac_std] = ...
            info_from_data_fraction_mean_std(fraction, wo_CA_ifraction);

        [w_CA_info_frac_mn, w_CA_info_frac_std] = ...
            info_from_data_fraction_mean_std(fraction, w_CA_ifraction);


        % Extrapolate the information values to get the final value for each spike train type
        [all_info_extrap] = info_extrapolate_from_mn_std(fraction, all_info_frac_mn);
        [wo_CA_info_extrap] = info_extrapolate_from_mn_std(fraction, wo_CA_info_frac_mn);
        [w_CA_info_extrap] = info_extrapolate_from_mn_std(fraction, w_CA_info_frac_mn);


        fprintf('\n');
        fprintf('Information for subsets:\n');
        fprintf('All: %.3f bits/spk\n', all_info_extrap);
        fprintf('WO : %.3f bits/spk\n', wo_CA_info_extrap);
        fprintf('W  : %.3f bits/spk\n\n', w_CA_info_extrap);


        % Save the data
        %---------------------------------------------------------------

        % Data from raw STA
        subsetinfo(i).fraction = fraction;
        subsetinfo(i).assembly = spksubset(i).assembly;
        subsetinfo(i).min_spikes = min_spikes;

        subsetinfo(i).all_info_frac_mn = all_info_frac_mn;
        subsetinfo(i).all_info_frac_std = all_info_frac_std;
        subsetinfo(i).all_info_extrap = all_info_extrap;

        subsetinfo(i).w_CA_info_frac_mn = w_CA_info_frac_mn;
        subsetinfo(i).w_CA_info_frac_std = w_CA_info_frac_std;
        subsetinfo(i).w_CA_info_extrap = w_CA_info_extrap;


        subsetinfo(i).wo_CA_info_frac_mn = wo_CA_info_frac_mn;
        subsetinfo(i).wo_CA_info_frac_std = wo_CA_info_frac_std;
        subsetinfo(i).wo_CA_info_extrap = wo_CA_info_extrap;


    else % if there were no common spikes, then save as empty matrices

        % Save the data
        %---------------------------------------------------------------


        % Data from raw STA
        subsetinfo(i).fraction = [];
        subsetinfo(i).assembly = [];
        subsetinfo(i).min_spikes = [];

        subsetinfo(i).all_info_frac_mn = [];
        subsetinfo(i).all_info_frac_std = [];
        subsetinfo(i).all_info_extrap = [];

        subsetinfo(i).w_CA_info_frac_mn = [];
        subsetinfo(i).w_CA_info_frac_std = [];
        subsetinfo(i).w_CA_info_extrap = [];


        subsetinfo(i).wo_CA_info_frac_mn = [];
        subsetinfo(i).wo_CA_info_frac_std = [];
        subsetinfo(i).wo_CA_info_extrap = [];


   end % if

%    if ( ~mod(i,25) )
%       save('temp-bicellular-sta-info-data-2', 'bcinfo', 'i');
%    elseif ( i == length(bcstr) )
%       save('temp-bicellular-sta-info-data-2', 'bcinfo', 'i');
%    end

end % (for i)

return;

















