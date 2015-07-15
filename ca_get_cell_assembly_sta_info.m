function [cainfo] = ca_get_cell_assembly_sta_info(cadata, stimulus)
% ca_get_cell_assembly_sta_info  Info from cell assemblies
%
% [cainfo] = ca_get_cell_assembly_sta_info(cadata, stimulus)
% ------------------------------------------------------------------------
%

%     cadata : struct holding calculations. Has the form:
%         cadata.spktrain = binned spike train, one row per neuron
%         cadata.fsdvd = sampling rate of sound stimulus
%         cadata.df = total downsampling factor of ripple noise envelope
%         cadata.position = recording depth. position(i) -> spktrain(i,:)
%         cadata.Patterns = cell assemblies
%         cadata.Activities = time course of cell assembly activity; one row per assembly
%         cadata.nf = number of frequencies in STA
%         cadata.nlags = number of time lags in STA
%         cadata.stamat = matrix of spike train STAs.
%         cadata.ca_stamat = matrix of cell assembly STAs.
%
% stimulus : the entire ripple stimulus envelope file as one matrix. 
% Each row is one frequency, and each column is one time bin.
%



% Analysis steps:
%
% 1. Calculate projection values for both sta's and the bicellular sta
% 2. Resample from the projection values for the sta's, so that you have the
% 3. Same number of projection values as in the bicellular case
% 4. Perform resampling 10 times
% 5. Estimate information



library('cell_assembly');
library('nonlinearity')
library('santosbox')
library('strfbox');
library('UtilitiesColormaps');


spktrain = cadata.spktrain;
nlags = cadata.nlags;
nf = cadata.nf;


% stim = #trials X #dims, #dims = nfreq * nlags; resp = # neurons X #trials
[stim, resp] = ca_create_stim_trial_from_stim_matrix(stimulus, spktrain, nlags);


% Threshold Activities to get cell assembly spiking
% ca_spktrain = # neurons X #trials
[ca_spktrain] = ca_threshold_cadata_activities(cadata);
ca_resp = ca_spktrain(:, nlags:end);


fraction = [90 95 97.5 99 100]; % previously had 80 and 85 - led to bad fits 
nsamples = 10;


[nCells, nBins] = size(spktrain);
[nAssemblies, nBins] = size(ca_spktrain);


% ----------- Information from Single Units -------------
tic
for i = 1:nCells
    fprintf('Processing neuron # %.0f of %.0f\n', i, nCells);
    
    locator = resp(i,:)'; % make locator a column vector
    sta = stim' * locator; % column vector
    
    % Get all the projection values; prior = p(x), posterior = p(x|spike)
    [xprior, xposterior] = ca_get_sta_stim_projection(sta, locator, stim);
    num_cell_spikes = length(xposterior);

    % Information over different data fractions
    spk_ifraction = cell( 1, length(fraction) );
    for k = 1:length(fraction)
        spk_ifraction{k} = ca_info_fraction( xprior, xposterior, fraction(k) );
    end % (for k)

    [info_mean, info_std] = ca_info_fraction_mean_std(fraction, spk_ifraction);
    [info_extrap] = ca_info_extrapolate(fraction, info_mean);
    
    spkinfo(i).nevents = num_cell_spikes;
    spkinfo(i).ifraction = spk_ifraction;
    spkinfo(i).info_mean = info_mean;
    spkinfo(i).info_std = info_std;
    spkinfo(i).info_extrap = info_extrap;
end % (for i)
toc




% ----------- Information from Cell Assemblies -------------
tic
for i = 1:nAssemblies
    fprintf('Processing assembly # %.0f of %.0f\n', i, nCells);
    
    locator = ca_resp(i,:)'; % make locator a column vector
    sta = stim' * locator; % column vector
    
    % Get all the projection values; prior = p(x), posterior = p(x|spike)
    [xprior, xposterior] = ca_get_sta_stim_projection(sta, locator, stim);
    num_assembly_spikes = length(xposterior);

    % Information over different data fractions
    ifraction = cell( 1, length(fraction) );
    for k = 1:length(fraction)
        ifraction{k} = ca_info_fraction( xprior, xposterior, fraction(k) );
    end % (for k)

    [info_mean, info_std] = ca_info_fraction_mean_std(fraction, ifraction);
    [info_extrap] = ca_info_extrapolate(fraction, info_mean);
    
    ca_spkinfo(i).nevents = num_assembly_spikes;
    ca_spkinfo(i).ifraction = ifraction;
    ca_spkinfo(i).info_mean = info_mean;
    ca_spkinfo(i).info_std = info_std;
    ca_spkinfo(i).info_extrap = info_extrap;
    
end % (for i)
toc    


% ----------- Save data to output struct -----------
cainfo.fsdvd = cadata.fsdvd;
cainfo.df = cadata.df;
cainfo.position = cadata.position;
cainfo.nf = cadata.nf;
cainfo.nlags = cadata.nlags;
cainfo.spkinfo = spkinfo;
cainfo.ca_spkinfo = ca_spkinfo;

return;
























% 
% %--------------------------------------------------------------------
% %   Info for raw STA for bicellular, and both neurons
% %--------------------------------------------------------------------
% 
% 
% 
% 
% % Steps:
% % 1. determine number of projections in bicellular distribution X
% % 2. get sample of projections for different data fractions X
% % 3. estimate information for the different fractions X
% % 4. get extrapolate information values
% % 
% % 5. obtain sample from neuron1 and neuron2 for the same number of projections as bicellular X
% % 6. perform 2-4 for neurons 1 and 2
% % 7. repeat 5-6 10 times to get an estimate of the mean extrapolated info values
% 
% 
% 
%       % Bicellular Calculations:
%       % -----------------------------------------------------------
%       % Calculate information for different data fractions.
%       % Do the estimate for the specified number of repetitions
% 
%       bicellular_ifraction = cell( 1, length(fraction) );
%       for k = 1:length(fraction)
%          [ifrac] = ca_info_fraction( xpriorbc, xposteriorbc, fraction(k) );
%          bicellular_ifraction{k} = ifrac;
%       end % (for k)
% 
% 
% 
%       % Neuron1 Calculations:
%       % -----------------------------------------------------------
%       % Calculate information for different data fractions.
%       % Do the estimate for the specified number of repetitions
% 
%       for j = 1:nsamples
% 
%          % Get indices into a subset of the neuron1 values. The number
%          % of neuron1 values is length(xposterior1). The number of values
%          % in the subset is num_bc_spikes.
%          num_sample_spikes = min([length(xposterior1) num_bc_spikes]);
%          xposterior1_sample = randsample(xposterior1, num_sample_spikes);
% 
%          for k = 1:length(fraction)
%             [ifrac] = ca_info_fraction( xprior1, xposterior1_sample, fraction(k) );
%             neuron1_ifraction{j}{k} = ifrac;
%          end % (for j)
% 
%       end % (for j)
% 
% 
% 
% 
%       % Get mean/std of information for the data fractions
%       [idataStr] = pairs_info_fraction_mean_std(fraction, ...
%          bicellular_ifraction, neuron1_ifraction, neuron2_ifraction);
% 
% 
%       % Extrapolate the information values to get the final value for each spike train type
%       [info_extrap_bicellular, info_extrap_neuron1, info_extrap_neuron2] = ...
%          calc_pairs_info_extrapolate(fraction, idataStr);
% 
% 
% 
% 
%       %--------------------------------------------------------------------
%       %   Significant STA Info
%       %--------------------------------------------------------------------
%       [xprior1_sig, xposterior1_sig] = get_pairs_sta_bc_projection(sta1sig, locator1, stimulus);
%       [xprior2_sig, xposterior2_sig] = get_pairs_sta_bc_projection(sta2sig, locator2, stimulus);
%       [xpriorbc_sig, xposteriorbc_sig] = get_pairs_sta_bc_projection(stabcsig, locatorbc, stimulus);
%       num_bc_spikes_sig = length(xposteriorbc_sig);
% 
%       
%       
% 
%       % Bicellular Calculations:
%       % -----------------------------------------------------------
%       bicellular_ifraction_sig = cell( 1, length(fraction) );
% 
%       for k = 1:length(fraction)
%          [ifrac] = pairs_info_fraction( xpriorbc_sig, xposteriorbc_sig, fraction(k) );
%          bicellular_ifraction_sig{k} = ifrac;
%       end % (for k)
% 
% 
% 
%       % Neuron1 Calculations:
%       % -----------------------------------------------------------
%       for j = 1:nsamples
%          num_sample_spikes = min([length(xposterior1_sig) num_bc_spikes_sig]);
%          xposterior1_sample_sig = randsample(xposterior1_sig, num_sample_spikes);
% 
%          for k = 1:length(fraction)
%             [ifrac] = pairs_info_fraction( xprior1_sig, xposterior1_sample_sig, fraction(k) );
%             neuron1_ifraction_sig{j}{k} = ifrac;
%          end % (for j)
%       end % (for j)
% 
%       % Get mean/std of information for the data fractions
%       [idataStr_sig] = pairs_info_fraction_mean_std(fraction, ...
%          bicellular_ifraction_sig, neuron1_ifraction_sig, neuron2_ifraction_sig);
% 
%       % Extrapolate the information values to get the final value for each spike train type
%       [info_extrap_bicellular_sig, info_extrap_neuron1_sig, info_extrap_neuron2_sig] = ...
%          calc_pairs_info_extrapolate(fraction, idataStr_sig);
% 
% 
% 
% 
%       % Save the data
%       %---------------------------------------------------------------
% 
%       % Data from raw STA
%       cainfo(i).fraction = fraction;
% 
%       cainfo(i).num_bc_spikes = num_bc_spikes;
% 
%       cainfo(i).bicellular_info_frac_mn = idataStr.bicellular_info_frac_mn;
%       cainfo(i).bicellular_info_frac_std = idataStr.bicellular_info_frac_std;
% 
%       cainfo(i).neuron1_info_frac_mn = idataStr.neuron1_info_frac_mn;
%       cainfo(i).neuron1_info_frac_std = idataStr.neuron1_info_frac_std;
% 
%       cainfo(i).bicellular_info_extrap = info_extrap_bicellular;
%       cainfo(i).neuron1_info_extrap = info_extrap_neuron1;
% 
% 
%       % Data from significant STA
%       cainfo(i).num_bc_spikes_sig = num_bc_spikes_sig;
% 
%       cainfo(i).bicellular_info_frac_mn_sig = idataStr_sig.bicellular_info_frac_mn;
%       cainfo(i).bicellular_info_frac_std_sig = idataStr_sig.bicellular_info_frac_std;
% 
%       cainfo(i).neuron1_info_frac_mn_sig = idataStr_sig.neuron1_info_frac_mn;
%       cainfo(i).neuron1_info_frac_std_sig = idataStr_sig.neuron1_info_frac_std;
% 
%       cainfo(i).bicellular_info_extrap_sig = info_extrap_bicellular_sig;
%       cainfo(i).neuron1_info_extrap_sig = info_extrap_neuron1_sig;
% 
% 
% 
% function [idataStr] = ca_info_fraction_mean_std(fraction, ...
%    bicellular_ifraction, neuron1_ifraction, neuron2_ifraction)
% %train_test_info_fraction_mean_std - mean, standard dev of information 
% % data fraction calculations.
% %
% % Information values were calculated for different fraction of the data for
% % each of 4 training or test data sets. For each fraction, 5 randomized 
% % estimates of information were calculated. This function takes the mean 
% % of the 5 estimates. It also takes another mean, this time across the 4 
% % data sets.
% %
% % [info_frac_mn, info_frac_std, info_mtx] = train_test_info_fraction_mean_std(fraction, info_fraction)
% % ----------------------------------------------------------------------------------------------------
% % fraction : proportions of the data used. A vector like [80 90 95 97.5
% % 100], or something similar.
% %
% % info_fraction : 4xlength(fraction) cell array. Information calculations 
% % for the 4 different data sets, with 5 information estimates for each
% % data fraction. Thus, each element of info_fraction is a 1x5 vector.
% %
% % info_frac_mn : mean information values across the 4 data sets for each
% % data fraction. A vector of length(fraction).
% %
% % info_frac_std : same as info_frac_mn, except for standard deviation.
% %
% % info_mtx : 4xlength(fraction) matrix of information values. Each element
% % represents the mean of 5 information values for a given data fraction.
% % Each row represent one of the 4 data sets.
% %
% % caa 3/13/09
% 
% fprintf('\nRunning train_test_info_fraction_mean_std ...\n');
% 
% 
% bicellular_info_frac_mn = zeros(size(fraction));
% bicellular_info_frac_std = zeros(size(fraction));
% 
% for i = 1:length(fraction)
% 
%    temp = bicellular_ifraction{i};
% 
%    bicellular_info_frac_mn(i) = mean(temp);
%    bicellular_info_frac_std(i) = std(temp);
% 
% end % (for i)
% 
% 
% neuron1_info_mtx = zeros( length(neuron1_ifraction), length(fraction) );
% 
% for i = 1:length(neuron1_ifraction) % for each data set
% 
%    for j = 1:length(fraction) % for each fraction of the data
% 
%       temp = neuron1_ifraction{i}{j};
%       neuron1_info_mtx(i, j) = mean(temp);
% 
%    end % (for j)
% 
% end % (for i)
% 
% 
% neuron1_info_frac_mn = mean( neuron1_info_mtx, 1 );
% neuron1_info_frac_std = std( neuron1_info_mtx, 0, 1 );
% 
% 
% neuron2_info_mtx = zeros( length(neuron2_ifraction), length(fraction) );
% 
% for i = 1:length(neuron2_ifraction) % for each data set
% 
%    for j = 1:length(fraction) % for each fraction of the data
% 
%       temp = neuron2_ifraction{i}{j};
%       neuron2_info_mtx(i, j) = mean(temp);
% 
%    end % (for j)
% 
% end % (for i)
% 
% 
% neuron2_info_frac_mn = mean( neuron2_info_mtx, 1 );
% neuron2_info_frac_std = std( neuron2_info_mtx, 0, 1 );
% 
% 
% 
% idataStr.bicellular_info_frac_mn = bicellular_info_frac_mn;
% idataStr.bicellular_info_frac_std = bicellular_info_frac_std;
% 
% idataStr.neuron1_info_frac_mn = neuron1_info_frac_mn;
% idataStr.neuron1_info_frac_std = neuron1_info_frac_std;
% 
% idataStr.neuron2_info_frac_mn = neuron2_info_frac_mn;
% idataStr.neuron2_info_frac_std = neuron2_info_frac_std;
% 
% 
% 
% return;
% 
% 













