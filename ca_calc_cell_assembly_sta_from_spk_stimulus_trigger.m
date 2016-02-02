function [cadata] = ca_calc_cell_assembly_sta_from_spk_stimulus_trigger(spk, stimulus, trigger, FsDVD, totalDF)
% ca_calc_cell_assembly_sta_from_spk_stimulus_trigger Estimate STAs from PCA/ICA cell assembly analysis
% 
%     [cadata] = ca_calc_cell_assembly_sta_from_spk_stimulus_trigger(spk, stimulus, trigger, FsDVD, totalDF)
%     ---------------------------------------------------------------------------------------------
%     spk : struct array of spike times. Obtaines from saveSpikeSort.m
%
%     stimulus : ripple stimulus in matrix format. Each row is one
%     frequency, and each column is one time bin.
%
%     trigger : trigger for ripple stimulus
%
%     FsDVD : sampling of DVD Audio system
%
%     totalDF : The total downsampling factor for the stimulus matrix. 
%     A scalar. This will be the downsampling factor for the original .spr 
%     file, as well as any additional downsampling. The original 
%     downsampling will usually be 48, and additional downsampling will be 
%     [5, 10, 20, 40, or 100]. 
%
%     Therefore, the totalDF would be 48 * [5, 10, 20, 40, or 100], which
%     is equivalent to 2.5 ms, 5 ms, etc. for 96 kHz sound stimulus sampling
%     rate.
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


ca;

% Make sure there are 5 input arguments.
narginchk(5,5);


% Get spike train matrix.
[spktrain_matrix, position] = ...
    ca_create_spktrain_matrix_from_stimulus_spk(spk, stimulus, trigger, FsDVD, totalDF);


% Find cell assemblies
fprintf('\nDetecting Cell Assemblies\n');
[Patterns, Activities] = ca_detect_cell_assemblies_data(spktrain_matrix);
ca_plot_cell_assemblies_data(spktrain_matrix, Patterns, Activities, position);


% Get stimulus-response matrices for later STA calculations
fprintf('\nForming Stimulus observation matrix\n');
if ( (totalDF / FsDVD * 1000) < 5 )
    nlags = ceil( 100 / (totalDF / FsDVD * 1000) ); % 100 ms / time bin size
else
    nlags = 20;
end

nf = size(stimulus, 1);



% Calculate STA for each neuron. Each row of spktrain_matrix holds the
% spike train for one neuron.
stamat = zeros(size(spktrain_matrix,1),size(stimulus,1)*nlags);
for i = 1:size(spktrain_matrix, 1)
    sta = ca_calc_sta_from_stimulus_spktrain(stimulus, spktrain_matrix(i,:), nlags);
    stamat(i,:) = sta(:)';
end


% Calculate STA for each cell assembly. Each row of Activities holds the
% activations for one cell assembly.
ca_stamat = zeros(size(Activities,1),size(stimulus,1)*nlags);
for i = 1:size(Activities, 1)
    sta = ca_calc_sta_from_stimulus_spktrain(stimulus, Activities(i,:), nlags);
    ca_stamat(i,:) = sta(:)';
end



% Assign data to struct for output argument
cadata.spktrain = spktrain_matrix;
cadata.fsdvd = FsDVD;
cadata.df = totalDF;
cadata.position = position;
cadata.Patterns = Patterns;
cadata.Activities = Activities;
cadata.nf = nf;
cadata.nlags = nlags;
cadata.stamat = stamat;
cadata.ca_stamat = ca_stamat;

return;




% Code for significance testing, if we want it.
% nreps = 4;
% pval = 0.05;
% for i = 1:size(stamat,1)
%     [sta_sig, siglevel, rand_dist] = ...
%         ca_sig_sta_from_stim_obs_resp(stamat(i,:), resp(i,:), stim, nreps, pval);
%     figure;
%     subplot(1,2,1);
%     imagesc(reshape(stamat(i,:),nf,nlags));
%     subplot(1,2,2);
%     imagesc(reshape(sta_sig,nf,nlags));
% end % (for i)




















