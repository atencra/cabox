function [nedata, nedata_jitter_total] = ca_calc_cell_assembly_spike_jitter(spk, stimulus, trigger, FsDVD, DF)
% ca_calc_cell_assembly_spike_jitter Estimate PCA/ICA cell assembly analysis with spike jitter control
% 
%     [cadata] = ca_calc_cell_assembly_spike_jitter(spk, stimulus, trigger, FsDVD, DF)
%     --------------------------------------------------------------
%     spk : struct array of spike times. Obtaines from saveSpikeSort.m
%     stimulus : ripple stimulus in matrix format
%     trigger : trigger for ripple stimulus
%     FsDVD : sampling of DVD Audio system
%     DF : total downsampling factor for stimulus matrix. Usually 480 or 
%       960 for 96kHz sampling rate, which is equivalent to 5ms, 10ms, etc.
%
%     cadata : struct holding the spike train matrix, the position of each
%     neuron, the cell assembly patterns, and the cell assembly activities.
%
%     This function does not calculate the STAs from the spike trains. To
%     calculate the STAs, use the function 
%        ca_calc_cell_assembly_sta_from_spk_stimulus_trigger

library('PatternJitter');
library('santosbox');

narginchk(5,5);

% Get spike train matrix.
[spkmat, position] = ca_create_spktrain_matrix_from_stimulus_spk(...
    spk, stimulus, trigger, FsDVD, DF);




% Find cell assemblies
fprintf('\nDetecting Data Cell Assemblies\n');
nedata = ca_spkmatrix_to_ensembles(spkmat);
nensembles = size(nedata.ensembles,2);
evals = nedata.eigenvalues;
lambda_min = nedata.lambda_min;
lambda_max = nedata.lambda_max;
index = find(evals<lambda_min | evals>lambda_max);
nensemble_neurons = length(index);


fprintf('\nDetecting Jittered Cell Assemblies\n');
nreps = 100;
nensembles_jitter = zeros(1,nreps);
nensemble_neurons_jitter = zeros(1,nreps);

for i = 1:nreps
    [spkmat_jitter, ~] = ca_jitter_spktrain_matrix_from_stimulus_spk(...
        spk, stimulus, trigger, FsDVD, DF);
    nedata_jitter = ca_spkmatrix_to_ensembles(spkmat_jitter);

    nensembles_jitter(i) = size(nedata_jitter.ensembles,2);

    evals = nedata_jitter.eigenvalues;
    lambda_min = nedata_jitter.lambda_min;
    lambda_max = nedata_jitter.lambda_max;
    index = find(evals<lambda_min | evals>lambda_max);
    nensemble_neurons_jitter(i) = length(index);

    nedata_jitter_total{i} = nedata_jitter;
end % (for i)






figure;
hold on;
plot(sort(nedata.eigenvalues), 'ko-', 'markersize', 2, 'markerfacecolor', 'k');
plot(sort(nedata_jitter.eigenvalues), 'ro-', 'markersize', 2, 'markerfacecolor', 'r');
plot([1 length(nedata.eigenvalues)], [nedata.lambda_max nedata.lambda_max], 'k--');
plot([1 length(nedata.eigenvalues)], [nedata.lambda_min nedata.lambda_min], 'k--');
tickpref;
xlabel('Eval #');
ylabel('Eval');


% How many ensembles?
% How many total neurons in the ensembles?
% How many neurons in each ensemble?
% Ensemble strength?


figure;
subplot(2,1,1);
hold on;
hist(nensembles_jitter, 20);
plot([nensembles nensembles], ylim, 'k-');
tickpref;
maxmax = max([nensembles_jitter(:)' nensembles(:)']);
xlim([0 1.2*maxmax]);

subplot(2,1,2);
hold on;
hist(nensemble_neurons_jitter, 20);
plot([nensemble_neurons nensemble_neurons], ylim, 'k-');
tickpref;
maxmax = max([nensemble_neurons_jitter(:)' nensemble_neurons(:)']);
xlim([0 1.2*maxmax]);




return;












