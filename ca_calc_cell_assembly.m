function [cadata] = ca_calc_cell_assembly(spk, stimulus, trigger, FsDVD, DF)
% cell_assembly_sta Estimate PCA/ICA cell assembly analysis
% 
%     [cadata] = cell_assembly_sta(spk, stimulus, trigger, FsDVD, DF)
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
%     calculate the STAs, use the function ca_calc_cell_assembly_sta.m

% Make sure there are 5 input arguments.
narginchk(5,5);

% Get spike train matrix.
[spktrain_matrix, position] = ...
    ca_create_spktrain_matrix_from_stimulus_spk(spk, stimulus, trigger, FsDVD, DF);


% Find cell assemblies
fprintf('\nDetecting Cell Assemblies\n');
[Patterns, Activities] = ca_detect_cell_assemblies_data(spktrain_matrix);

ca_plot_cell_assemblies_data(spktrain_matrix, Patterns, Activities, position);



% Assign data to struct for output argument
cadata.spktrain = spktrain_matrix;
cadata.position = position;
cadata.Patterns = Patterns;
cadata.ActivitiesPatterns = Activities;


return;





















