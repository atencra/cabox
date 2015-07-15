function fio = ca_calc_cell_assembly_nonlinearity_from_cadata_stimulus(cadata, stimulus)
% ca_calc_cell_assembly_nonlinearity_from_cadata_stimulus Estimate STA nonlinearities for PCA/ICA cell assemblies
% 
%     fio = ca_calc_cell_assembly_nonlinearity_from_cadata_stimulus(cadata, stimulus)
%     --------------------------------------------------------------------
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
%     stimulus : ripple stimulus in matrix format. Each row is one frequency,
%     and each column is one time bin.
% 
%     fio : struct holding the nonlinearities for teh spike train and the cell
%     assembly activations.
% 
%     fio.position = cadata.position;
%     fio.Patterns = cadata.Patterns;
%     fio.fsdvd = cadata.fsdvd;
%     fio.df = cadata.df;
%     fio.nf = nf;
%     fio.nlags = nlags;
%     fio.spk_fio = spk_fio;
%     fio.ca_fio = ca_fio;
%
%     The fio fields of fio are structs and have the following form:
% 
%     fio.Nbins = Nbins;
%     fio.sta = sta;
%     fio.x = x;
%     fio.px = px;
%     fio.pxspk = pxspk;
%     fio.pspk = pspk;
%     fio.pspkx = pspkx;
% 
%     fio.index_train = index_train_jack;
%     fio.sta_train = sta_train_jack;
%     fio.x_train = x_train_jack;
%     fio.px_train = px_train_jack;
%     fio.pxspk_train = pxspk_train_jack;
%     fio.pspk_train = pspk_train_jack;
%     fio.pspkx_train = pspkx_train_jack;




library('cell_assembly')
library('nonlinearity')
library('santosbox')
library('strfbox');
library('UtilitiesColormaps');


spktrain = cadata.spktrain; % spike train matrix for single units
nlags = cadata.nlags; % # time bins in STA
nf = cadata.nf; % # freq's in STA/Stimulus


% Threshold Activities to get cell assembly spiking
% ca_spktrain = # neurons X #trials
[ca_spktrain] = ca_threshold_cadata_activities(cadata);


% Nonlinearity over entire data
Nbins = 15;
for i = 1:size(spktrain,1) % go through each single unit spike train
    fprintf('Processing neuron #%.0f of %.0f\n', i, size(spktrain,1));
tic
    spk_fio(i) = ca_calc_fio_from_stimulus_spktrain(stimulus, spktrain(i,:), nlags, Nbins);
toc
end



for i = 1:size(ca_spktrain,1) % go through each cell assembly spike train
    fprintf('Processing assembly #%.0f of %.0f\n', i, size(ca_spktrain,1));
    ca_fio(i) = ca_calc_fio_from_stimulus_spktrain(stimulus, ca_spktrain(i,:), nlags, Nbins);
end


fio.position = cadata.position;
fio.Patterns = cadata.Patterns;
fio.fsdvd = cadata.fsdvd;
fio.df = cadata.df;
fio.nf = nf;
fio.nlags = nlags;
fio.spk_fio = spk_fio;
fio.ca_fio = ca_fio;


return;


















