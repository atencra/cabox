function [spktrain_matrix, position] = ...
    ca_create_spktrain_matrix_from_stimulus_spk(spk, stimulus, trigger, ...
        FsDVD, totalDF)
% ca_create_spktrain_matrix_from_stimulus_spk Neuron spike time activity matrix
% 
%     [spktrain_matrix, position] = ...
%         ca_create_spktrain_matrix_from_stimulus_spk(spk, stimulus, trigger, ...
%               FsDVD, totalDF)
% 
%     Creates a matrix of binned spike trains. Each row represents one neuron,
%     and each column is a time bin. Spike times are obtained from the spk
%     struct array, and the number of rows is equal to the number of 
%     elements in spk.
% 
%     stimulus is a matrix, representing the ripple that was presented to the
%     neurons. It is Nf X Nt, where Nf is the number of frequencies and Nt is
%     the number of time bins.
% 
%     trigger is a vector of trigger pulse times, in sample number.
% 
%     FsDVD is the sampling rate of the DVD-Audio presentation system.
% 
%     totalDF is the total downsampling factor for stimulus. In most cases, 
%     the raw .spr file, a version downsampled in time is used. totalDF
%     describes how much the original ripple stimulus envelope was downsampled.
%     In usual cases, the original .spr file was downsampled 48 times.
%     Then, for STA and spike train analysis, this envelope file was
%     further downsampled by a factor of [5, 10, 20, 40, or 100].
%     Therefore, the totalDF would be 48 * [5, 10, 20, 40, or 100], which 
%     corresponds to 2.5, 5, 10, 20, and 50 ms bin sizes.
%
%     Remember bin size in ms = 48 * [5, 10, ..., 100] /96kHz * 1000.
%
%     spktrain_matrix is the matrix of binned spike times.
% 
%     position is the anatomical of each neuron in spk. Elements in position
%     correspond to rows in spktrain_matrix.


spktrain_matrix = zeros(length(spk), size(stimulus,2));

for i = 1:length(spk)
    spktimes = spk(i).spiketimes; % spike times from SpikeSort - in ms
    FsAD = spk(i).fs; % A/D system sampling rate
    spet = spktimes / 1000 * FsAD; % convert to sample number
    %spktimes = spktimes - trigger(1)/FsAD*1000; % align to trigger, convert to ms
    [locator] = ca_ripple_stim_spktrain(stimulus, spet, trigger, FsAD, FsDVD, totalDF);
    spktrain_matrix(i,:) = locator;
    position{i} = num2str(spk(i).position);
end % (for i)

return;



