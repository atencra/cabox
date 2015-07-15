function [locator] = ca_ripple_stim_spktrain(stimulus, spet, trigger, FsAD, FsDVD, DF)
% ca_ripple_stim_spktrain Binned spike train from stimulus, triggers.
%
%  [locator] = ca_ripple_stim_spktrain(stimulus, spet, trigger, FsAD, FsDVD, DF)
%  stimulus       :  NF X Ntrials matrix. A spectrogram with NF frequencies 
%                       and Ntrials time bins.
%  spet           :  array of spike times, in sample number.
%  trigger        :  Array of Trigger Times, in sample number
%  FsAD           :  Sampling Rate of A/D system.
%  FsDVD          :  Sampling rate of stimulus presentation system.
%  DF             : total downsampling factor of ripple stimulus. Usually
%                   48 * DFt.
%
%  locator    : spike train binned at the resolution of the stimulus.
%

numtrig = length(trigger);

% Setting 0 to be initial spike or trigger, whichever is smallest
mintime = min([trigger(:)' spet(:)']);
maxtime = max([trigger(:)' spet(:)']);
spet = spet - mintime + 1;
trigger = trigger - mintime + 1;

% Ntrials = NT*numtrig; % newnt=320;
Ntrials = size(stimulus, 2);
NT = Ntrials / length(trigger);
locator = zeros(Ntrials,1);

for trigcount = 2:length(trigger)-1 % skips spikes between 1st and 2nd trigger

    % Finding SPET in between triggers and resampling spet relative to the 
    index = find( spet>=trigger(trigcount) & spet<trigger(trigcount+1) );
    spettrig = ceil( (spet(index)-trigger(trigcount)+1) * FsDVD / FsAD / DF );

    for k = 1:length(spettrig)
        spike = spettrig(k);
        locator(spike+(trigcount-1)*NT) = locator(spike+(trigcount-1)*NT)+1;
    end
    
end % (while ~feof(fid) & trigcount<length(trigger)-1 )

% Closing all opened files
fclose('all');


% Error check locator length and stimulus length
%------------------------------------------------------------
% [n, nt, nf] = ripple_stim_length(specfile);
% 
% index_min = min( [length(locator) n/nf] );
% locator = locator( 1:index_min );
% locator_ipsi = locator_ipsi( 1:index_min );
% 
% 
% numspikes = sum(locator);
% stim_duration = ( max(trigger) - min(trigger) ) / FsAD;
% averate = numspikes / stim_duration;
% 
% % get time and frequency axes for the receptive field
% faxis = faxis;
% taxis = (-n1:n2-1) / (Fs/DF);


return;


