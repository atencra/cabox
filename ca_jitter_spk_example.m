function ca_jitter_spk_example(spk,n)
% ca_jitter_spk_example Plot jittered spk spike train
%
%    ca_jitter_spk_example(spk,n)
%    
%    spk : struct array holding spike trains, one neurons
%    per spk element
%    
%    n : element in spk to plot. Default == 1.
%    
%    Takes the nth spike train in spk, and jitters individual
%    spike times by 25 ms. 15 surrogates are created, and then
%    plotted along with the original spike train.

narginchk(1,2);

if ( nargin == 1 )
    n = 1;
end
  
% create a random spike train (one 1-second trial, 25 spikes)
spktimes = spk(1).spiketimes; % in ms
spktimes = spktimes(:);


L = 2*.025; % jitter length (25 ms)
N = 15; % number of jitter resamples

L1000 = round(L*1000);


s1000 = SpikeJitter(N,spktimes,L1000,'centered','real');

allspikes = [spktimes(:)' s1000(:)'];
minmin = min(min(allspikes));
maxmax = max(max(allspikes));

% plot the rasters
figure;
hold on;
for i = 1:length(spktimes)
    plot([spktimes(i) spktimes(i)], [0 1.05*N], '-', 'color', 0.75*ones(1,3));
end % (for i)
plot([spktimes]',repmat([1.051*N],1,numel(spktimes)),'or','markersize', 2,'markerfacecolor', 'r')
plot([s1000]',repmat([(1:N)'],1,numel(spktimes)),'ok','markersize', 2,'markerfacecolor', 'k');
tickpref;
xlim([minmin maxmax]);
for i = 1:N, yt(i)=i;, ytlabel{i} = num2str(i);, end;
yt(end+1) = 1.051*N;
ytlabel{end+1} = 'Data';
set(gca,'ytick', yt, 'yticklabel', ytlabel);
xlabel('Time (ms)');
ylabel('Spike Train #');

return;


