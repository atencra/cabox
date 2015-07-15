function [centers, Px, Pxspk, Pspk, Pspkx] = ca_proj2fio(xprior, xposterior, Nbins)
% stim_resp_sta_fio Nonlinearity from stimulus trial matrix and spike train vector
% 
%     [x, Px, Pxspk, Pspk, Pspkx] = stim_resp_sta_fio(v, stim, spikes, Nbins)
% 
%     v : filter. Usually the STA. A matrix: #freq's X #time lags
%     stimulus : stimulus matrix. #freq's X #time bins
%     spikes : Ntrials X 1 vector of 0's and 1's.
%     Nbins : number of bins in the nonlinearity. 15 is standard.
% 
%     x : projection value bin centers in units of SD.
%     Px : probability distribution of projections without regard to a spike.
%     Pxspk : prob dist of projections with respect to a spike.
%     Pspk : prob of a spike. sum(spikes) / length(spikes)
%     Pspkx : nonlinearity. Prob(spike | projection)



%Compute histogram edges
xmin = min(xprior);
xmax = max(xprior);
step = (xmax-xmin)/Nbins;

edges = xmin:step:xmax;
centers = ( edges(1:end-1) + edges(2:end) ) / 2;

%Compute Px
[Px,x_ind] = histc(xprior, edges);
Px(end-1) = Px(end-1) + Px(end); % histc puts all entries that match the 
Px(end) = [];                    % upper limit of the last bin in an extra category
x_ind(x_ind==Nbins+1) = Nbins;
Px = Px / sum(Px); 
Px = Px(:);


%Compute Pv(x|spike)
[Pxspk, x_ind] = histc(xposterior, edges);
Pxspk(end-1) = Pxspk(end-1) + Pxspk(end); % histc puts all entries that match the 
Pxspk(end) = [];                    % upper limit of the last bin in an extra category
x_ind(x_ind==Nbins+1) = Nbins;
Pxspk = Pxspk / sum(Pxspk); 
Pxspk = Pxspk(:);

Pspk = length(xposterior) ./ length(xprior); % prob(spk)

Pspkx = Pspk * Pxspk ./ Px; % prob(spk|x) = prob(spk) * prob(x|spk) / prob(x)



return;











% figure;
% 
% x = x .* 200;
% 
% subplot(2,3,1);
% sd = std(x);
% hist(x,50);
% xlim([-1000 1000]);
% ylim([0 10000]);
% tickpref;
% box off;
% ylabel('Count');
% title('Projection Values');
% 
% subplot(2,3,2);
% hist(x ./ sd,50);
% ylim([0 10000]);
% tickpref;
% box off;
% title('Normalized Projection Values');
% 
% subplot(2,3,4);
% hist(x(spikes>0),50)
% xlim([-1000 1000]);
% ylim([0 125]);
% tickpref;
% box off;
% ylabel('Count');
% xlabel('Projection Value (Arb Units)');
% 
% subplot(2,3,5);
% hist(x(spikes>0) ./ sd,50);
% ylim([0 125]);
% tickpref;
% box off;
% xlabel('Projection Value (SD)');
% 
% subplot(233)
% % plot_strf_symmetric_colormap(reshape(v,25,20))
% hold on;
% plot(centers, Pvx, 'ko-', 'markerfacecolor', 'k', 'markersize', 2);
% plot(centers, Px_spike, 'o-', 'color', [0.6 0.6 0.6], ...
%     'markerfacecolor', [0.6 0.6 0.6], 'markersize', 2);
% plot(centers, Pspike .* Px_spike ./ Pvx, 'ro-', ...
%     'markerfacecolor', 'r', 'markersize', 2)
% plot([centers(1) centers(end)], [Pspike Pspike], 'k-');
% ylim([0 0.25]);
% tickpref;
% title('Probability Distributions');
% legend('P(x)', 'P(x|spike)', 'P(spike|x');
% 
% subplot(236)
% hold on;
% plot(centers, Pvx, 'ko-', 'markerfacecolor', 'k', 'markersize', 2);
% plot(centers, Px_spike, 'o-', 'color', [0.6 0.6 0.6], ...
%     'markerfacecolor', [0.6 0.6 0.6], 'markersize', 2);
% plot(centers, Pspike .* Px_spike ./ Pvx, 'ro-', ...
%     'markerfacecolor', 'r', 'markersize', 2)
% plot([centers(1) centers(end)], [Pspike Pspike], 'k-');
% ylim([0 0.25]);
% tickpref;
% title('Nonlinearity');
% xlabel('Projection Value (SD)');
% 
% 
% drawnow
% 
% set(gcf,'position', [1015         675         747         294]);
% 
% 
% figure;
% plot_strf_symmetric_colormap(reshape(v,25,20))



