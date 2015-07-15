function ca_plot_cadata_stc_projection(sigstcdata, proj_exc)

nargchk(2,2,nargin);

library('nonlinearity');

y = filter(ones(1,20)/20,1,proj_exc);

t = (0:length(y)-1) * 5; % time in ms

evecs_exc = sigstcdata.evecs_exc;

figure;
subplot(2,1,1);
plot(t, abs(hilbert(proj_exc(:,1))));
ylim([-1 5]);
xlim([0 30000]);

subplot(2,1,2);
plot(t, abs(hilbert(proj_exc(:,3))));
ylim([-1 5]);
xlim([0 30000]);


% figure;
% subplot(2,1,1);
% plot(((proj_exc(:,1))));
% 
% subplot(2,1,2);
% plot(((proj_exc(:,2))));



return;




% 
% 
% 
% % Plot eigenvalue spectrum
% figure;
% hold on;
% plot([1-0.05*nft nft+0.05*nft], [0 0], 'k-');
% plot(evals_sort,'ko', 'markersize', 2); %, 'markerfacecolor', 'k');
% % plot((1:length(revals_sort)) / (length(revals_sort)/length(evals_sort)), revals_sort,'r-');
% plot([1-0.05*nft nft+0.05*nft], [max_reval max_reval], 'r-');
% plot([1-0.05*nft nft+0.05*nft], [min_reval min_reval], 'r-');
% xlim([1-0.05*nft nft+0.05*nft]);
% tickpref;
% ylim([minmin maxmax]);
% xlabel('Eigenvalue #');
% title(sprintf('#Ex: %.0f, #Inh: %.0f', nexc, ninh));
% 
% 
% return;
% 
% ca_plot_thresh_stc_evals_evecs(sta, evals, evecs, nf, nlags, evals_rand, '');
% 
% 
% return;
% 





function hf = plot_partial_distributions(stim, resp, prior_base, nfilt)
%PLOT_PARTIAL_DISTRIBUTIONS(stim,resp, prior_base) plots the prior and spike-triggered partial
%distributions along all possible couples of vectors of prior_base.
%
%
%STIM is the array of stimulus, each column corresponds to a stimulus
%dimension, each raw corresponds to a time
%
%RESP is the column array of spike times, with 1 for spikes, 0 elsewhere
%
%PRIOR_BASE is the matrix of vectors of the basis. Each raw corresponds to
%a stimulus dimension, each column corresponds to a vector of the base.


if size(stim,2)~=size(prior_base,1)
    error('stim and prior_base must be in the same space')
end

%project the prior and spike-triggered stimuli on the prior_base
s_prior = stim*prior_base;

spikes = find(resp);
s_spk = s_prior(spikes,:);

Nthres=size(prior_base,2);
if Nthres==2
    N=1;
    M=1;
elseif Nthres==3
    N=2;
    M=2;
elseif Nthres==4
    N=2;
    M=3;
elseif Nthres==5
    N=2;
    M=5;
elseif Nthres==6
    N=3;
    M=5;
elseif Nthres<=9
    N=3;
    M=6;
elseif Nthres<=20   
    N=4;
    M=6;
else
    error('prior_base has too many vectors, remove some of them');
end

figure;
set(gcf,'position',[20    22   953   693]);
print_mfilename(mfilename);
k=1;
for ii=1:Nthres-1
    for jj=ii+1:Nthres
        if k==N*M+1   %change window if there are too many plots to display
            k=1;
            figure;
            set(gcf,'position',[20    22   953   693]);
            print_mfilename(mfilename);
        end
        subplot(N,M,k);
        hold on;
        plot(s_prior(:,ii)',s_prior(:,jj)','bx');
        plot(s_spk(:,ii)',s_spk(:,jj)','rx');
        title(sprintf('X: %.0f , Y = %.0f', nfilt(ii), nfilt(jj) ));
        k=k+1;
    end
end




% % Plot Each filter and corresponding 1D projection distributions
% figure;
% for i = 1:size(prior_base,2)
% 
%    subplot(2,size(prior_base,2),i);
%    v = reshape(prior_base(:,i),25,20);
%    imagesc(v);
%    minmin = min(min(v));
%    maxmax = max(max(v));
%    boundary = max([abs(minmin) abs(maxmax)]);
%    set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
%    cmap = brewmaps('rdbu', 21);
%    colormap(cmap);
% 
%    subplot(2,size(prior_base,2),i+size(prior_base,2));
% 
%    prior_mn = mean(s_prior(:,i));
%    prior_sd = std(s_prior(:,i));
% 
%    prior = (s_prior(:,i) - prior_mn ) / prior_sd;
%    posterior = (s_spk(:,i) - prior_mn ) / prior_sd;
% 
%    edges = linspace(min(prior), max(prior), 15 );
%    [px] = histc(prior, edges);
%    [pxspk] = histc(posterior, edges);
%    pspkx = pxspk ./ px;
%    hold on;
%    plot(edges, pspkx, 'ko-');
% 
% end % (for i)
% print_mfilename(mfilename);
% 
return;





















% [evals] = diag(evals);
% [evals_sort,index] = sort(evals); % sort from low to high
% evecs_sort = evecs(:,index);
% 
% prior_base = evecs_sort(:,nfilt);
% 
% plot_partial_distributions(stim, resp, prior_base, nfilt);
% 
% 
% 
% figure;
% subplot(2,1,1);
% imagesc(Cspike)
% subplot(2,1,2);
% imagesc(Cprior)
% 
% size(priorstim)
% size(resp)
% 
% return;


% % STA over entire stimulus
% sta = stim' * resp; %   (Ntrials X Ndims)' *(Ntrials X 1)
%                     % = (Ndims X Ntrials) * (Ntrials X 1) = Ndims X 1
% xraw = stim * sta; % all projection values
% xprior = xraw ./ std(xraw); % scale to get units of standard deviation
% 
% 
% % Check for correct number of inputs
% narginchk(3,5);
% 
% close all;
% 
% for k = 1:length(threshstc)
% 
%    % Get STC results
%    nf = threshstc(k).nf; % #frequencies
%    nlags = threshstc(k).nlags; % #time lags
%    x0 = threshstc(k).x0;
% 
%    Cspike = threshstc(k).Cspike;
%    Cprior = threshstc(k).Cprior;
%    delC = Cspike-Cprior; % difference matrix
%    [evecs, evals] = eig(delC);
%    evals = diag(evals);
% 
%    % Get randomized STC results
%    revals = [];
%    for i = 1:length(threshstcrand(k).evals)
%       revals = [revals; threshstcrand(k).evals{i}];
%    end % (for i)
% %    rCspike = threshstcrand(k).Cspike;
% %    rCprior = threshstcrand(k).Cprior;
% %    rdelC = Cspike-Cprior; % difference matrix
% %    [revecs, revals] = eig(rdelC);
% %    revals = diag(revals);
% 
% 
%    % Specify which eigenvectors to plot
%    nft = nf * nlags;
%    nfilt = [nft:-1:(nft-5) 1 round(nft/2)];
% 
%    % Get STA filter
%    sta = threshsta(k).sta;
% 
%    str.exp = threshstc(k).exp;
%    str.chan = threshstc(k).chan;
%    str.position = threshstc(k).position;
%    str.stim = threshstc(k).stim;
% 
%    plot_thresh_stc_evals_evecs(sta, evals, evecs, nf, nlags, revals, str);
% 
%    if ( nargin == 4 )
% 
%       % Get stimulus
%       dsfile = threshloc(k).dsfile; % same stim file works for entire struct array
%       full_file_name = fullfile(stimfolder, dsfile);
%       load(full_file_name, 'stimulus'); % load stimulus matrix
% 
%       locator = threshloc(k).locator;
% 
%       % Parameters for the trial-based stimulus. Use all the
%       % frequencies and 20 of the time bins.
%       [stim, resp] = get_mne_stim_resp(stimulus, locator, nf, nlags, x0);
% 
%       [evals] = diag(evals);
%       [evals_sort,index] = sort(evals); % sort from low to high
%       evecs_sort = evecs(:,index);
% 
%       prior_base = evecs_sort(:,nfilt);
% 
%       plot_partial_distributions(stim, resp, prior_base, nfilt);
% 
%    end % (if)
% 
%    pause
%    close all;
% 
% end % (for k)
% 
% return;


