function ca_error_check_fio(cadata, stimstr)


stimulus = stimstr.stimulus;
spktrain = cadata.spktrain(2,:);
nlags = cadata.nlags;
nf = cadata.nf;


Nbins = 15;


sta1 = ca_get_sta_from_locator(spktrain, stimulus, nlags);
tic
[xprior1, xposterior1] = ca_get_sta_projection_from_stimulus_matrix(sta1, spktrain, stimulus);
[x1, px1, pxspk1, pspk1, pspkx1] = ca_proj2fio(xprior1, xposterior1, Nbins);
toc

[stim, resp] = ca_create_stim_trial_from_stim_matrix(stimulus, spktrain, nlags);
tic
sta2 = stim' * resp';
[x2, px2, pxspk2, pspk2, pspkx2] = stim_resp_sta_fio(sta2, stim, resp, Nbins);
toc

xraw2 = stim * sta2; % all projection values
xprior2 = xraw2 ./ std(xraw2); % scale to get units of standard deviation



% size(stim)
% size(xprior1)
% size(xprior2)

% sum(sum(xprior1-xprior2))
% max(xprior1)
% max(xprior2)


figure;
subplot(2,1,1);
imagesc(sta1);

subplot(2,1,2);
imagesc(reshape(sta2,nf, nlags));

figure;
subplot(2,1,1);
hist(xprior1);

subplot(2,1,2);
hist(xprior2);


figure;
plot(sort(xprior1), sort(xprior2),'.');



figure;
hold on;
plot(x1, pspkx1, 'ro-');
plot(x2, pspkx2, 'ko', 'markerfacecolor', 'k');







return;

