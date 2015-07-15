function ca_plot_fio_prob(fio)


spk_fio = fio.spk_fio;
ca_fio = fio.ca_fio;
nf = fio.nf;
nlags = fio.nlags;


for i = 1:length(spk_fio)

    sta = spk_fio(i).sta;
    x = spk_fio(i).x;
    pspkx = spk_fio(i).pspkx;
    pspk = spk_fio(i).pspk;

    sta_train = spk_fio(i).sta_train;
    x_train = spk_fio(i).x_train;
    pspkx_train = spk_fio(i).pspkx_train;
    pspk_train = spk_fio(i).pspk_train;

    figure;

    for ii = 1:length(sta_train)

        subplot(2, length(sta_train)+1, ii);
        imagesc(reshape(sta_train{ii},nf,nlags));
        tickpref;
        title(sprintf('Jackknife #%.0f', ii));

        subplot(2, length(sta_train)+1, ii + length(sta_train)+1 );
        plot(x_train{ii}, pspkx_train{ii}, 'ko-');
        tickpref;

    end


    subplot(2, length(sta_train)+1, length(sta_train)+1);
    imagesc(reshape(sta,nf,nlags));
    tickpref;
    title(sprintf('Total STA'));

    subplot(2, length(sta_train)+1, 2*(length(sta_train)+1) );
    plot(x, pspkx, 'ko-');
    tickpref;

pause

end % (for i)


return;


figure;

x = x .* 200;

subplot(2,3,1);
sd = std(x);
hist(x,50);
xlim([-1000 1000]);
ylim([0 10000]);
tickpref;
box off;
ylabel('Count');
title('Projection Values');

subplot(2,3,2);
hist(x ./ sd,50);
ylim([0 10000]);
tickpref;
box off;
title('Normalized Projection Values');

subplot(2,3,4);
hist(x(spikes>0),50)
xlim([-1000 1000]);
ylim([0 125]);
tickpref;
box off;
ylabel('Count');
xlabel('Projection Value (Arb Units)');

subplot(2,3,5);
hist(x(spikes>0) ./ sd,50);
ylim([0 125]);
tickpref;
box off;
xlabel('Projection Value (SD)');

subplot(233)
% plot_strf_symmetric_colormap(reshape(v,25,20))
hold on;
plot(centers, Pvx, 'ko-', 'markerfacecolor', 'k', 'markersize', 2);
plot(centers, Px_spike, 'o-', 'color', [0.6 0.6 0.6], ...
    'markerfacecolor', [0.6 0.6 0.6], 'markersize', 2);
plot(centers, Pspike .* Px_spike ./ Pvx, 'ro-', ...
    'markerfacecolor', 'r', 'markersize', 2)
plot([centers(1) centers(end)], [Pspike Pspike], 'k-');
ylim([0 0.25]);
tickpref;
title('Probability Distributions');
legend('P(x)', 'P(x|spike)', 'P(spike|x');

subplot(236)
hold on;
plot(centers, Pvx, 'ko-', 'markerfacecolor', 'k', 'markersize', 2);
plot(centers, Px_spike, 'o-', 'color', [0.6 0.6 0.6], ...
    'markerfacecolor', [0.6 0.6 0.6], 'markersize', 2);
plot(centers, Pspike .* Px_spike ./ Pvx, 'ro-', ...
    'markerfacecolor', 'r', 'markersize', 2)
plot([centers(1) centers(end)], [Pspike Pspike], 'k-');
ylim([0 0.25]);
tickpref;
title('Nonlinearity');
xlabel('Projection Value (SD)');


drawnow

set(gcf,'position', [1015         675         747         294]);


figure;
plot_strf_symmetric_colormap(reshape(v,25,20))
