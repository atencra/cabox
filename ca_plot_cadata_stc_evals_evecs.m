function ca_plot_cadata_stc_evals_evecs(stcdata)

sta = stcdata.sta;
evals = stcdata.evals;
evecs = stcdata.evecs;
nf = stcdata.nf;
nlags = stcdata.nlags;
revals = stcdata.evals_rand;
str = '';

% Order the spike-triggered eigenvalues
[evals_sort,index] = sort(evals); % sort from low to high
evecs_sort = evecs(:,index); % order the eigenvectors

% Used later for plotting
erange = max(evals_sort) - min(evals_sort);
eminmin = min(evals_sort)-0.05*erange;
emaxmax = max(evals_sort)+0.05*erange;


% Order the random spike-triggered eigenvalues
[revals_sort,index] = sort(revals); % sort from low to high

rerange = max(revals_sort) - min(revals_sort);
reminmin = min(revals_sort)-0.05*rerange;
remaxmax = max(revals_sort)+0.05*rerange;

% max/min axis limits for later use
minmin = min([eminmin reminmin]);
maxmax = max([emaxmax remaxmax]);

[crit, ncrit] = randcrit(revals, 0.0001, 2);

% Determine the number of significant eigenvalues/eigenvectors

% % Determine significance level for different critical levels
% nrevals = length(revals_sort);
% nsig = ceil(nrevals / 1000);
% min_reval = revals_sort(nsig);
% max_reval = revals_sort(end-nsig);

% Most conservative approach - only count as significant if evals exceed
% every random eval
min_reval = totalmin(revals_sort);
max_reval = totalmax(revals_sort);


% index_exc = find(evals_sort > max_reval);
% index_inh = find(evals_sort < min_reval);

index_exc = find(evals_sort > crit(2) );
index_inh = find(evals_sort < crit(1) );


nexc = length(index_exc);
ninh = length(index_inh);
nfilt = [sort(index_exc(:)',2,'descend') sort(index_inh(:)',2,'descend') round(length(evals_sort)/2)];

Nthres = nexc + ninh + 4; % num of significant, plus plots for sta, 
                          % eigenvalues, and randomized eigenvale distribution
if Nthres==2
    N=2;
    M=3;
elseif Nthres==3
    N=2;
    M=3;
elseif Nthres==4
    N=2;
    M=3;
elseif Nthres==5
    N=2;
    M=3;
elseif Nthres==6
    N=2;
    M=3;
elseif Nthres>6 && Nthres <= 9
    N=3;
    M=3;
elseif Nthres>9 && Nthres <= 12
    N=4;
    M=3;
elseif Nthres>12 && Nthres <= 16
    N=4;
    M=4;
elseif Nthres>16 && Nthres <= 20
    N=5;
    M=4;
else
    error('Too many vectors to plot. Remove some of them');
end

nft = nf * nlags; % total number of dimensions


figure;

h = fspecial('gaussian',2,2);

subplot(N,M,1);
sta = reshape(sta, nf, nlags);
rfsig = imfilter(sta, h);
plot_strf_symmetric_colormap(rfsig);
set(gca, 'xtick', [], 'xticklabel', '');
set(gca, 'ytick', [], 'yticklabel', '');
cmap = brewmaps('rdbu', 21);
colormap(cmap);
title('STA');


% Plot eigenvalue spectrum
subplot(N,M,2);
hold on;
plot([1-0.05*nft nft+0.05*nft], [0 0], 'k-');
plot(evals_sort,'ko', 'markersize', 2); %, 'markerfacecolor', 'k');
plot((1:length(revals_sort)) / (length(revals_sort)/length(evals_sort)), revals_sort,'r-');
plot([1-0.05*nft nft+0.05*nft], [max_reval max_reval], 'r-');
plot([1-0.05*nft nft+0.05*nft], [min_reval min_reval], 'r-');
xlim([1-0.05*nft nft+0.05*nft]);
tickpref;
ylim([minmin maxmax]);
xlabel('Eigenvalue #');
title(sprintf('#Ex: %.0f, #Inh: %.0f', nexc, ninh));


% Plot eigenvalue spectrum
subplot(N,M,3);
hold on;
edges = linspace(min(revals_sort),max(revals_sort), 100);
count = histc(revals_sort, edges);
barh(edges, count, 'histc');
plot([0 max(count)], [max_reval max_reval], 'r-');
plot([0 max(count)], [min_reval min_reval], 'r-');
xlim([0 1.05*max(count)]);
ylim([minmin maxmax]);
tickpref;
title('Rand Eigenvalues');


 % Plot eigenvectors specified in nfilt
 for i = 1:length(nfilt)
    subplot(N,M,i+3);
    evec = reshape(evecs(:,nfilt(i)),nf,nlags);
    plot_strf_symmetric_colormap(evec);
    set(gca, 'xtick', [], 'xticklabel', '');
    set(gca, 'ytick', [], 'yticklabel', '');
    tickpref;
    colormap(cmap);
    title(sprintf('E%.0f = %.2f', nfilt(i), evals(nfilt(i))));
 end % (for i)


if ( isfield(stcdata, 'exp') )
    suptitle(sprintf('%s site%.0f chan%.0f model%.0f %.0f um', ...
        stcdata.exp, stcdata.site, stcdata.chan, stcdata.model, stcdata.position));
end

set(gcf,'position', [80   188   581   811]);

% suptitle(sprintf('%s chan%.0f %.0fum %s', str.exp, str.chan, str.position, str.stim));

return;


