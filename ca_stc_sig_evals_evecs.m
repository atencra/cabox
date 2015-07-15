function sigstcdata = ca_stc_sig_evals_evecs(stcdata)

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


% Order the random spike-triggered eigenvalues
[revals_sort,index] = sort(revals); % sort from low to high

[crit, ncrit] = randcrit(revals, 0.0001, 2);


% Determine the number of significant eigenvalues/eigenvectors

% index_exc = find(evals_sort > max_reval);
% index_inh = find(evals_sort < min_reval);

index_exc = find(evals_sort > crit(2) );
index_inh = find(evals_sort < crit(1) );

evals_exc = evals_sort(index_exc);
evals_exc = fliplr(evals_exc(:)');
evecs_exc = fliplr(evecs(:,index_exc));


evals_inh = evals_sort(index_inh);
evals_inh = fliplr(evals_inh(:)');
evecs_inh = fliplr(evecs(:,index_inh));


nexc = length(index_exc);
ninh = length(index_inh);
nfilt = [sort(index_exc(:)',2,'descend') sort(index_inh(:)',2,'descend') round(length(evals_sort)/2)];

nft = nf * nlags; % total number of dimensions


sigstcdata.nf = nf;
sigstcdata.nlags = nlags;
sigstcdata.sta = sta;
sigstcdata.evals = evals;
sigstcdata.crit = crit;
sigstcdata.evals_exc = evals_exc;
sigstcdata.evecs_exc = evecs_exc;
sigstcdata.evals_inh = evals_inh;
sigstcdata.evecs_inh = evecs_inh;


return;


