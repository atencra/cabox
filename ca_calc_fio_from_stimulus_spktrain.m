function fio = ca_calc_fio_from_stimulus_spktrain(stimulus, spktrain, nlags, Nbins)
% ca_calc_fio_from_stimulus_spktrain Calculate nonlinearity from STA
%
% fio = ca_calc_fio_from_stimulus_spktrain(stimulus, spktrain, nlags, Nbins)
%
% stimulus : stimulus matrix, one frequency per row, one time bin per column.
% spktrain: binned spike times. length(spktrain) == size(stimulus,2)
% nlags : # time bins in stim matrix
% Nbins : # bins in the nonlinearity
%
% fio : output struct holding the results


% Check and process input arguments
[nfreqs, ntimebins] = size(stimulus);

[nrows, ncols] = size(spktrain);
if ( nrows > 1 && ncols > 1 )
    error('spikes should be a vector');
end

if ( ntimebins ~= length(spktrain) )
    error('Need same number of time bins in stimulus and spktrain');
end




sta = ca_calc_sta_from_stimulus_spktrain(stimulus, spktrain, nlags);

[xprior, xposterior] = ca_calc_projection_from_sta_stimulus_spktrain(sta, stimulus, spktrain);

[x, px, pxspk, pspk, pspkx] = ca_calc_fio_prob_from_proj_spktrain(xprior, spktrain, Nbins);


% Now do the calculations for different training/test set data
time_bin_index = 1:size(stimulus,2); % number columns is number of time bins
Njackknifes = 4;


for n = Njackknifes:-1:1
    
    test_index = time_bin_index(floor(size(stimulus,2)*(n-1)/Njackknifes)+1:floor(size(stimulus,2)*n/Njackknifes));

    train_index = setdiff(time_bin_index, test_index); % training set indices
    stimulus_train = stimulus(:,train_index); % training set stimulus
    spktrain_train = spktrain(train_index); % training set response


    % Estimate training set STA
    sta_train = ca_calc_sta_from_stimulus_spktrain(stimulus_train, spktrain_train, nlags);

    [xprior_train, ~] = ca_calc_projection_from_sta_stimulus_spktrain(sta_train, stimulus_train, spktrain_train);

    [x_temp, px_temp, pxspk_temp, pspk_temp, pspkx_temp] = ...
        ca_calc_fio_prob_from_proj_spktrain(xprior_train, spktrain_train, Nbins);


    index_train_jack{n} = train_index;
    sta_train_jack{n} = sta_train;
    x_train_jack{n} = x_temp;
    px_train_jack{n} = px_temp;
    pxspk_train_jack{n} = pxspk_temp;
    pspk_train_jack{n} = pspk_temp;
    pspkx_train_jack{n} = pspkx_temp;

    clear('train_index', 'stimulus_train', 'spktrain_train', 'sta_train', ...
    'xprior_train', 'x_temp', 'px_temp', 'pxspk_temp', 'pspk_temp', 'pspkx_temp');

end % (for n)


% Save STA/nonlinearity results
fio.Nbins = Nbins;
fio.sta = sta;
fio.x = x;
fio.px = px;
fio.pxspk = pxspk;
fio.pspk = pspk;
fio.pspkx = pspkx;

fio.index_train = index_train_jack;
fio.sta_train = sta_train_jack;
fio.x_train = x_train_jack;
fio.px_train = px_train_jack;
fio.pxspk_train = pxspk_train_jack;
fio.pspk_train = pspk_train_jack;
fio.pspkx_train = pspkx_train_jack;

return;



