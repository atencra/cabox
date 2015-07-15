function fio = ca_calc_fio_from_stim_resp(stim, resp, Nbins)
% ca_sta_fio Calculate nonlinearity from STA
%
% fio = ca_sta_fio(stim, resp, Nbins, nf, nlags)
%
% stim : stimulus matrix, one trial per row. size(stim,2) = nf*nlags
% resp : vector of spike counts. length(resp) = size(stim,1)
% Nbins : # bins in the nonlinearity
% nf : # frequencies in stim matrix
% nlags : # time bins in stim matrix
%
% nf and nlags will also be the # of frequencies and time lags in the STA.
%
% fio : output struct holding the results


% Check and process input arguments
[nrows, ncols] = size(stim);
if ( ncols > nrows )
    error('Reshape stim to have one observation per row.');
end


[nrows, ncols] = size(resp);
if ( nrows > 1 && ncols > 1 )
    error('spikes should be a vector');
end

if ( nrows < ncols )
    resp = resp(:);
end



% STA over entire stimulus
sta = stim' * resp; %   (Ntrials X Ndims)' *(Ntrials X 1)
                    % = (Ndims X Ntrials) * (Ntrials X 1) = Ndims X 1
xraw = stim * sta; % all projection values
xprior = xraw ./ std(xraw); % scale to get units of standard deviation


% Nonlinearity over entire data
[x, px, pxspk, pspk, pspkx] = ca_calc_fio_from_stim_sta_proj(xprior, resp, Nbins);


% Now do the calculations for different training/test set data
r_perm = 1:size(stim,1); % number of rows is number of trials
Njackknifes = 4;


for n = Njackknifes:-1:1
    
    test_ind = r_perm(floor(size(stim,1)*(n-1)/Njackknifes)+1:floor(size(stim,1)*n/Njackknifes));
    %stim_test = stim(test_ind, :); % Don't need for this calculation
    %resp_test = resp(test_ind); % Don't need for this calculation

    train_ind = setdiff(r_perm, test_ind); % training set indices
    stim_train = stim(train_ind, :); % training set stimulus
    resp_train = resp(train_ind); % training set response


    % Estimate training set STA
    sta_train = stim_train'*resp_train; % training set STA
    xraw_train = stim_train * sta_train; % projection values for training set
    xprior_train = xraw_train ./ std(xraw_train); % scale to get units of standard deviation


    % Estimate training set nonlinearity
    [x_temp, px_temp, pxspk_temp, pspk_temp, pspkx_temp] = ...
        ca_calc_fio_from_stim_sta_proj(xprior_train, resp_train, Nbins);

    index_train_jack{n} = train_ind;
    sta_train_jack{n} = sta_train;
    x_train_jack{n} = x_temp;
    px_train_jack{n} = px_temp;
    pxspk_train_jack{n} = pxspk_temp;
    pspk_train_jack{n} = pspk_temp;
    pspkx_train_jack{n} = pspkx_temp;

    clear('stim_test', 'stim_train', 'resp_test', ...
    'resp_train', 'train_ind', 'sta_train', 'x_temp', ...
    'px_temp', 'pxspk_temp', 'pspk_temp', 'pspkx_temp');

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




