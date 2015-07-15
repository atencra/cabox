function cafiofit = ca_calc_cell_assembly_nonlinearity_fit(cafio)
% ca_calc_cell_assembly_nonlinearity_fit  HVV nonlinearity curve fit to cell assembly nonlinearities
%
%     fiofit = ca_calc_cell_assembly_nonlinearity_fit(cafio)
%     --------------------------------------------------------------------
%     cafio : struct holding the nonlinearities for teh spike train and the cell
%     assembly activations. cafio is obtained from:
%     ca_calc_cell_assembly_sta_fio.m
% 
%     cafio has the following fields:
%
%             cafio.position
%             cafio.Patterns
%             cafio.fsdvd
%             cafio.df
%             cafio.nf
%             cafio.nlags
%             cafio.spk_fio
%             cafio.ca_fio
%
%     The fio fields of cafio are structs and have the following form:
% 
%             fio.Nbins = Nbins;
%             fio.sta = sta;
%             fio.x = x;
%             fio.px = px;
%             fio.pxspk = pxspk;
%             fio.pspk = pspk;
%             fio.pspkx = pspkx;
% 
%             fio.index_train = index_train_jack;
%             fio.sta_train = sta_train_jack;
%             fio.x_train = x_train_jack;
%             fio.px_train = px_train_jack;
%             fio.pxspk_train = pxspk_train_jack;
%             fio.pspk_train = pspk_train_jack;
%             fio.pspkx_train = pspkx_train_jack;

if ( isfield(cafio, 'cell_fio') )
    spk_fio = cafio.cell_fio; % nonlinearities for single neurons
else
    spk_fio = cafio.spk_fio; % nonlinearities for single neurons
end

ca_fio = cafio.ca_fio; % nonlinearities for cell assemblies
fsdvd = cafio.fsdvd; % sampling rate of audio system
df = cafio.df; % total downsampling factor of spectrotemporal envelope
dt = df / fsdvd; % temporal resolution of analysis

cafiofit = cafio; % reassign struct so that important fields are retained.

nCells = length(spk_fio); % # of single neurons
nAssemblies = length(ca_fio); % # of cell assemblies

fprintf('\nca_sta_fio_fit: spkfit\n');
spkfit = ca_sta_fio_fit(spk_fio, dt); % fit curves to single cell nonlinearities

fprintf('ca_sta_fio_fit: cafit\n');
cafit = ca_sta_fio_fit(ca_fio, dt); % fits to cell assembly nonlinearities

cafiofit.spk_fiofit = spkfit; % save to output struct
cafiofit.ca_fiofit = cafit; % save to output struct

return;




function fiostr = ca_sta_fio_fit(fio, dt)
% ca_sta_fio_fit Fit parametric curve to nonlinearity
% 
%     ca_sta_fio_fit(fio, dt)
% 
%     fio: struct holding nonlinearity data. fio has the following fields:
% 
%         Nbins
%         sta
%         x
%         px
%         pxspk
%         pspk
%         pspkx
%         index_train
%         sta_train
%         x_train
%         px_train
%         pxspk_train
%         pspk_train
%         pspkx_train
% 
%     dt: time bin size for which teh nonlinearity data was calculated, 
%     in seconds

for i = 1:length(fio)
    
    fprintf('fio #%.0f of %.0f\n', i, length(fio));

    % Entire spike train data
    x = fio(i).x;
    px = fio(i).px;
    pxspk = fio(i).pxspk;
    pspk = fio(i).pspk;
    pspkx = fio(i).pspkx;

    
    % Training set data
    x_train = fio(i).x_train;
    px_train = fio(i).px_train;
    pxspk_train = fio(i).pxspk_train;
    pspk_train = fio(i).pspk_train;
    pspkx_train = fio(i).pspkx_train;

    
    % Nonlinearity for the complete spike train
    fprintf('Global nonlinearity\n');
    fx = pspkx;
    fx = fx(:)';
    pspkx = fx;
    fx = fx ./ dt; % converty nonlinearity to firing rate
    [fiofit] = hvv_fio_fit(x, fx, pspkx, pspk, dt); % fit curve to nonlinearity

    fiofit_train = [];

%     % Nonlinearity for each training set nonlinearity
%     fprintf('Training set nonlinearity\n');
%     for ii = 1:length(x_train)
%         x = x_train{ii};
%         fx = pspkx_train{ii};
%         fx = fx(:)';
%         x = x(1:end-1); % last bin is most under-sampled, so ignore
%         fx = fx(1:end-1);
%         pspkx = fx;
%         pspk = pspk_train{ii};
%         fx = fx ./ dt; % converty nonlinearity to firing rate
%         [fiofit_train{ii}] = hvv_fio_fit(x, fx, pspkx, pspk, dt); % fit curve to nonlinearity
%     end % (for ii)

    fiostr(i).fiofit = fiofit;
    fiostr(i).fiofit_train = fiofit_train;

end





function [fiofit] = hvv_fio_fit(xbins, fx, pspkx, pspk, dt)
% hvv_fio_fit Hansel-van Vreewijk curve fit to nonlinearity
% 
%     [fiofit] = hvv_fio_fit(xbins, fx, pspkx, pspk, dt)
% 
%     xbins : normalized bins at which the nonlinearity was estimated.
% 
%     fx : nonlinearity, in spikes/s
% 
%     pspkx : nonlinearity, probability.
% 
%     pspk : probability of a spike in a bin.
% 
%     dt : temporal resolution, in seconds.
% 
%     fiofit : struct with curve fitting parameters.


% Need different optimization function - crashing for bad initial point ...
% Fit function from Ringach and Malone (2007)
a0 = [1 1 1]; % starting guess of optimization
[fitParams, resNorm] = lsqcurvefit(@ca_ringach_malone_func, a0, xbins, fx);

% a0 = [1 1 1 1]; % starting guess of optimization
% [fitParams, resNorm] = lsqcurvefit(@ringach_malone_func_with_baseline, a0, xbins, fx);


% interpolate to make it easier to find parameters
xFit = linspace(min(xbins), max(xbins), 100);
fxFit = ca_ringach_malone_func(fitParams, xFit);
% fxFit = ringach_malone_func_with_baseline(fitParams, xFit);


% Goodness of fit: normalized mean square error, r2 value
fxTemp = ca_ringach_malone_func(fitParams, xbins);
% fxTemp = ringach_malone_func_with_baseline(fitParams, xbins);
nmse = gfit2(fx, fxTemp, '2'); % gets goodness of fit measure, nmse
r2 = gfit2(fx, fxTemp, '8'); % gets goodness of fit measure, nmse


% Asymmetry of nonlinearity
right = pspkx(xbins>0.1);
left = pspkx(xbins<-0.1);
asi = ( sum(right) - sum(left) ) / ( sum(right) + sum(left) );   


% peak and skewness of nonlinearity
peakRate = max(fx);
skew = skewness(fx);

% find threshold of nonlinearity; this is the point where the nonlinearity
% exceeds the average firing rate of the neuron
index = find(xbins > -1);
bins_right = xbins(index);
pspkx_right = pspkx(index);
bins_right = bins_right(1:end-2);
pspkx_right = pspkx_right(1:end-2);
xi = linspace(min(bins_right), max(bins_right), 1000); % projection values
yi = interp1(bins_right, pspkx_right, xi, 'spline'); % interpolated nonlinearity

index = find(yi > pspk, 1);
threshold = xi(index);


clf;
hold on;
plot(xbins, fx, 'ko', 'markerfacecolor', 'k');
plot(xFit, fxFit, 'r-');
plot([fitParams(2) fitParams(2)], [0 max(fx)], 'k-');
xmax = max([max(xbins) abs(min(xbins))]);
xmax = xmax + 2*xmax*0.05;
xlim([-xmax xmax]);
ylimit = get(gca,'ylim');
set(gca,'ylim', [-0.1*max(fx) max(ylimit)]);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
legend('Data', 'HVV Fit', 'location', 'northwest');
title(sprintf('NMSE = %.3f, R2 = %.3f', nmse, r2));
pause(1);



% Save the data to a struct
fiofit.dt = dt;
fiofit.bins = xbins;
fiofit.pspkx = pspkx;
fiofit.pspk = pspk;
fiofit.x = xbins;
fiofit.fx = fx;
fiofit.xFit = xFit;
fiofit.fxFit = fxFit;
fiofit.fitParams = fitParams;
fiofit.nmse = nmse;
fiofit.r2 = r2;
fiofit.asi = asi;
fiofit.peakRate = peakRate;
fiofit.skew = skew;
fiofit.threshold = threshold;


return;






