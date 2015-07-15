function [ifrac] = ca_info_fraction(xprior, xposterior, frac, nreps)
%ca_info_fraction - filter information for different data fractions
%
% [ifrac] = ca_info_fraction(locator, x, xbins, frac, nreps)
% -------------------------------------------------------------
%
% x : projection values onto a filter. A vector.
%
% xbins : center of bins at which projection probability distributions
% will be calculated
%
% frac : data fraction for which information will be estimated. A scalar, 
% usually something like 80 or 92.5
%
% ifrac : data fraction information estimates. 1x5 vector. 5 estimates 
% 
% caa 3/10/09

if ( nargin == 3 )
   nreps = 5;
end


frac = frac / 100;

ntrials = length(xposterior);

ntrials_subset = round(frac * ntrials); % number of reduced trials

ifrac = zeros(1,nreps);

rand('twister',sum(100*clock));

% Now normalize by the mean and std of the prior
xmn = mean(xprior);
xstd = std(xprior);

xprior_scaled = (xprior - xmn) ./ xstd;

xbins_edges = linspace(min(xprior_scaled), max(xprior_scaled), 15);
xbins_centers = edge2center(xbins_edges);

% length(xposterior)
% ntrials
% ntrials_subset

for m = 1:nreps

   xspktemp = randsample(xposterior, ntrials_subset); % get a subset of the posterior distribution


   % Now normalize by the mean and std of the prior
   xspktemp = (xspktemp - xmn) ./ xstd;


   % Form the probability distributions
   nx = hist(xprior_scaled, xbins_centers);
   px = nx ./ sum(nx); % p(x)
   px = px(:);

   nxspk = hist(xspktemp, xbins_centers);
   pxspk = nxspk ./ sum( nxspk ); % p(x|spk)
   pxspk = pxspk(:);

   % Get and assign the information
   iplugin = ca_info_px_pxspk(px, pxspk);

   ifrac(m) = iplugin;

% [frac length(xspktemp) ntrials iplugin]
% 
% pause

end % (for m)

return;

