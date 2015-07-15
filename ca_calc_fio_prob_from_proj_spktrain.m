function [centers, Px, Pxspk, Pspk, Pspkx] = ca_calc_fio_prob_from_proj_spktrain(xprior, spktrain, Nbins)
% ca_calc_fio_prob_from_proj_spktrain Nonlinearity probabilities from projections
% 
%     [x, Px, Pxspk, Pspk, Pspkx] = ca_calc_fio_prob_from_proj_spktrain(xprior, spktrain, Nbins)
% 
%     xprior : projections of stimulus onto filter. Prior dist - w/o regard
%     to a spike.
%
%     spktrain : Ntrials X 1 vector of 0's and 1's. 
%     length(xprior) == length(spktrain).
%
%     Nbins : number of bins in the nonlinearity. 15 is standard.
% 
%     x : projection value bin centers in units of SD.
%     Px : probability distribution of projections without regard to a spike.
%     Pxspk : prob dist of projections with respect to a spike.
%     Pspk : prob of a spike. sum(spikes) / length(spikes)
%     Pspkx : nonlinearity. Prob(spike | projection)



% Check and process input arguments
if ( length(xprior) ~= length(spktrain) )
    error('#projections should equal length of spike train.');
end


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
Pxspk = zeros(1,Nbins);
for i = 1:Nbins
   Pxspk(i) = sum(spktrain(x_ind==i));
end
Pxspk = Pxspk(:);

Pxspk = Pxspk / sum(Pxspk); % prob(x|spk)
Pspk = sum(spktrain) ./ length(spktrain); % prob(spk)
Pspkx = Pspk * Pxspk ./ Px; % prob(spk|x) = prob(spk) * prob(x|spk) / prob(x)

return;



