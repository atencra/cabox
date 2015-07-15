function [info_mean, info_std] = ca_info_fraction_mean_std(fraction, ifraction)
%ca_info_fraction_mean_std - mean, standard dev of information data fraction calculations.
%
% [info_mean, info_std] = ca_info_fraction_mean_std(fraction, ifraction)
%
% fraction : proportions of the data used. A vector like [80 90 95 97.5
% 100], or something similar.
%
% info_fraction : cell array for the data fractions in fraction. Each 
% element of the cell array usually holds multiple information estimates 
% for each data fraction(i).
%
% info_mean : mean information values for each data fraction(i). A vector.
%
% info_std : standard deviation of information values for each data fraction(i). A vector.
%

info_mean = zeros(size(fraction));
info_std = zeros(size(fraction));

for i = 1:length(fraction)
   temp = ifraction{i};
   info_mean(i) = mean(temp);
   info_std(i) = std(temp);
end % (for i)


return;

