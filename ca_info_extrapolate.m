function [info_extrap] = ca_info_extrapolate(fraction, info_mean)
% ca_info_extrapolate Asymptotic info estimate from mean info values
% 
% [info_extrap] = ca_info_extrapolate(fraction, info_mean)
% 
% We fit a line to the information values versus the inverse of the data
% fraction. The y-intercept of the line represents the information value
% that would be present for infinite data set size. This y-intercept value
% is the extrapolated information value. 
% 
% fraction : the fractions of the data that were used. Usually something
% like [80 85 90 92.5 95 97.5 100]
% 
% info_mean : a vector of mean information values for each data fraction.
% 
% info_extrap : extrapolated information value.


info_extrap = 0;
x = 1./fraction;
beta = polyfit(x,info_mean,1);
info_extrap = beta(2);

return;

