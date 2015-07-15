function sta = ca_calc_sta_from_stimulus_spktrain(stimulus, spktrain, numtbins)
% ca_get_sta_from_stimulus_spktrain STA from binned spike train vector and stimulus matrix
%
% sta = ca_calc_sta_from_stimulus_spktrain(stimulus, spktrain, numtbins)
% ----------------------------------------------------------------------
%
% spktrain : vector of integers, where values greater than one imply a 
% spike, and values of 0 imply no spike. This is a binned spike train.
%
% stimulus : the entire ripple stimulus envelope file as one matrix. Is 
% contained in a .mat file such as:
%
%           dmr-50flo-40000fhi-4SM-500TM-40db-48DF-21min_DFt2_DFf8-matrix.mat
%
% sta : spike triggered average. The sta has the same number of frequencies 
% as the stimulus matrix while the number of time bins is 20.
%
% numtbins : how much memory to include in the sta. Default is 20 time
% bins. numtbins is the same as nlags used in other programs.
% 


if ( nargin == 2 )
   numtbins = 20;
end

if ( size(stimulus,2) ~= length(spktrain) )
    error('Stimulus and spike train must have same number of time bins.');
end

sta = zeros(size(stimulus,1), numtbins);


for i = numtbins:length(spktrain)
   if ( spktrain(i) )
      sta = sta + spktrain(i) * stimulus(:,i-numtbins+1:i);
   end
end


return;





