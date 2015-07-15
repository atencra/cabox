function [sta, itrain] = ca_get_sta_train_from_locator(locator, stimulus, numtbins)
% get_sta_train_from_locator - calculate STA from spike train vector and
%    a stimulus matrix. Do for 4 different training sets.
%
% [sta, itrain] = get_sta_train_from_locator(locator, stimulus, numtbins)
% -----------------------------------------------------------------------
%
% locator : vector of integers, where values greater
%           than one imply a spike, and values of 0
%           imply no spike
%
% stimulus : the entire ripple stimulus envelope file
%            as one matrix. Is contained in a .mat
%            file such as:
%
%           dmr-50flo-40000fhi-4SM-500TM-40db-48DF-21min_DFt2_DFf8-matrix.mat
%
% numtbins : number of time bins in sta. Default is 20.
%
% sta : 1x4 cell array. spike triggered average. The sta has the
%       same number of frequencies as the stimulus
%       matrix while the number of time bins is equal to numtbins.
%
% itrain : index for data that was used to construct the sta.
% 
% caa 4/12/11

[nf, ntrials] = size(stimulus);

if ( length(locator) ~= ntrials )
   error('Length of stimulus and length of locator must match.');
end


if ( isempty( numtbins ) )
   numtbins = 20;
end


[itrain1, itrain2, itrain3, itrain4] = ca_training_set_indices(ntrials);

% Put in cell array so we can loop easily
itrain{1} = itrain1;
itrain{2} = itrain2;
itrain{3} = itrain3;
itrain{4} = itrain4;


% Loop through each training set and get the sta
for n = 1:length(itrain)

   statemp = zeros(size(stimulus,1), numtbins);

   loc = locator( itrain{n} );
   stim = stimulus( :, itrain{n} );


   for i = numtbins:length(loc)
      if ( loc(i) )
         statemp = statemp + loc(i) * stim(:,i-numtbins+1:i);
      end
   end

   sta{n} = statemp; % save the sta to a cell array

end % (for n)



return;





