function [irep1, irep2, irep3, irep4] = ca_training_set_indices(ntrials)
% ca_training_set_indices - indices for training set repetition 
%
% [irep1, irep2, irep3, irep4] = ca_training_set_indices(ntrials)
% --------------------------------------------------------------
%
% ntrials : total number of trials in data set. Usually a number like
% 1115226, or something similar and of that magnitude.
%
% We separate data into training and test sets. Training sets are made
% up from 3/4 of the data. We can divide the data four ways:
%
% training set 1 = [2 3 4]
% training set 2 = [1 3 4]
% training set 3 = [1 2 4]
% training set 4 = [1 2 3]
%
% where 1, 2, 3, 4 represent the first 1/4, the second 1/4, the third 1/4,
% and the fourth 1/4 of the data set.
%
% The test set then is comprised of 1/4 of the data. We can notate this as
%
% test set 1 = [1]
% test set 2 = [2]
% test set 3 = [3]
% test set 4 = [4]
%
% So this function returns the indices, irep1, ..., irep4, 
% into the 4 different training sets.
%
% caa 3/12/09

ntrials_testrep = floor( ntrials / 4 );

index1 = [1:ntrials_testrep];
index2 = [ntrials_testrep+1:2*ntrials_testrep];
index3 = [2*ntrials_testrep+1:3*ntrials_testrep];
index4 = [3*ntrials_testrep+1:ntrials];

irep1 = [index2 index3 index4];
irep2 = [index1 index3 index4];
irep3 = [index1 index2 index4];
irep4 = [index1 index2 index3];

return;
