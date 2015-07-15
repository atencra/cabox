function [CI] = ca_calc_ICA_threshold(spktrain,num_iter,alpha,opts)

% Get thresholds for cell assembly data to determine which units are
% 'considered part of an assembly'. Makes use of a random permutation test 
% to get the full range of ICA values by overriding the Marcenko-Pastur
% test and taking the first N PCs, where N = number of assemblies obtained
% from running the Marcenko-Pastur on the original dataset.
%
% Description of inputs:
%
%   spktrain:   Matrix of spike counts per time bin. Rows are single units;
%               columns are time bins.
%   num_iter:   Number of iterations for random permutation test. Default
%               set to 100.
%   alpha:      Confidence interval. For example, 95% returns 2.5 and 97.5
%               percentiles. Default set to 95.
%   opts:       Struct array containing options for assembly_patterns.m.
%               Optional.
%
%   Default number of iterations for fast_ICA for original data is 500.
%   However, this is 100 for the random permutations for speed. This 
%   can be changed in opts. 


if nargin == 1
    num_iter = 100;
    alpha = 95;
    opts = [];
elseif nargin == 2
    alpha = 95;
    opts = [];
elseif nargin == 3
    opts = [];
else
    error('Need at least one input');
end

%process original data
if isempty(opts)
    patterns = assembly_patterns(spktrain);
else
    patterns = assembly_patterns(spktrain,opts);
end

num_assemblies = size(patterns,2);

%for all the ICA values from permutation.
bigpermmat = zeros(size(patterns,1),size(patterns,2),num_iter);

for i = 1:num_iter
    
    fprintf('%d of %d iterations...\n',i,num_iter);
    spktrainperm = zeros(size(spktrain));
    N = size(spktrainperm,2);
    
    %permutation
    for j = 1:size(spktrain,1)
        X = rand(N,1);
        [~,idx] = sort(X);
        spktrainperm(j,:) = spktrain(j,idx);      
    end
    
    %run each iteration through the modified assembly_patterns code
    if isempty(opts)
        perm_patterns = assembly_patterns_permutation(spktrainperm, ...
            num_assemblies);
    else
         perm_patterns = assembly_patterns_permutation(spktrainperm, ...
            num_assemblies,opts);
    end
    
    bigpermmat(:,:,i) = perm_patterns;
    
end

%get lower and upper percentiles based on alpha
lowerprc = (100 - alpha)/2;
upperprc = 100 - lowerprc;
CI = prctile(bigpermmat(:),[lowerprc upperprc]);



    
