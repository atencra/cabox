function [ca_train] = ca_threshold_cadata_activities(cadata, ndev)
% ca_get_cell_assembly_spktrain Thresholded Cell Assembly Activations
% 
%     [assem_train] = ca_get_cell_assembly_spktrain(cadata, thresh)
% 
%     Takes cell assembly activations, stored in the struct cadata,
%     and returns threshold activations, where each activation
%     is the response of an assembly. The threshold is set to:
%     median + 5 mad.
% 
%     cadata : struct holding cell assembly analysis results. Usually
%     stored in a file such as:
%     
%     141215_234435-site6-800um-20db-rn1-fs20000-A-spk-strfcmb-ca-10dft.mat
%
%     ndev : optional deviation value to estimate threshold. The threshold
%     to determine an activity event is: thresh = med + ndev * dev. Default
%     value is ndev = 5.
% 
%     assemb_train : matrix of threshold activations, same size as
%     cadata.Activities

if ( nargin == 1 )
    ndev = 5;
end


%convert all spktrain to logical
spktrain = cadata.spktrain;
spktrain = spktrain >= 1; 

%convert all ca_train to logical

activities = cadata.Activities;
dev = mad(activities,1,2);
med = median(activities,2);
thresh = med + ndev * dev;
ca_train = zeros(size(activities));
assem_train = cell(size(activities,1),1);

for i = 1:length(thresh)
    ca_train(i,:) = activities(i,:)>=thresh(i); % convert event values to 0 or 1
end 






