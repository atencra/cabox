function [assem_train] = ca_get_cell_assembly_spktrain(cadata)
% ca_get_cell_assembly_spktrain Thresholded Cell Assembly Activations
% 
%     [assem_train] = ca_get_cell_assembly_spktrain(cadata)
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
%     assemb_train : matrix of threshold activations, same size as
%     cadata.Activities



%convert all spktrain to logical
spktrain = cadata.spktrain;
spktrain = spktrain>=1; 

%convert all ca_train to logical

activities = cadata.Activities;
dev = mad(activities,1,2);
med = median(activities,2);
thresh = med + 5*dev;
ca_train = zeros(size(activities));
assem_train = cell(size(activities,1),1);

for i = 1:length(thresh)
    ca_train(i,:) = activities(i,:)>=thresh(i);
    
    %get cell array for all assemblies and their members
    assem_train{i} = ca_train(i,:);
    
    for j = 1:length(cadata.assembly_members{i})
        idx = cadata.assembly_members{i}(j);
        assem_train{i} = [assem_train{i};spktrain(idx,:)];
    end
        
end 






