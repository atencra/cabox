function assembly_members = ca_identify_assembly_members(cadata)

% Identifies assembly members based on thresholds calculated in
% ca_calc_ICA_threshold.m.
%   
%   cadata: must include cadata.CI

if ~isfield(cadata,'CI')
    error('Please calculate threshold first');
end

CI = cadata.CI;
patterns = cadata.Patterns;
assembly_members = cell(size(patterns,2),1);

for i = 1:size(patterns,2)
    
    idx1 = find(patterns(:,i) <= CI(1));
    idx2 = find(patterns(:,i) >= CI(2));
    assembly_members{i} = sort([idx1; idx2]);

end

