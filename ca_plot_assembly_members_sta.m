function ca_plot_assembly_members_sta(cadata, nf, nlags)

% Plots the STA of each cell assembly with the STAs of its member neurons 
% as determined via the permutation tests in ca_calc_ICA_threshold. 
%
%   cadata:     needs to have cadata.assembly_members. To obtain that, get
%               threshold from ca_calc_ICA_threshold.m first, and then run
%               ca_identify_assembly_members.m.
%   
%   nf:         #frequencies in the STA
% 
%   nlags:      #time bins in the STA        

if ~isfield(cadata,'assembly_members')
    error('Please calculate threshold and identify members of cell assembly first');
end

ca_stamat = cadata.ca_stamat;
stamat = cadata.stamat;

for i = 1:size(ca_stamat,1)

    hf = figure;
   
    num_mem = length(cadata.assembly_members{i});
    
    if num_mem == 1
        nrows = 2;
        ncols = 1;
    elseif num_mem == 2
        nrows = 2;
        ncols = 2;
    elseif num_mem <= 4
        nrows = 3;
        ncols = 2;
    elseif num_mem <= 6
        nrows = 3;
        ncols = 3;
    elseif num_mem <= 9
        nrows = 4;
        ncols = 3;
    elseif num_mem <= 12
        nrows = 4;
        ncols = 4;
    elseif num_mem <= 16
        nrows = 5;
        ncols = 4;
    elseif num_mem <= 20
        nrows = 6;
        ncols = 5;
    elseif num_mem <= 25
        nrows = 7;
        ncols = 5;
    else
        error('Too many members right now.');
    end
      
    subplot(nrows, ncols, ceil(ncols/2));
    rfmat = reshape(ca_stamat(i,:), nf, nlags);
    plot_strf_symmetric_colormap(rfmat);
    title(sprintf('CA #%.0f',i));
    
    for j = 1:length(cadata.assembly_members{i})
        subplot(nrows, ncols, ncols+j);
        rfmat = reshape(stamat(cadata.assembly_members{i}(j),:), nf, nlags);
        plot_strf_symmetric_colormap(rfmat);
        title(sprintf('Neuron #%.0f',cadata.assembly_members{i}(j)));
    end
    
    
end % (for i)



% % Plot STAs for the cell assemblies.
% fprintf('Single Neuron STA:\n');
% hf = figure;
% for i = 1:size(stamat,1)
%     clf(hf);
%     rfmat = reshape(stamat(i,:), nf, nlags);
%     plot_strf_symmetric_colormap(rfmat);
%     pause
% end % (for i)



end

