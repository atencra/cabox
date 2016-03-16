function ca_plot_cell_assembly_stamat(stamat, nf, nlags, datatype)
% ca_plot_cell_assembly_sta Cell assembly STA plot
% 
%     ca_plot_cell_assembly_sta(ca_stamat, nf, nlags, datatype)
% 
%     stamat or ca_stamat : stamat is a matrix of STAs for individual
%     neurons, while ca_stamat is matrix of cell assembly STAs. 
%     Each row is one STA.
% 
%     nf : #frequencies in the STA
% 
%     nlags : #time bins in the STA
%
%     datatype : 'STA' for single neuron STAs, or 'CA' for cell assembly
%     data. Default is 'CA'.
%
%     stamat and ca_stamat may be obtained from:
%
%       [cadata] = ca_calc_cell_assembly_sta(spk, stimulus, trigger, FsDVD, totalDF);
%
%     stamat and ca_stamat are fields of cadata.
%

library('strfbox');

narginchk(3,4);

if ( nargin == 3 )
    datatype = 'CA';
end

if ( nargin == 4 )
    if ( isempty(datatype) )
        datatype = 'CA';
    elseif ( ~strcmp(datatype, 'STA') )
        datatype = 'CA';
    end

end


ncells = size(stamat,1);

if ncells == 1
    nrows = 1;
    ncols = 1;
elseif ncells == 2
    nrows = 1;
    ncols = 2;
elseif ncells <= 4
    nrows = 2;
    ncols = 2;
elseif ncells <= 6
    nrows = 2;
    ncols = 3;
elseif ncells <= 9
    nrows = 3;
    ncols = 3;
elseif ncells <= 12
    nrows = 3;
    ncols = 4;
else % ncells <= 16
    nrows = 4;
    ncols = 4;
end



if ( strcmp(datatype, 'CA') ) % Plot STAs for the cell assemblies.
    fprintf('Cell Assembly STAs:\n');
else
    fprintf('Single Neuron STAs:\n');
end


figure;
n = 1;
for i = 1:size(stamat,1)
    subplot(nrows, ncols, n);
    rfmat = reshape(stamat(i,:), nf, nlags);
    plot_strf_symmetric_colormap(rfmat);
    if ( strcmp(datatype, 'CA') ) % Plot STAs for the cell assemblies.
        title(sprintf('CA #%.0f',i));
    else
        title(sprintf('STA #%.0f',i));
    end
    n = n + 1;
    if ( n > 16 )
        if ( mod(i, size(stamat,1)) )
            figure;
            n = 1;
        end
    end
end % (for i)



return;


