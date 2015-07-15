function ca_batch_calc_cell_assembly_sta_fio_fit
% ca_batch_calc_cell_assembly_sta_fio_fit Parametric curve fit to cell assembly nonlinearities
% 
% ca_batch_calc_cell_assembly_sta_fio_fit() looks inside the current
% folder and searches for files that end in *-ca-*dft-fio.mat. These files 
% contain a cafio variable that holds the nonlinearities for the
% PCA/ICA cell assembly STA analysis.
%
% This function searches for all files in the current directory that end
% in *-*dft-fio.mat. Thus, it will process all cell assembly results for each
% stimulus and temporal downsampling factor.
%
% The function will then estimate the nonlinearities for each single unit
% and each cell assembly. The data will be stored in a struct called
% cafiofit, and then saved to a file ending in *dft-fio-fit.mat
%


d = dir('*-site*-*um-*db-*-spk-strfcmb-ca-*dft-fio.mat'); % single units


% Process every spk train file and every stimulus for every df(i) value.
for i = 1:1 %length(d)

    fprintf('Processing %.0f of %.0f\n', i, length(d));

    filename = d(i).name;
    data = load(filename, 'cafio');
    
    outfile = sprintf('%s-fio-fit.mat', filename(1:end-4) );

    fprintf('Outfile: %s\n', outfile);

    if ( exist(outfile, 'file') ) % skip if output already exists
        fprintf('%s already processed.', outfile);
    else
        
        cafiofit = ca_calc_cell_assembly_sta_fio_fit(cafio);
%         save(outfile, 'cafio');
    end

end % (for i)


return;

            
            








