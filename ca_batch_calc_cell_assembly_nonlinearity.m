function ca_batch_calc_cell_assembly_nonlinearity
% ca_batch_calc_cell_assembly_sta_fio Batch cell assembly nonlinearities
% 
% ca_batch_calc_cell_assembly_sta_fio() looks inside the current
% folder and searches for files that end in *-ca-*dft.mat. These files 
% contain a cadata variable that holds the results of the PCA/ICA 
% cell assembly analysis.
%
% This function searches for all files in the current directory that end
% in *-*dft.mat. Thus, it will process all cell assembly results for each
% stimulus and temporal downsampling factor.
%
% The function will then estimate the nonlinearities for each single unit
% and each cell assembly. The data will be stored in a struct called
% cafio, and then saved to a file ending in *dft-fio.mat
%


d = dir('*-site*-*um-*db-*-spk-strfcmb-ca-*dft.mat'); % single units

if ( isempty(d) )
    fprintf('\nNo *-ca-*dft.mat files in current directory.\n\n');
    return;
end



% Process every spk train file and every stimulus for every df(i) value.
for i = 1:length(d)

    fprintf('\nProcessing %.0f of %.0f\n', i, length(d));

    filename = d(i).name;
    data = load(filename, 'cadata');
    
    rn = ca_get_ripple_noise_number_from_spk_filename(filename);

    outfile = sprintf('%s-fio.mat', filename(1:end-4) );

    fprintf('Outfile: %s\n', outfile);

    if ( exist(outfile, 'file') ) % skip if output already exists
        fprintf('\n%s already processed.\n', outfile);
    else
        % Get temporal downsampling factor so we can use correct stimulus
        index1 = findstr(filename, 'ca-');
        index2 = findstr(filename, 'dft.mat');
        dft = str2num( filename( (index1+3):(index2-1) ) ); %#ok<*ST2NM>
        
        [stimstr] = ca_get_ripple_noise_stimulus([], rn, dft );
        cafio = ca_calc_cell_assembly_nonlinearity_from_cadata_stimulus(data.cadata, stimstr.stimulus); %#ok<*NASGU>
        save(outfile, 'cafio');
        clear('cafio', 'stimstr', 'dft');
    end

end % (for i)


return;

            
            
            