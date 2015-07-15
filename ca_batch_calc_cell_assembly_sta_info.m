function ca_batch_calc_cell_assembly_sta_info
% ca_batch_calc_cell_assembly_sta_info Batch cell assembly STA information
% 
% ca_batch_calc_cell_assembly_sta_info() looks inside the current
% folder and searches for files that end in *-ca-*dft.mat. These files 
% contain a cadata variable that holds the results of the PCA/ICA 
% cell assembly analysis.
%
% This function searches for all files in the current directory that end
% in *-*dft.mat. Thus, it will process all cell assembly results for each
% stimulus and temporal downsampling factor.
%
% The function will then estimate the STA information for each single unit
% and each cell assembly. The data will be stored in a struct called
% cainfo, and then saved to a file ending in *dft-info.mat
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

    outfile = sprintf('%s-info.mat', filename(1:end-4) );

    fprintf('Outfile: %s\n', outfile);

    if ( exist(outfile, 'file') ) % skip if output already exists
        fprintf('%s already processed.', outfile);
    else
        % Get temporal downsampling factor so we can use correct stimulus
        index1 = findstr(filename, 'ca-');
        index2 = findstr(filename, 'dft.mat');
        dft = str2num( filename( (index1+3):(index2-1) ) ); %#ok<*ST2NM>
        
        [stimstr] = ca_get_ripple_noise_stimulus([], rn, dft );
        [cainfo] = ca_get_cell_assembly_sta_info(data.cadata, stimstr.stimulus);
        save(outfile, 'cainfo');
        clear('cainfo', 'data', 'stimstr', 'dft');
    end

end % (for i)


return;

            
            
            