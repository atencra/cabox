function ca_batch_calc_cell_assembly_nonlinearity_fit
% ca_batch_calc_cell_assembly_nonlinearity_fit Batch cell assembly nonlinearity curve fits
% 
% ca_batch_calc_cell_assembly_nonlinearity_fit() looks inside the current
% folder and searches for files that end in *-ca-*dft-fio.mat. These files 
% contain nonlinearities for the single unit and cell assembly STAs.
%
% This function then takes each nonlinearity and fits a parametric curve
% to it. The curve fit allows us to estimate the threshold and transition
% of the nonlinearity.
%
% The data are then save in files that end in *-fio-fit.mat


d = dir('*-site*-*um-*db-*-spk-strfcmb-ca-100dft-fio.mat'); % single units

if ( isempty(d) )
    fprintf('\nNo *-ca-*dft-fio.mat files in current directory.\n\n');
    return;
end



% Process every spk train file and every stimulus for every df(i) value.
for i = 3:3 %length(d)

    fprintf('\nProcessing %.0f of %.0f\n', i, length(d));

    filename = d(i).name;
    data = load(filename, 'cafio');
    outfile = sprintf('%s-fit.mat', filename(1:end-4) );

    fprintf('Infile: %s\n', filename);
    fprintf('Outfile: %s\n', outfile);

%     if ( exist(outfile, 'file') ) % skip if output already exists
%         fprintf('\n%s already processed.\n', outfile);
%     else

        cafiofit = ca_calc_cell_assembly_nonlinearity_fit(data.cafio);

        save(outfile, 'cafiofit');
        clear('cafio');
%     end

end % (for i)


return;

            
            
            