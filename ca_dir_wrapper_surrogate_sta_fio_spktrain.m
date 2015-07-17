function ca_dir_wrapper_surrogate_sta_fio_spktrain
% ca_batch_calc_cell_assembly_sta_stim_binsize Batch estimate PCA/ICA cell assemblies and STAs
% 
% ca_batch_calc_cell_assembly_sta_stim_binsize() looks inside the current
% folder and searches for files that end in *-strfcmb.mat or
% *-thresh-strf.mat. These files contain either spk/trigger/strf or
% thresh/trigger/strf variables.
%
% The spk and trigger variables are used to create binned spike trains.
% The spike trains are then processed by the cell assembly detection 
% algortithm of Vitor Santos. 
%
% The spk and trigger variables will have been attained after presenting a
% ripple noise stimulus. The stimuli may be rn1, rn4, rn8, or rn16. Each
% stimulus is approximately 10 minutes long.
%
% Additionally, the code will be run for different downsampling factors of
% the ripple files. The code assumes that each ripple noise envelope file
% will have been downsampled in time by a factor of 5, 10, 20, 40, and 100.
%
% Please remember that the 5-100 factors do not represent the total
% downsampling factor. The original ripple was created at 96kHz. The
% envelope file was not saved at this resolution. It was downsampled by a
% factor of 48. The 5-100 factors describe how the downsampled envelope
% file was further downsampled. To get the total downsampling factor of the
% stimulus, you multiply 48 times the factor: 48 * [5, 10, 20, 40, 100].
%
% The time resolution, in ms, of the binned spike trains will then be
%    48 *[5, 10, 20, 40, or 100] / 96 kHz * 1000
% 
% Note that this function will run the procedure on all *-strfcmb.mat and
% *-thresh-strf.mat files in the direction. Thus, it will process all
% ripple noise files.
% 
% At the end, the STAs for the binned spike trains, and for the cell
% assemblies, are computed and saved to files ending in *-ca-*dft.mat,
% where the * before the dft indicates the downsampling factor of the
% ripple envelope file, i.e. [5, 10, 20, 40, or 100].

library('santosbox');


d1 = dir('*-site*-*um-*db-rn16-fs*-*-strfcmb.mat'); % for single units
d2 = dir('*-site*-*um-*db-rn16-fs*-thresh-strf.mat'); % for multi-units

d = [d1(:)' d2(:)'];

% Downsampling factors of the .spr file. The total downsampling factor
% will be 48 * df(i), since the original stimulus was downsampled by
% a factor of 48 to achieve 0.5 ms resolution.
df = [10 20 40];
df = 10;

% Process every spk train file and every stimulus for every df(i) value.
for i = 1:length(d)

    fprintf('Processing %.0f of %.0f\n', i, length(d));

    filename = d(i).name;
    if ( isempty(findstr(filename,'thresh')) )
        str = load(filename, 'spk', 'trigger');
    else
        str = load(filename, 'thresh', 'trigger');
        str.spk = str.thresh;
    end
    rn = ca_get_ripple_noise_number_from_spk_filename(filename);
    stimtype = str.spk.stim;

    for j = 1:length(df)

        close all;

        outfile = sprintf('%s-ca-surrogate-%.0fdft.mat', filename(1:end-4), df(j));
        fprintf('Outfile: %s\n', outfile);

%         if ( exist(outfile, 'file') ) % skip if output already exists
%             fprintf('%s already processed.\n', outfile);
%         else
            [stimstr] = ca_get_ripple_noise_stimulus([], rn, df(j), []);
            [cadata] = ca_calc_surrogate_spktrain_from_sta_fio(...
                str.spk, stimstr.stimulus, str.trigger, stimstr.Fs, stimstr.DF);
            cadata.stimfile = stimstr.matrixfile;
            cadata.paramfile = stimstr.paramfile;
            cadata.faxis = stimstr.faxis;
            cadata.taxis = stimstr.taxis;
%             save(outfile, 'cadata');
%         end
    end % (for j)

end % (for i)


return;






