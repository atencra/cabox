function rn = ca_get_ripple_noise_number_from_spk_filename(filename)
% ca_get_ripple_noise_number_from_spk_filename Ripple Noise # from spk filename
% 
%     rn = ca_get_ripple_noise_number_from_spk_filename(filename)
% 
%     Determines which ripple noise is indicated in filename.
% 
%     Filenames have the form: 
%
%         141215_222459-site5-800um-20db-rn1-fs20000-A-spk-strfcmb.mat
%
%     or
%
%        141215_222459-site5-800um-20db-rn1-fs20000-A-thresh-strf.mat
% 
%     where rn* indicates which type of ripple noise was used.
% 
%     filename : a string. Must have 'rn1*', 'rn4*', 'rn8*', or 'rn16*' in
%     the filename.
% 
%     rn : Ripple noise number. 1, 4, 8, or 16


if ( ~isempty(findstr(filename, '-rn1-')) || ~isempty(findstr(filename, 'dmr')) )
    rn = 1;
elseif ( ~isempty(findstr(filename, '-rn4-')) )
    rn = 4;
elseif ( ~isempty(findstr(filename, '-rn8-')) )
    rn = 8;
elseif ( ~isempty(findstr(filename, '-rn16-')) )
    rn = 16;
else
    error('Wrong file. RN must be 1, 4, 8, or 16');
end

return;

