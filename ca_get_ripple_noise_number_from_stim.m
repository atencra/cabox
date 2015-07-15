function rn = ca_get_ripple_noise_number_from_stim(stim)
% ca_get_ripple_noise_number_from_stim Get ripple noise from stim string
% 
%     rn = ca_get_ripple_noise_number_from_spk_filename(stim)
% 
%     stim : a string. Must have 'rn1*', 'rn4*', 'rn8*', or 'rn16*'.
% 
%     rn : Ripple noise number. 1, 4, 8, or 16


if ( strcmp(stim, 'rn1')  )
    rn = 1;
elseif ( strcmp(stim, 'rn4')  )
    rn = 4;
elseif ( strcmp(stim, 'rn8')  )
    rn = 8;
elseif ( strcmp(stim, 'rn16')  )
    rn = 16;
else
    error('Wrong file. RN must be 1, 4, 8, or 16');
end

return;

