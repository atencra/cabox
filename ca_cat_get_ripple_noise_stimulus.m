function [stimstr] = ca_cat_get_ripple_noise_stimulus(abspath, basefile, dft, dff)
% ca_cat_get_ripple_noise_stimulus - Cat experiment Ripple Noise stimulus matrix and parameters
%
% For cat experiments using the Neuralynx system, loads stimulus matrix 
% and parameters from *-matrix.mat and *_param.mat files.
%
% [stimstr] = ca_get_ripple_noise_stimulus(abspath, basefile, dft, dff)
% ===============================================================
%
% abspath : absolute path to the folder holding the stimulus files.
%
% basefile : initial file name for ripple stimulus. 
% Example : basefile = 'rn1-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-15min';
% basefile is obtained from: [basefile] = ca_get_base_downsampled_spr_filename(exp, stim)
%
%
% dft : temporal downsampling factor of the .spr file. Note that this is
% not the total downsampling factor. The original ripple was created at 
% 96kHz. The envelope file was not saved at this resolution. It was 
% downsampled by a factor of 48. The dft factor describes how the 
% downsampled envelope file was further downsampled. To get the total 
% downsampling factor of the stimulus, you multiply 48 times the dft factor: 
% 48 * dft.
% 
% The default value of dft is 10, or 5 ms bin size for 96 kHz spr files.
% (48*dft = 480; 480/96kHz * 1000 = 5 ms)
%
% dff : spectral downsampling factor. Default = 5. 
%
% dft and dff are needed because they included as part of the stimulus
% envelope filename. They help identify the correct file name.
%
% stimstr : struct holding the stimulus matrix, the sampling rate of the
% DVD-Audio player, the overall temporal downsampling factor, the 
% temporal axis of the envelope file, and the frequency axis of the 
% envelope file.
%
% stimstr also holds the parameter and matrix file names that were used
% to obtain the stimulus data.
%
% caa 1/27/15


narginchk(0,4);

if ( nargin == 0 )
    abspath1 = 'C:\Users\craig\Data\Stimuli';
    abspath2 = 'C:\Users\craig\Stimuli';

    if ( exist(abspath1, 'dir') )
        abspath = abspath1;
    elseif ( exist(abspath2, 'dir') )
        abspath = abspath2;
    else
        error('Wrong absolute path to ripple noise files.');
    end
    
    basefile = 'rn1-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-15min';
    dft = 10;
    dff = 5;
end

if ( nargin > 0 && isempty(abspath) )
    abspath1 = 'C:\Users\craig\Data\Stimuli';
    abspath2 = 'C:\Users\craig\Stimuli';

    if ( exist(abspath1, 'dir') )
        abspath = abspath1;
    elseif ( exist(abspath2, 'dir') )
        abspath = abspath2;
    else
        error('Wrong absolute path to ripple noise files.');
    end
end




if ( nargin == 1 )
    basefile = 'rn1-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-15min';
    dft = 10;
    dff = 5;
end


if ( nargin == 2 )
    dft = 10;
    dff = 5;
end

if ( nargin > 2 && isempty(dft) )
    dft = 10;
end


if ( nargin == 3 )
    dff = 5;
end


if ( nargin == 4 && isempty(dff) )
    dff = 5;
end



file_middle = sprintf('%s_DFt%.0f_DFf%.0f', basefile, dft, dff);

paramfile = fullfile(abspath, sprintf('%s_param.mat', file_middle));
sparams = load(paramfile, 'Fs', 'DF', 'taxis', 'faxis');

matrixfile = fullfile(abspath, sprintf('%s-matrix.mat', file_middle));
sstim = load(matrixfile, 'stimulus');

stimstr.matrixfile = matrixfile;
stimstr.paramfile = paramfile;
stimstr.Fs = sparams.Fs;
stimstr.DF = sparams.DF;
stimstr.taxis = sparams.taxis;
stimstr.faxis = sparams.faxis;
stimstr.stimulus = sstim.stimulus;
stimstr.stimfile = matrixfile;
stimstr.paramfile = paramfile;

clear('stimulus', 'taxis', 'faxis', 'Fs', 'DF');

return






