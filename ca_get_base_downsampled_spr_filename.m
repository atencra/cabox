function [basefile] = ca_get_base_downsampled_spr_filename(exp, stim)
% get_base_downsampled_spr_filename - base filename of envelope stimulus for each given experiment.
%
% [basefile] = ca_get_base_downsampled_spr_filename(exp, stim)
% ===============================================================
%
% exp : experiment name.  A string. Example: '2002-7-30'
%
% The date of the experiment is used to get the correct
% stimulus file.
%
% stim : type of stimulus. Either 'dmr' or 'rn'. Used to get the 
% correct stimulus file.
%
% specfile : filename for downsampled ripple envelope file.
%

if ( nargin == 1 )
    stim = 'rn1';
end


switch exp

case {'2002-6-19'}
   basefile = 'dmr-500flo-20000fhi-4SM-40TM-40db-44khz-10DF-10min';

case {'2002-7-30', '2002-8-26', '2002-8-27', '2003-3-5', '2003-4-8'}
   if ( ~isempty(findstr(stim, 'dmr') ) )
      basefile = 'dmr-500flo-20000fhi-4SM-40TM-40db-44khz-10DF-15min';
   elseif ( ~isempty(findstr(stim, 'rn') ) )
      basefile = 'rn-500flo-20000fhi-4SM-40TM-40db-44khz-10DF-15min';
   end

case {'2003-5-6', '2003-5-27'}
   basefile = 'rn1-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-15min';

case {'2003-10-21', '2003-10-28', '2003-11-12', '2004-1-14'}
   basefile = 'dmr-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-21min';

case {'2013-9-18', '2013-09-18', '2013-9-19', '2013-09-19', '2013-10-17', '2013-11-25', '2013-12-19', ...
'2014-3-6', '2014-3-25'}
   
   if ( ~isempty(findstr(stim, 'rn1') ) )
      basefile = 'rn1-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-10min';
   elseif ( ~isempty(findstr(stim, 'rn4') ) )
      basefile = 'rn4-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-10min';
   elseif ( ~isempty(findstr(stim, 'rn8') ) )
      basefile = 'rn8-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-10min';
   elseif ( ~isempty(findstr(stim, 'rn8') ) )
      basefile = 'rn16-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-10min';
   end

case {'2003-11-19'}
   basefile = 'dmr-50flo-40000fhi-4SM-500TM-40db-48DF-21min';

case {'2003-11-24'}
   basefile = 'dmr-50flo-40000fhi-4SM-500TM-40db-48DF-21min';

case {'2004-2-17'}
   basefile = 'dmr-50flo-20000fhi-4SM-500TM-40db-44khz-22DF-11min';

otherwise
   basefile = [];
end


return






