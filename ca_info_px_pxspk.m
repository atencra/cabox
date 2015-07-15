function iplugin = ca_info_px_pxspk(px, pxspk)
% ca_info_px_pxspk - information from probability distributions for filters
% 
% iplugin = ca_info_px_pxspk(px, pxspk)
% ----------------------------------
% 
% px : probability of a projection without respect to a spike
% 
% pxspk : probability of a projection given a spike
% 
% iplugin : information value. Estimated as 
% 
%             sum( pxspk(index) .* log2( pxspk(index) ./ px(index) ) )
%
% caa 3/10/09

px = px(:);
pxspk = pxspk(:);

index = find( px>0 & pxspk>0 ); % dividing by zero is bad!

iplugin = sum( pxspk(index) .* log2( pxspk(index) ./ px(index) ) );

return;









