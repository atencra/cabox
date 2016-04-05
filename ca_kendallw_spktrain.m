function [w, chi2, pval] = ca_kendallw_spktrain(spkmat)

[nrows, ncols] = size(spkmat);

% Put variables in each column, observations on each row
if ( ncols > nrows )
    spkmat = spkmat';
end

[w, chi2, pval] = kendallw(spkmat);

return;