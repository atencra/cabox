function cell_assembly_analysis_outline
% 1. Load stimulus matrix
% Determine number of times
% Determine how they correspond to times
% 
% 2. Create spike train(s) at resolution indicated by matrix.
% Align spike train relative to trigger beginning.
% Make spike train bins according to bins in stimulus matrix
% 
% 3. Make stimulus observation matrix and response vector
% 
% 4. Estimate cell assemblies.
% 
% 5. Calculate cell assembly STA using cell assembly activation pattern and
% the stimulus observation matrix.
% 
% 6. For each cell assembly, extract stimuli that precede activation.
% 
% 7. Cluster stimuli from 6.


1. Run ca_batch_calc_cell_assembly_sta_stim_binsize
from inside the folder with the spk-strfcmb files



2. Run ca_batch_calc_cell_assembly_nonlinearity
from inside folder where the results from 1. are stored.


3. Run ca_batch_calc_cell_assembly_nonlinearity_fit
after you have fun 2.










