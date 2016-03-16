function ca_plot_cadata_corrmat_ensembles(cadata, startx)
% ca_plot_cadata_corrmat_ensembles Plot cell assembly analysis from processed data
% 
%     ca_plot_cadata_corrmat_ensembles(cadata)
%     -----------------------------------------------------------------
%     Plots the correlation matrix, neural ensembles, part of the spike
%     train matrix, part of the predicted activities, and the
%     eigenvalue spectrum.
%
%     cadata : struct holding the spike train matrix, the position of each
%     neuron, the cell assembly patterns, and the cell assembly activities.
%
%     cadata already contains the cell assembly analysis. This function
%     recomputes the analysis, and additionally plots the eigenvalue
%     spectrum, which is not included in cadata.
%

narginchk(1,2);

if ( nargin == 1 )
    startx = 0;
end

if ( nargin == 2 )
    if ( isempty(startx) )
        startx = 0;
    end
end


spkmat = cadata.spktrain;
position = cadata.position;
dt = cadata.df / cadata.fsdvd * 1000; % spike train bin size, in ms


spkmat = cadata.spktrain;
spkmat = spkmat(1:25,:);
nedata = ca_spkmatrix_to_ensembles(spkmat);
ensembles = nedata.ensembles;
evals = nedata.eigenvalues;
lambda_max = nedata.lambda_max;
lambda_min = nedata.lambda_min;

Activities = assembly_activity(ensembles, spkmat);

corrmat = nedata.corrmat;
corrmat = diag2zero(corrmat);



figure;

% Plot the pairwise correlations
subplot(2,3,1);
imagesc(corrmat);
if ( isempty(position) )
    xlabel('Neuron #');
    ylabel('Neuron #');
else
    tick = 1:length(position);
    pos = zeros(size(tick));
    for i = 1:length(position)
        pos(i) = str2num(position{i});
    end
    
    unique_pos = unique(pos);
    unique_tick = zeros(size(unique_pos));
    for i = 1:length(unique_pos)
        unique_tick(i) = find(pos == unique_pos(i), 1);
    end
    
    set(gca,'xtick', unique_tick, 'xticklabel', unique_pos);
    set(gca,'ytick', unique_tick, 'yticklabel', unique_pos);
%     set(gca,'xtick', tick, 'xticklabel', pos);
%     set(gca,'ytick', tick, 'yticklabel', pos);
    xlabel('Position (um)');
    ylabel('Position (um)');
end

tickpref;


% Plot the spike train matrix
subplot(2,3,[2 3]); 
zSpikeCount = zscore(spkmat')';
index_raster = round(linspace(1,size(zSpikeCount,1),25));
zSpikeCount = zSpikeCount(index_raster,:);
position_vec = position(index_raster);
imagesc(zSpikeCount);

% imagesc(Activitymatrix);
xlim([startx startx+500]);
if ( isempty(position_vec) )
    ylabel('Neuron #');
else
    tick = 1:length(position_vec);
    set(gca,'ytick', unique_tick, 'yticklabel', unique_pos);
    %set(gca,'ytick', tick, 'yticklabel', position);
    ylabel('Position (um)');
end
xtick = get(gca,'xtick');
xticklabel = dt * xtick;
set(gca, 'xtick', xtick, 'xticklabel', xticklabel);
tickpref;
mn = totalmin(zSpikeCount);
mx = 0.75*totalmax(zSpikeCount);
cmap = brewmaps('reds', 19);
cmap = flipud(cmap);
cmap = [1 1 1; cmap];
colormap(cmap);



% Plot the cell assemblies / independent components
subplot(2,3,4);
temp = abs(ensembles);
temp(temp<0.2) = 0;
%temp(temp>=0.2) = 1;
imagesc(temp);
xlabel('Assembly #');
if ( nargin == 3 )
    ylabel('Neuron #');
else
    tick = 1:length(position);
    set(gca,'ytick', unique_tick, 'yticklabel', unique_pos);
    %set(gca,'ytick', tick, 'yticklabel', position);
    ylabel('Position (um)');
end
tickpref;
set(gca,'xtick', 1:size(nedata.ensembles,2), 'xticklabel', 1:size(nedata.ensembles,2));




% Plot the activities of the cell assemblies / time course of cell assembly activity
subplot(2,3,[5 6]);
plot(Activities');
xlim([startx startx+500]);
% minmin = min( min(Activities,[], 1) )
% maxmax = max( max(Activities,[], 1) )
% ylim([minmin maxmax]);

xlabel('Time (ms)');
tickpref;
[~, nc] = size(Activities');
if ( nc > 0 )
    for i = 1:nc
        leg{i} = num2str(i);
    end
%     legend(leg,'Location', 'Best');
    hl = legend(leg);
    set(hl,'position',[0.91 0.15 0.1 nc*(0.15/6)]); % this position seems to work
end
xtick = get(gca,'xtick');
xticklabel = dt * xtick;
set(gca, 'xtick', xtick, 'xticklabel', xticklabel);
set(gcf,'position', [100 100 1121 620]);

set(0,'defaulttextinterpreter','none')
suptitle(sprintf('Bin size = %.1f ms', cadata.df / cadata.fsdvd * 1000));



figure;
hold on;
evalsort = sort(evals);
evalnum = 1:length(evalsort);

index = find(evalsort > lambda_max);
plot(evalnum(index), evalsort(index), 'ko', ...
'markerfacecolor', 'k');

index = find(evalsort <= lambda_max & evalsort >= lambda_min);
plot(evalnum(index), evalsort(index), 'o', ...
'markerfacecolor', [0.6 0.6 0.6], ...
'markeredgecolor', [0.6 0.6 0.6]);

index = find(evalsort < lambda_min);
plot(evalnum(index), evalsort(index), 'ko', ...
'markerfacecolor', 'k');


plot(1:length(evals), lambda_max * ones(size(evals)), 'r-');
title('E-val dist and M-P boundary');

plot(evalnum, lambda_min * ones(size(evalnum)), 'r-');

xlabel('Eigenvalue #');
ylabel('Eigenvalue');

tickpref;



return;

















