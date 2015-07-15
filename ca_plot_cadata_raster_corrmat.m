function ca_plot_cadata_raster_corrmat(cadata, startx)
% ca_plot_cell_assemblies_cadata Display cell assemblies, activities
% 
%     ca_plot_cadata(cadata, startx)
%
%
%     cadata : struct holding cell assembly results. 
%     Has the following fields:
%
%     cadata = 
% 
%           spktrain: spike train matrix
%              fsdvd: sampling rate of dvd audio system
%                 df: overall envelope downsampling factor
%           position: unit position, in um
%           Patterns: cell assemblies. Each column is an assembly
%         Activities: time course of assemblies. Each row is an assembly
%                 nf: # freq's in STA
%              nlags: # time bins in STA
%             stamat: STA of each neuron. Each row is a neuron.
%          ca_stamat: STA of each assembly. Each row is an assembly
%
%     startx : Optional. Beginning bin of raster and activity plots. Default = 0. 
% 

% Check input arguments
narginchk(1,2);

if ( nargin == 1 )
    startx = 0;
end

if ( nargin == 2 )
    if ( isempty(startx) )
        startx = 0;
    end
end


Activitymatrix = cadata.spktrain;
Patterns = cadata.Patterns;
Activities = cadata.Activities;
position = cadata.position;

pos = zeros(size(position));
for i = 1:length(position)
    pos(i) = str2num(position{i});
end

index_position = find(pos < 1300);
Activitymatrix = Activitymatrix(index_position,:);
position = position(index_position);

index_raster = round(linspace(4,size(Activitymatrix,1),25));
Activitymatrix = Activitymatrix( index_raster, : );
position_vec = position(index_raster);



dt = cadata.df / cadata.fsdvd * 1000; % spike train bin size, in ms


% Estimate pairwise channel correlations
correlationmat = corr(Activitymatrix');

% Set diagonal to zero
correlationmat = diag2zero(correlationmat);


figure;

% Plot the pairwise correlations
subplot(1,3,3);
imagesc(correlationmat);
if ( isempty(position) )
    xlabel('Neuron #');
    ylabel('Neuron #');
else
    tick = 1:length(position_vec);
    pos = zeros(size(tick));
    for i = 1:length(position_vec)
        pos(i) = str2num(position_vec{i});
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
mn = totalmin(correlationmat);
mx = totalmax(correlationmat);
cmap = brewmaps('reds', 19);
cmap = flipud(cmap);
cmap = [1 1 1; cmap];
colormap(cmap);
set(gca, 'clim', [mn mx]);
colorbar;
xticklabel_rotate([],-90,[]);
tickpref;


% Plot the spike train matrix
subplot(1,3,[1 2]); 
zSpikeCount = zscore(Activitymatrix')';
zSpikeCount(zSpikeCount <= 0) = 0;
% zSpikeCount(zSpikeCount > 0) = 1;
imagesc(zSpikeCount);

% imagesc(Activitymatrix);
xlim([startx startx+500]);
if ( isempty(position_vec) )
    ylabel('Neuron #');
else
    set(gca,'ytick', unique_tick, 'yticklabel', unique_pos);
    %set(gca,'ytick', tick, 'yticklabel', position);
    ylabel('Position (um)');
end
xtick = get(gca,'xtick');
xticklabel = dt * xtick;
set(gca, 'xtick', xtick, 'xticklabel', xticklabel);
tickpref;
mn = totalmin(zSpikeCount);
mx = 0.15*totalmax(zSpikeCount);
cmap = brewmaps('reds', 19);
cmap = flipud(cmap);
cmap = [1 1 1; cmap];
colormap(cmap);
set(gca, 'clim', [mn mx]);



set(gcf,'position', [497 486 1230 345]);


set(0,'defaulttextinterpreter','none')
suptitle(sprintf('Bin size = %.1f ms', cadata.df / cadata.fsdvd * 1000));


return;



