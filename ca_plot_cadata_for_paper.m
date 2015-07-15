function ca_plot_cadata_for_paper(cadata, startx, nbins)
% ca_plot_cell_assemblies_cadata Display cell assemblies, activities
% 
%     ca_plot_cadata(cadata, startx, endx)
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
%     nbins : Optional. Number of bins to show in the raster. Default = 250. 
% 

% Check input arguments
narginchk(1,3);

if ( nargin == 1 )
    startx = 1;
    endx = 250;
end

if ( nargin == 2 )
    if ( isempty(startx) )
        startx = 1;
    end
    endx = startx + 250 - 1;
end

if ( nargin == 3 )
    if ( isempty(startx) )
        startx = 1;
    end
    if ( isempty(nbins) )
        endx = startx + 250 - 1;
    end
    endx = startx + nbins - 1;
end
if (startx == 0)
    startx = 1;
end



Activitymatrix = cadata.spktrain;
% index = setdiff(1:size(Activitymatrix,1), [20 27 32 33 37]);
% Activitymatrix = Activitymatrix(index,:);
Patterns = cadata.Patterns;
Activities = cadata.Activities;
position = cadata.position;
% position = position(index);

pos = zeros(size(position));
for i = 1:length(position)
    pos(i) = str2num(position{i});
end

pos


zSpikeCount = zscore(Activitymatrix')';
zSpikeCount(zSpikeCount <= 0) = 0;

pop_rate = sum(zSpikeCount,1);
spk_vec = sum(Activitymatrix,1);

index = startx:endx;
zSpikeCount = zSpikeCount(:,index);
raster = Activitymatrix(:,index);


% index_raster = round(linspace(4,size(Activitymatrix,1),50));
% Activitymatrix = Activitymatrix( index_raster, : );
% position_vec = position(index_raster);

position_vec = position;



dt = cadata.df / cadata.fsdvd * 1000; % spike train bin size, in ms

% Estimate pairwise channel correlations
correlationmat = corr(Activitymatrix');

% Set diagonal to zero
correlationmat = diag2zero(correlationmat);





figure;

position

% Plot the pairwise correlations
% subplot(1,3,3);
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
mx = 0.65*totalmax(correlationmat);
cmap = brewmaps('reds', 19);
cmap = flipud(cmap);
cmap = [1 1 1; cmap];
colormap(cmap);
set(gca, 'clim', [mn mx]);
% colorbar;
xticklabel_rotate([],-90,[]);
tickpref;

close all;



figure;

subplot(5,1,[1 3]);
imagesc(zSpikeCount);
ytick = 1:2:size(zSpikeCount,1);
yticklabel = position(ytick);
ytick
set(gca,'ytick', ytick, 'yticklabel', yticklabel);

xtick = get(gca,'xtick');
xticklabel = dt * (xtick+startx);
set(gca, 'xtick', xtick, 'xticklabel', xticklabel);
tickpref;
mn = totalmin(zSpikeCount);
mx = 0.2*totalmax(zSpikeCount);
cmap = brewmaps('reds', 19);
cmap = flipud(cmap);
cmap = [1 1 1; cmap];
colormap(cmap);
set(gca, 'clim', [mn mx]);

set(0,'defaulttextinterpreter','none')
suptitle(sprintf('Bin size = %.1f ms', dt));



subplot(5,1,4);
time = index * dt;
plot(time, spk_vec(index), 'k-');
ylim([0 max(spk_vec(index))]);
tickpref;
ylabel('#Spikes');
box off;

subplot(5,1,5);
time = index * dt;
plot(time, pop_rate(index), 'k-');
ylim([min(pop_rate(index)) max(pop_rate(index))]);
tickpref;
xlabel('Time (s)');
ylabel('Summed z-scores');
box off;

set(gcf,'position', [680 336 439 591]);
% colormap(brewmaps('rdbu',21))


return;










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




















figure;

subplot(4,1,[1 3]);
% Plot the raster as an image matrix
imagesc(raster);

% % Put labels on the rows to denote the neuron number or neuron depth
% if ( isempty(position) )
%     ylabel('Neuron #');
% else
%     tick = 1:length(position);
%     set(gca,'ytick', unique_tick, 'yticklabel', unique_pos);
%     %set(gca,'ytick', tick, 'yticklabel', position);
%     ylabel('Position (um)');
% end

ytick = 1:2:size(zSpikeCount,1);
set(gca,'ytick', ytick, 'yticklabel', ytick);


% Set x-axis labels, ticks
xtick = get(gca,'xtick');
xticklabel = dt * (xtick+startx);

% xticklabel = dt * xtick;
set(gca, 'xtick', xtick, 'xticklabel', xticklabel);
tickpref;
set(0,'defaulttextinterpreter','none')
suptitle(sprintf('Bin size = %.1f ms', cadata.df / cadata.fsdvd * 1000));


subplot(4,1,4);
time = index * dt;
plot(time, spk_vec(index), 'k-');
ylim([0 max(spk_vec(index))]);
tickpref;
ylabel('#Spikes');
box off;


set(gcf,'position', [680 336 439 591]);
colormap(brewmaps('rdbu',21))



return;



figure;
% Plot the cell assemblies / independent components
% subplot(2,3,4);
imagesc(Patterns);
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
set(gca,'xtick', 1:size(Patterns,2), 'xticklabel', 1:size(Patterns,2));
set(0,'defaulttextinterpreter','none')
suptitle(sprintf('Bin size = %.1f ms', cadata.df / cadata.fsdvd * 1000));



figure;
% Plot the activities of the cell assemblies / time course of cell assembly activity
% subplot(2,3,[5 6]);
plot(Activities');
xlim([startx endx]);
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


% % Place a superior title over all the plots
% exp = exp_site_cadata.exp;
% site = exp_site_cadata.site;
% stim = exp_site_cadata.stim;

set(0,'defaulttextinterpreter','none')
suptitle(sprintf('Bin size = %.1f ms', cadata.df / cadata.fsdvd * 1000));


% for i = 1:size(Patterns,2)
%     figure;
%     stem(Patterns(:,i));
%     title(sprintf('Ind Comp #%.0f', i));
% end


return;



