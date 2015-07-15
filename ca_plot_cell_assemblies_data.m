function ca_plot_cell_assemblies_data(Activitymatrix, Patterns, Activities, position, fsdvd, totaldf, startx)
% ca_plot_cell_assemblies_data Plot cell assemblies, correlation, activations
% 
%     ca_plot_cell_assemblies_data(Activitymatrix, Patterns, Activities, position, fsdvd, totaldf, startx)
%
%     Activitymatrix : spike train matrix. m x n matrix. m = number of neurons, n = number of
%     time bins. Activitymatrix(i,j) = # of spikes by neuron i at time bin
%     j. 
%
%     Patterns : independent components from Activitymatrix. Each Pattern
%     is a cell assembly.
%
%     Activities : time course of cell assembly activity.
%
%     position : vector of electrode recording depths for each spike train
%     in the Activitymatrix. position(i) corresponds to the spike train for
%     Activitymatrix(i,:), the spike train for neuron i.
%
%     fsdvd : sampling rate of Audio-DVD player.
%
%     totaldf : total downsampling 
%
%     startx : Optional. beginning bin of raster and activity plots. Default = 0. 
% 

% Check input arguments
narginchk(3,7);



% Estimate pairwise channel correlations
correlationmat = corr(Activitymatrix');



figure;

subplot(2,3,1);
imagesc(correlationmat);
if ( nargin == 3 )
    xlabel('Neuron #');
    ylabel('Neuron #');
    startx = 0;
    position = [];
    dt = 1;
end

if ( nargin == 4 )
    if ( isempty(position) )
        xlabel('Neuron #');
        ylabel('Neuron #');
    else
        tick = 1:length(position);
        set(gca,'xtick', tick, 'xticklabel', position);
        set(gca,'ytick', tick, 'yticklabel', position);
        xlabel('Position (um)');
        ylabel('Position (um)');
    end
    startx = 0;
    dt = 1;
end

if ( nargin == 5 )
    if ( isempty(position) )
        xlabel('Neuron #');
        ylabel('Neuron #');
    else
        tick = 1:length(position);
        set(gca,'xtick', tick, 'xticklabel', position);
        set(gca,'ytick', tick, 'yticklabel', position);
        xlabel('Position (um)');
        ylabel('Position (um)');
    %     rotateticklabel(gca);
    end
    if ( isempty(fsdvd) )
        dt = 1;
    end
    startx = 0;
end

if ( nargin == 6 )
    if ( isempty(position) )
        xlabel('Neuron #');
        ylabel('Neuron #');
    else
        tick = 1:length(position);
        set(gca,'xtick', tick, 'xticklabel', position);
        set(gca,'ytick', tick, 'yticklabel', position);
        xlabel('Position (um)');
        ylabel('Position (um)');
    %     rotateticklabel(gca);
    end
    if ( isempty(fsdvd) || isempty(totaldf) )
        dt = 1;
    end
    startx = 0;
end

if ( nargin == 6 )
    if ( isempty(position) )
        xlabel('Neuron #');
        ylabel('Neuron #');
    else
        tick = 1:length(position);
        set(gca,'xtick', tick, 'xticklabel', position);
        set(gca,'ytick', tick, 'yticklabel', position);
        xlabel('Position (um)');
        ylabel('Position (um)');
    %     rotateticklabel(gca);
    end
    if ( ~isempty(fsdvd) && ~isempty(totaldf) )
        dt = totaldf / fsdvd * 1000; % spike train bin size, in ms
    else
        dt = 1;
    end
    if ( isempty(startx) )
        startx = 0;
    end
end



tickpref;



% subplot(2,3,[2 3]); 
% zSpikeCount = zscore(Activitymatrix');
% imagesc(zSpikeCount');
% % imagesc(Activitymatrix);
% xlim([startx startx+500]);
% if ( nargin == 3 )
%     ylabel('Neuron #');
% else
%     tick = 1:length(position);
%     set(gca,'ytick', tick, 'yticklabel', position);
%     ylabel('Position (um)');
% end
% tickpref;


% Plot the spike train matrix
subplot(2,3,[2 3]); 
zSpikeCount = zscore(Activitymatrix');
imagesc(zSpikeCount');
% imagesc(Activitymatrix);
xlim([startx startx+500]);
if ( isempty(position) )
    ylabel('Neuron #');
else
    tick = 1:length(position);
    set(gca,'ytick', tick, 'yticklabel', position);
    ylabel('Position (um)');
end
xtick = get(gca,'xtick');
xticklabel = dt * xtick;
set(gca, 'xtick', xtick, 'xticklabel', xticklabel);
tickpref;





subplot(2,3,4);
imagesc(Patterns);
xlabel('Assembly #');
if ( nargin == 3 )
    ylabel('Neuron #');
else
    tick = 1:length(position);
    set(gca,'ytick', tick, 'yticklabel', position);
    ylabel('Position (um)');
end
tickpref;
set(gca,'xtick', 1:size(Patterns,2), 'xticklabel', 1:size(Patterns,2));


% Plot the activities of the cell assemblies / time course of cell assembly activity
subplot(2,3,[5 6]);
plot(Activities');
xlim([startx startx+500]);
xlabel('Time (ms)');
tickpref;
[nr, nc] = size(Activities');
if ( nc > 0 )
    for i = 1:nc
        leg{i} = num2str(i);
    end
    legend(leg,'Location','Best');
end
xtick = get(gca,'xtick');
xticklabel = dt * xtick;
set(gca, 'xtick', xtick, 'xticklabel', xticklabel);
set(gcf,'position', [432 334 1121 620]);





% for i = 1:size(Patterns,2)
%     figure;
%     stem(Patterns(:,i));
%     title(sprintf('Ind Comp #%.0f', i));
% end




return;



