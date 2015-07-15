function ca_plot_cell_assemblies_exp_site_cadata(exp_site_cadata, startx)
% ca_plot_cell_assemblies_cadata Display cell assemblies, activities
% 
%     ca_plot_cell_assemblies_data(exp_site_cadata, startx)
%
%     exp_site_cadata : struct of data from cell assembly calculations for one
%     experiment and one site.
%
%     Has the form:
%
%     exp_site_cadata = 
% 
%            exp: experiment date
%           site: penetration number
%           stim: 'rn1', 'rn4', etc.
%             df: downsampling factor. Not overall downsampling factor
%         cadata: struct of cell assembly data
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


cadata = exp_site_cadata.cadata;
Activitymatrix = cadata.spktrain;
Patterns = cadata.Patterns;
Activities = cadata.Activities;
position = cadata.position;

dt = cadata.df / cadata.fsdvd * 1000; % spike train bin size, in ms


% Estimate pairwise channel correlations
correlationmat = corr(Activitymatrix');


figure;

% Plot the pairwise correlations
subplot(2,3,1);
imagesc(correlationmat);
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
zSpikeCount = zscore(Activitymatrix');
imagesc(zSpikeCount');
% imagesc(Activitymatrix);
xlim([startx startx+500]);
if ( isempty(position) )
    ylabel('Neuron #');
else
    tick = 1:length(position);
    set(gca,'ytick', unique_tick, 'yticklabel', unique_pos);
    %set(gca,'ytick', tick, 'yticklabel', position);
    ylabel('Position (um)');
end
xtick = get(gca,'xtick');
xticklabel = dt * xtick;
set(gca, 'xtick', xtick, 'xticklabel', xticklabel);
tickpref;



% Plot the cell assemblies / independent components
subplot(2,3,4);
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



% Plot the activities of the cell assemblies / time course of cell assembly activity
subplot(2,3,[5 6]);
plot(Activities');
xlim([startx startx+500]);
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


% Place a superior title over all the plots
exp = exp_site_cadata.exp;
site = exp_site_cadata.site;
stim = exp_site_cadata.stim;

set(0,'defaulttextinterpreter','none')
suptitle(sprintf('Exp %s site%.0f %s %.1f ms', ...
    exp, site, stim, ...
    cadata.df / cadata.fsdvd * 1000));



% for i = 1:size(Patterns,2)
%     figure;
%     stem(Patterns(:,i));
%     title(sprintf('Ind Comp #%.0f', i));
% end




return;



