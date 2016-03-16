function ca_plot_multi_cadata_corrmat_ensembles(cadatacell, startx)
% ca_plot_multi_cadata_corrmat_ensembles Plot multiple ensembles from processed data
% 
%     ca_plot_multi_cadata_corrmat_ensembles(cadata)
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


mx = 0;
for i = 1:length(cadatacell)
    cadata = cadatacell{i};
    [nr, nc] = size(cadata.spktrain);
    [nr nc]

    if ( mx == 0 )
        mx = nr;
    else
        mx = min([mx nr]);
    end
end % (for i)

mx = min([mx 25]);


nedatacell = cell(size(cadatacell));
NumberEnsembles = zeros(1,length(cadatacell));

for i = 1:length(cadatacell)

    cadata = cadatacell{i};

    spkmat = cadata.spktrain;
    position = cadata.position;
    
%     pos = zeros(size(position));
%     for j = 1:length(position)
%         pos(j) = str2double(position{j});
%     end % (for j)

    dt = cadata.df / cadata.fsdvd * 1000; % spike train bin size, in ms

    spkmat = cadata.spktrain;
    index_raster = round(linspace(1,size(spkmat,1),25));

%     index = find(pos<1100)
%     index_raster = index;
    

    spkmat = spkmat(index_raster,:);
    nedata = ca_spkmatrix_to_ensembles(spkmat);
    nedata.position = position(index_raster);
    nedatacell{i} = nedata;
    NumberEnsembles(i) = size(nedata.ensembles,2);
end % (for i)


figure;

for i = 1:length(cadatacell)

    nedata = nedatacell{i};

    position = nedata.position;

    ensembles = nedata.ensembles;
    evals = nedata.eigenvalues;
    lambda_max = nedata.lambda_max;
    lambda_min = nedata.lambda_min;

    Activities = assembly_activity(ensembles, spkmat);

    corrmat = nedata.corrmat;
    corrmat = diag2zero(corrmat);



    % Plot the pairwise correlations
    subplot(length(cadatacell),3, (i-1)*3+1);
    corrmat(corrmat>0.25) = 0.25;
    imagesc(corrmat);
    if ( isempty(position) )
        xlabel('Neuron #');
        ylabel('Neuron #');
    else
        tick = 1:length(position);
        pos = zeros(size(tick));
        for n = 1:length(position)
            pos(n) = str2num(position{n});
        end

        unique_pos = unique(pos);
        unique_tick = zeros(size(unique_pos));
        for n = 1:length(unique_pos)
            unique_tick(n) = find(pos == unique_pos(n), 1);
        end

        set(gca,'xtick', unique_tick, 'xticklabel', unique_pos);
        set(gca,'ytick', unique_tick, 'yticklabel', unique_pos);
    %     set(gca,'xtick', tick, 'xticklabel', pos);
    %     set(gca,'ytick', tick, 'yticklabel', pos);
        xlabel('Position (um)');
        ylabel('Position (um)');
    end

    tickpref;

%     mn = totalmin(zSpikeCount);
%     mx = 0.75*totalmax(zSpikeCount);
    cmap = brewmaps('reds', 19);
    cmap = flipud(cmap);
    cmap = [1 1 1; cmap];
    colormap(cmap);




    subplot(length(cadatacell),3,(i-1)*3+2);

    hold on;
    evalsort = sort(evals);
    evalnum = 1:length(evalsort);


    plot(evalnum, evalsort, 'ko', ...
    'markerfacecolor', 'k');


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

    %plot(evalnum, lambda_min * ones(size(evalnum)), 'r-');

    xlabel('Eigenvalue #');
    ylabel('Eigenvalue');

    tickpref;




    % Plot the cell assemblies / independent components
    subplot(length(cadatacell), 3, (i-1)*3+3);

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
    xlim([0 max(NumberEnsembles)+1]);

end


set(gcf,'position', [627 198 779 793]);


return;

















