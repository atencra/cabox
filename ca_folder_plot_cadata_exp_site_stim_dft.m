function ca_folder_plot_cadata_exp_site_stim_dft(exp, site, stim, dft)
% ca_plot_cell_assembly_stim_df_sta Plot stim, bin size cell assembly results
% 
%     ca_plot_cell_assembly_stim_df_sta(exp, site, stim, dft)
%     --------------------------------------------------------------
%     exp : experiment, a string. For neuralynx A/D recordings, it will have
%     the form '2013-11-25'. For Intan A/D recordings, it will be similar
%     to '141214_222459'
%
%     site : penetration number. An integer.
%
%     stim : cell array of strings. Either 'rn1', 'rn4', 'rn8', or 'rn16'.
%     Default = 'rn1'. To process multiple stimuli, set stim to a cell
%     array of strings. Example: stim = {'rn1', 'rn1y'}.
%
%     dft : temporal downsampling factor for stimulus matrix. An integer or
%     vector of integers.
%
%     This is not the total downsampling factor, but is the downsampling 
%     factor used after the original stimulus was created. 
%
%     During the creation process, the sound envelope was downsampled 
%     by a factor of 48, giving a time resolution of 48/96000 = 0.5 ms. 
%     Then the envelope was further downsampled by a factor of dft, 
%     usually [5, 10, 20, 40, or 100]. Thus, the time resolution 
%     was [2.5, 5, 10, 20, or 50] ms.
%

library('santosbox');
library('UtilitiesColormaps');
library('dataviz');


narginchk(2,4);

if ( nargin == 2 )
    stim = {'rn1'};
    dft = 10;
end

if ( nargin == 3 )
    dft = 10;
end


if ( ~iscell(exp) && ~isempty(exp) ) % it's a string
    sexp{1} = exp; 
    exp = sexp;
end


if ( ~iscell(site) && ~isempty(site) ) % it's a vector
    for i = 1:length(site)
        s{i} = num2str(site);
    end
    site = s;
end


if ( ~iscell(stim) && ~isempty(stim) )
    stimulus_type = stim;
    clear('stim');
    stim{1} = stimulus_type;
end


if ( isempty(stim) )
    disp('here')
    stim = {'rn1'};
end

if ( isempty(dft) )
    dft = 10;
end



dtotal = [];
if ( isempty(exp) )
    if ( isempty(site) )
        for i = 1:length(stim)
            for  ii = 1:length(dft)
                file = sprintf('*-site*-*um-*db-%s-fs20000-A-spk-strfcmb-ca-%.0fdft.mat', stim{i}, dft(ii) );
                d = dir( file );
                dtotal = [dtotal d(:)'];
            end
        end
    else
        for i = 1:length(site)
            for ii = 1:length(stim)
                for  iii = 1:length(dft)
                    file = sprintf('*-site%.0f-*um-*db-%s-fs*-*-ca-%.0fdft.mat', ...
                    site{i}, stim{ii}, dft(iii) ); 
                    d = dir( file ); 
                    dtotal = [dtotal d(:)'];
                end
            end
        end
    end

else
    if ( isempty(site) )
        for i = 1:length(exp)
            for ii = 1:length(stim)
                for  iii = 1:length(dft)
                    file = sprintf('%s-site*-*um-*db-%s-fs*-*-ca-%.0fdft.mat', ...
                    exp{i}, stim{ii}, dft(iii) ); 
                    d = dir( file ); 
                    dtotal = [dtotal d(:)'];
                end
            end
        end
    else
        for i = 1:length(exp)
            for ii = 1:length(site)
                for iii = 1:length(stim)
                    for  iv = 1:length(dft)
                        file = sprintf('%s-site*-*um-*db-%s-fs*-*-ca-%.0fdft.mat', ...
                            exp{i}, site{ii}, stim{iii}, dft(iv) ); 
                        d = dir( file ); 
                        dtotal = [dtotal d(:)'];
                    end % (for iv)
                end % (for iii)
            end % (for ii)
        end % (for i)
    end % (if-else)
end % (if-else)



for i = 1:length(dtotal)
    files{i} = dtotal(i).name;
end
cafiles = unique(files);


exp_site_cadata = [];
for i = 1:length(cafiles)
    fprintf('%s\n', cafiles{i});
    load(files{i},'cadata');
    [exp, site, stim, dft] = ca_get_exp_site_stim_dft_from_file(cafiles{i});
    exp_site_cadata(i).exp = exp;
    exp_site_cadata(i).site = site;
    exp_site_cadata(i).stim = stim;
    exp_site_cadata(i).dft = dft;
    exp_site_cadata(i).cadata = cadata;
end % (for i)                
            
if ( isempty(exp_site_cadata) )
    error('No data loaded');
end

% for i = 1:length(exp_site_cadata)
%     exp_site_cadata(i)
% end



% Plot the data you have loaded
% for i = 1:length(exp_site_cadata)
%     ca_plot_cell_assemblies_exp_site_cadata( exp_site_cadata(i) );
% end % (for i)


% [rn# df #assemblies]
data = zeros(length(exp_site_cadata),3);
for i = 1:length(exp_site_cadata)
    data(i,1) = exp_site_cadata(i).stim;
    data(i,2) = exp_site_cadata(i).cadata.df;
    data(i,3) = size(exp_site_cadata(i).cadata.Patterns,2);
end

FsDVD = exp_site_cadata(1).cadata.fsdvd;

data(:,2) = data(:,2) ./ FsDVD * 1000;


figure;
subplot(1,2,1);
hold on;
unique_bin = unique(data(:,2));
unique_stim = unique(data(:,1));
cmap = brewmaps('blues',length(unique_stim)+1);
for i = 1:length(unique_stim)
    index = find(data(:,1) == unique_stim(i));
    hp = plot(data(index,2), jitter(data(index,3),1), ...
        'ko-', 'markerfacecolor', cmap(i,:) );
    set(hp, 'color', cmap(i,:) );
end % (for i)
for i = 1:length(unique_stim)
    leg{i} = sprintf('RN%.0f', unique_stim(i));
end % (for i)
legend(leg);
set(gca,'xscale', 'log');
set(gca,'xtick', unique_bin, 'xticklabel', unique_bin);
tickpref;
xlabel('Bin Size (ms)');
ylabel('# Assemblies');
xlim([0.75*min(unique_bin) 1.1*max(unique_bin)]);
limit = ylim;
ylim([0 limit(2)]);



subplot(1,2,2);
hold on;
unique_bin = unique(data(:,2));
cmap = brewmaps('blues',length(unique_bin)+1);
for i = 1:length(unique_bin)
    index = find(data(:,2) == unique_bin(i));
    hp = plot(data(index,1), jitter(data(index,3),1), ...
        'ko-', 'markerfacecolor', cmap(i,:) );
    set(hp, 'color', cmap(i,:) );
end % (for i)
for i = 1:length(unique_bin)
    leg{i} = sprintf('%.1f ms', unique_bin(i));
end % (for i)
legend(leg);
set(gca,'xscale', 'log');
set(gca,'xtick', unique_stim, 'xticklabel', unique_stim);
tickpref;
xlabel('RN #');
ylabel('# Assemblies');
limit = xlim;
xlim([0.75*min(unique_stim) 1.1*max(unique_stim)]);
limit = ylim;
ylim([0 limit(2)]);

set(gcf,'position', [65 342 983 420]);
print_mfilename(mfilename);


return;





function [exp, site, stim, dft] = ca_get_exp_site_stim_dft_from_file(file)

index = findstr(file, '-');
index_end = findstr(file, 'dft.mat');

exp = file(1:index(1)-1);
site = str2double( file(index(1)+5:index(2)-1) );
stim = str2double( file(index(4)+3:index(5)-1) );
dft = str2double( file(index(end)+1:index_end-1) );

return;


















