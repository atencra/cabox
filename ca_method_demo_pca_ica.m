function ca_method_demo_pca_ica
% ca_method_demo_pca_ica Illustrates PCA/ICA neural ensemble analysis


close all;
clear all % just clearing workspace

% define mean firing rate
Network_opts.meanspikebin = 0.1;

% define number of neurons
Network_opts.nneurons = 10;

% define number of bins
Network_opts.nbins = 100000;

% above define assembly membership

% line below sets neurons 1,2,3 and 4 to be in assembly 1
Assembly_opts.assembly_neurons{1} = [1 2 3 4]; 
% line below sets neurons 5,6 and 7 to be in assembly 2
Assembly_opts.assembly_neurons{2} = [5 6 7]; 

% defines number of activation bins
Assembly_opts.number_of_activations = 1000;

% defines mean rate in activation bins
Assembly_opts.meanspikerate_activations = 5000;

% running the function
Activitymatrix = toy_simulation(Network_opts,Assembly_opts);


% correlationmat = corr(Activitymatrix');
% figure;
% clf;
% imagesc(correlationmat)


clear all % just clearing workspace




% Create simulated spike train matrix with two neurons shared

% mean firing rate
Network_opts.meanspikebin = 1;

% Number of neurons
Network_opts.nneurons = 10;

% Number of bins
Network_opts.nbins = 10000;


% Ensemble membership - no shared neurons
Assembly_opts.assembly_neurons{1} = [1 2 3 4]; 
Assembly_opts.assembly_neurons{2} = [5 6 7 8]; 


% Number of activation bins
Assembly_opts.number_of_activations = 300;

% Rate in activated bins
Assembly_opts.meanspikerate_activations = 3;


% Simulate, analyze, and plot
% simulated_spktrain_analysis(Network_opts, Assembly_opts);




% Create simulated spike train matrix with two neurons shared
Network_opts.meanspikebin = 0.1;
Network_opts.nneurons = 8;
Network_opts.nbins = 10000;

Assembly_opts.assembly_neurons{1} = [1 2 3 4 5]; 
Assembly_opts.assembly_neurons{2} = [4 5 6 7 8]; 

Assembly_opts.number_of_activations = 300;
Assembly_opts.meanspikerate_activations = 2.5;

simulated_spktrain_analysis(Network_opts, Assembly_opts);

return;





function simulated_spktrain_analysis(Network_opts, Assembly_opts)


Activitymatrix = toy_simulation(Network_opts,Assembly_opts);


corrmat = corr(Activitymatrix');




% Analyze using PCA
opts.threshold.method = 'MarcenkoPastur';
opts.AssemblyTemplate.method = 'PCA';
opts.Patterns.method = 'PCA';
[pcaAssemblyTemplates,pcadata] = assembly_patterns_caa(Activitymatrix,opts);
[pcaActivities] = assembly_activity(pcaAssemblyTemplates,Activitymatrix);



% Analyze using ICA
opts.threshold.method = 'MarcenkoPastur';
opts.AssemblyTemplate.method = 'ICA';
opts.Patterns.method = 'ICA';
opts.AssemblyTemplate.number_of_iterations = 500;
opts.Patterns.number_of_iterations = 500;
[icaAssemblyTemplates,icadata] = assembly_patterns_caa(Activitymatrix,opts);
[icaActivities] = assembly_activity(icaAssemblyTemplates,Activitymatrix);






% figure;
% clf;
% subplot(2,1,1);
% imagesc(corrmat)
% tickpref;
% xlabel('Neuron #');
% ylabel('Neuron #');
% title('Correlation Matrix');
% 
% subplot(2,1,2);
% corrmat(corrmat>0.5*max(max(corrmat))) = 0.65*max(max(corrmat));
% imagesc(corrmat)
% tickpref;
% xlabel('Neuron #');
% ylabel('Neuron #');
% 
% 
% 
% figure;
% clf;
% subplot(2,3,1);
% x = 1:Network_opts.nneurons;
% index = Assembly_opts.assembly_neurons{1};
% y = zeros(size(x));
% y(index) = 1;
% stem(x,y,'filled', 'k');
% box off;
% xtick = 1:Network_opts.nneurons;
% set(gca,'xtick',xtick,'xticklabel', xtick);
% ytick = 0:0.5:1;
% set(gca,'ytick',ytick,'yticklabel', ytick);
% xlim([0 Network_opts.nneurons+1]);
% ylim([-0.25 1.1]);
% tickpref;
% title(sprintf('Model\nEnsemble #1'));
% 
% 
% subplot(2,3,4);
% x = 1:Network_opts.nneurons;
% index = Assembly_opts.assembly_neurons{2};
% y = zeros(size(x));
% y(index) = 1;
% stem(x,y,'filled', 'k');
% box off;
% xtick = 1:Network_opts.nneurons;
% set(gca,'xtick',xtick,'xticklabel', xtick);
% ytick = 0:0.5:1;
% set(gca,'ytick',ytick,'yticklabel', ytick);
% xlim([0 Network_opts.nneurons+1]);
% ylim([-0.25 1.1]);
% tickpref;
% title('Ensemble #2');
% 
% 
% 
% subplot(2,3,2)
% x = pcaAssemblyTemplates(:,1);
% if ( abs(max(x)) < abs(min(x)) )  
%     x = -x;
% end
% stem(x,'filled', 'k');
% box off;
% xlim([0 Network_opts.nneurons+1]);
% range = max(x) - min(x);
% ylim([min(x)-0.05*range max(x)+0.05*range]);
% tickpref;
% xtick = 1:Network_opts.nneurons;
% set(gca,'xtick',xtick,'xticklabel', xtick);
% ylabel('PC Value');
% title(sprintf('PCA Analysis\nEnsemble #1'));
% 
% 
% 
% subplot(2,3,5);
% x = pcaAssemblyTemplates(:,2);
% if ( abs(max(x)) < abs(min(x)) )  
%     x = -x;
% end
% stem(x,'filled', 'k');
% box off;
% xlim([0 Network_opts.nneurons+1]);
% range = max(x) - min(x);
% ylim([min(x)-0.05*range max(x)+0.05*range]);
% tickpref;
% set(gca,'xtick',xtick,'xticklabel', xtick);
% xlabel('Neuron #');
% ylabel('PC Value');
% title(sprintf('Ensemble #2'));
% 
% 
% 
% subplot(2,3,3)
% x = icaAssemblyTemplates(:,1);
% if ( abs(max(x)) < abs(min(x)) )  
%     x = -x;
% end
% stem(x,'filled', 'k');
% box off;
% xlim([0 Network_opts.nneurons+1]);
% range = max(x) - min(x);
% ylim([min(x)-0.05*range max(x)+0.05*range]);
% tickpref;
% xtick = 1:Network_opts.nneurons;
% set(gca,'xtick',xtick,'xticklabel', xtick);
% ylabel('IC Value');
% title(sprintf('ICA Analysis\nEnsemble #1'));
% 
% 
% 
% subplot(2,3,6)
% x = icaAssemblyTemplates(:,2);
% if ( abs(max(x)) < abs(min(x)) )  
%     x = -x;
% end
% stem(x,'filled', 'k');
% box off;
% xlim([0 Network_opts.nneurons+1]);
% range = max(x) - min(x);
% ylim([min(x)-0.05*range max(x)+0.05*range]);
% tickpref;
% set(gca,'xtick',xtick,'xticklabel', xtick);
% xlabel('Neuron #');
% ylabel('IC Value');
% title(sprintf('Ensemble #2'));
% 
% 
% 
% figure;
% clf;
% subplot(3,1,1)
% imagesc(Activitymatrix);
% xlim([0 100]);
% tickpref;
% ylabel('Neuron #');
% title('Spike Train Matrix');
% 
% 
% subplot(3,1,2)
% plot(pcaActivities');
% box off;
% xlim([0 100]);
% mxmx = max(max(pcaActivities(:,1:100)));
% mnmn = min(min(pcaActivities(:,1:100)));
% range = mxmx - mnmn;
% ylim([mnmn-0.1*range mxmx+0.1*range]);
% tickpref;
% ylabel('Response');
% title('PCA Activities');
% 
% 
% subplot(3,1,3)
% plot(icaActivities');
% box off;
% xlim([0 100]);
% mxmx = max(max(icaActivities(:,1:100)));
% mnmn = min(min(icaActivities(:,1:100)));
% range = mxmx - mnmn;
% ylim([mnmn-0.1*range mxmx+0.1*range]);
% tickpref;
% title('ICA Activities');
% xlabel('Time bin (X)');
% ylabel('Response');




% All plots on one figure





figure;
clf;
subplot(6,7,[1 8]);
x = 1:Network_opts.nneurons;
index = Assembly_opts.assembly_neurons{1};
y = zeros(size(x));
y(index) = 1;
stem(x,y,'filled', 'k');
box off;
xtick = 1:Network_opts.nneurons;
set(gca,'xtick',xtick,'xticklabel', xtick);
ytick = 0:0.5:1;
set(gca,'ytick',ytick,'yticklabel', ytick);
xlim([0 Network_opts.nneurons+1]);
ylim([-0.25 1.1]);
tickpref;
title(sprintf('Ensemble #1'));
ylabel('Model');


subplot(6,7,[15 22])
x = pcaAssemblyTemplates(:,1);
if ( abs(max(x)) < abs(min(x)) )  
    x = -x;
end
stem(x,'filled', 'k');
box off;
xlim([0 Network_opts.nneurons+1]);
range = max(x) - min(x);
ylim([min([0; x])-0.05*range max(x)+0.05*range]);
tickpref;
xtick = 1:Network_opts.nneurons;
set(gca,'xtick',xtick,'xticklabel', xtick);
ylabel('PCA');



subplot(6,7,[29 36])
x = icaAssemblyTemplates(:,1);
if ( abs(max(x)) < abs(min(x)) )  
    x = -x;
end
stem(x,'filled', 'k');
box off;
xlim([0 Network_opts.nneurons+1]);
range = max(x) - min(x);
ylim([min(x)-0.05*range max(x)+0.05*range]);
tickpref;
xtick = 1:Network_opts.nneurons;
set(gca,'xtick',xtick,'xticklabel', xtick);
ylabel('ICA');
xlabel('Neuron #');






%     % Plot the cell assemblies / independent components
%     subplot(length(cadatacell), 3, (i-1)*3+3);
% 
%     temp = abs(ensembles);
%     temp(temp<0.2) = 0;
%     %temp(temp>=0.2) = 1;
%     imagesc(temp);
%     xlabel('Assembly #');
%     if ( nargin == 3 )
%         ylabel('Neuron #');
%     else
%         tick = 1:length(position);
%         set(gca,'ytick', unique_tick, 'yticklabel', unique_pos);
%         %set(gca,'ytick', tick, 'yticklabel', position);
%         ylabel('Position (um)');
%     end
%     tickpref;
%     set(gca,'xtick', 1:size(nedata.ensembles,2), 'xticklabel', 1:size(nedata.ensembles,2));
%     xlim([0 max(NumberEnsembles)+1]);








subplot(6,7,[2 9]);
x = 1:Network_opts.nneurons;
index = Assembly_opts.assembly_neurons{2};
y = zeros(size(x));
y(index) = 1;
stem(x,y,'filled', 'k');
box off;
xtick = 1:Network_opts.nneurons;
set(gca,'xtick',xtick,'xticklabel', xtick);
ytick = 0:0.5:1;
set(gca,'ytick',ytick,'yticklabel', ytick);
xlim([0 Network_opts.nneurons+1]);
ylim([-0.25 1.1]);
tickpref;
title('Ensemble #2');





subplot(6,7,[16 23]);
x = pcaAssemblyTemplates(:,2);
if ( abs(max(x)) < abs(min(x)) )  
    x = -x;
end
stem(x,'filled', 'k');
box off;
xlim([0 Network_opts.nneurons+1]);
range = max(x) - min(x);
ylim([min([0; x])-0.05*range max(x)+0.05*range]);
tickpref;
set(gca,'xtick',xtick,'xticklabel', xtick);




subplot(6,7,[30 37])
x = icaAssemblyTemplates(:,2);
if ( abs(max(x)) < abs(min(x)) )  
    x = -x;
end
stem(x,'filled', 'k');
box off;
xlim([0 Network_opts.nneurons+1]);
range = max(x) - min(x);
ylim([min(x)-0.05*range max(x)+0.05*range]);
tickpref;
set(gca,'xtick',xtick,'xticklabel', xtick);
xlabel('Neuron #');





subplot(6,7,[3 10]);
imagesc(corrmat)
tickpref;
xlabel('Neuron #');
ylabel('Neuron #');
title('Correlation Matrix');
xtick = 1:Network_opts.nneurons;
ytick = xtick;
set(gca,'xtick',xtick,'xticklabel', xtick);
set(gca,'ytick',ytick,'yticklabel', ytick);


subplot(6,7,[31 38]);
% corrmat(corrmat>0.5*max(max(corrmat))) = 0.5*max(max(corrmat));
% imagesc(corrmat)
% tickpref;
% xlabel('Neuron #');
% ylabel('Neuron #');

hold on;
evalsort = sort(icadata.eigenvalues);
evalnum = 1:length(evalsort);

index = find(evalsort > icadata.lambda_max);
plot(evalnum(index), evalsort(index), 'ko', ...
'markerfacecolor', 'k');

index = find(evalsort <= icadata.lambda_max & evalsort >= icadata.lambda_min);
plot(evalnum(index), evalsort(index), 'o', ...
'markerfacecolor', [0.6 0.6 0.6], ...
'markeredgecolor', [0.6 0.6 0.6]);

index = find(evalsort < icadata.lambda_min);
plot(evalnum(index), evalsort(index), 'ko', ...
'markerfacecolor', 'k');


plot(evalnum, icadata.lambda_max * ones(size(icadata.eigenvalues)), 'r-');
title('E-val dist and M-P boundary');

% plot(evalnum, icadata.lambda_min * ones(size(evalnum)), 'r-');

tickpref;
xtick = 1:Network_opts.nneurons;
set(gca,'xtick',xtick,'xticklabel', xtick);
xlim([0 Network_opts.nneurons+1]);
ylabel('Eigenvalue Mag');
xlabel('Eigenvalue #');








nbins = 100;

subplot(6,7,[5 14])
imagesc(Activitymatrix);
xlim([0 nbins]);
tickpref;
ylabel('Neuron #');
title('Spike Train Matrix');
set(gca,'xticklabel', []);

ytick = 1:Network_opts.nneurons;
set(gca,'ytick',ytick,'yticklabel', ytick);



subplot(6,7,[19 28])
plot(pcaActivities');
box off;
xlim([0 nbins]);
mxmx = max(max(pcaActivities(:,1:100)));
mnmn = min(min(pcaActivities(:,1:100)));
range = mxmx - mnmn;
ylim([mnmn-0.1*range mxmx+0.1*range]);
tickpref;
ylabel('Response');
title('PCA Activities');
set(gca,'xticklabel', []);



subplot(6,7,[33 42])
plot(icaActivities');
box off;
xlim([0 nbins]);
mxmx = max(max(icaActivities(:,1:100)));
mnmn = min(min(icaActivities(:,1:100)));
range = mxmx - mnmn;
ylim([mnmn-0.1*range mxmx+0.1*range]);
tickpref;
title('ICA Activities');
xlabel('Time bin (X)');
ylabel('Response');
legend('Ensemble #1', 'Ensemble #2');



cmap = brewmaps('reds', 21);
cmap = cmap(size(cmap,1):-1:1,:);
cmap = [1 1 1; cmap];
colormap(cmap);


set(gcf,'position', [59 317 1707 623]);




figure;
temp = abs(icaAssemblyTemplates);
temp(temp<0.2) = 0;

imagesc(temp);
cmap = brewmaps('reds', 39);
cmap = flipud(cmap);
cmap = [1 1 1; cmap];
colormap(cmap);
hc = colorbar;
set(hc,'tickdir', 'out');
tickpref;
ylabel('Neuron #');
xlabel('Ensemble #');




return;

    