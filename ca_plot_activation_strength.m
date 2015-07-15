function [Patterns,contribution] = ca_plot_activation_strength(varargin)
% ca_plot_activation_strength Neuron contribution to CA activation strength
% 
%     [Patterns,contribution] = ca_plot_activation_strength(cadata)
%     --------------------------------------------------------------
%
%     cadata : struct holding the spike train matrix, the position of each
%     neuron, the cell assembly patterns, and the cell assembly activities.
%
%     Example cadata:
% 
%     cadata = 
% 
%           spktrain: [23x120265 double]
%              fsdvd: 96000
%                 df: 480
%           position: {1x23 cell}
%           Patterns: [23x6 double]
%         Activities: [6x120265 double]
%                 nf: 64
%              nlags: 20
%             stamat: [23x1280 double]
%          ca_stamat: [6x1280 double]
% 
%     Since it takes time to calculate the activation strengths of each neuron,
%     this function can also be called by:
% 
%     ca_plot_activation_strength(Patterns, contribution);
% 
%     Here Patterns and contribution would have been output variables from a 
%     previous call to this function.
% 
%     This function will know that if there is one input argument, then the input
%     is cadata.
% 
%     If there are two input arguments, then they are Patterns and contribution.




if nargin == 1
    cadata = varargin{1};
end

if nargin == 2
    Patterns = varargin{1};
    contribution = varargin{2};
end


if nargin == 1;
    [Patterns, contribution] = ca_calc_assembly_contribution(cadata);
end


ca_plot_assembly_contribution(Patterns, contribution);

return;




function [Patterns, contribution] = ca_calc_assembly_contribution(cadata)

spktrain_matrix = cadata.spktrain;
Patterns = cadata.Patterns;
Activities = cadata.Activities;

nCells = size(spktrain_matrix, 1);
nPatterns = size(Patterns, 2);

contribution = zeros(nCells, nPatterns);
mean_activity_total = mean(Activities,2);

for i = 1:nCells
tic
    fprintf('Cell #%.0f\n', i);

    spktrain_matrix_subset = spktrain_matrix;
    spktrain_matrix_subset(i,:) = 0;

    activity_subset = assembly_activity(Patterns, spktrain_matrix_subset);
    mean_activity_subset = mean(activity_subset,2);


    ratio = mean_activity_subset ./ mean_activity_total;
    strength = 1 - ratio;

    contribution(i,:) = strength(:)';
toc
end % (for ii)

return;




function ca_plot_assembly_contribution(Patterns, contribution)

nPatterns = size(Patterns, 2);

% Plot activation strength for each cell against ICA weights
% Make plot for each cell assembly
if nPatterns == 1
    nrows = 1;
    ncols = 1;
elseif nPatterns == 2
    nrows = 1;
    ncols = 2;
elseif nPatterns <= 4
    nrows = 2;
    ncols = 2;
elseif nPatterns <= 6
    nrows = 2;
    ncols = 3;
elseif nPatterns <= 9
    nrows = 3;
    ncols = 3;
elseif nPatterns <= 12
    nrows = 3;
    ncols = 4;
elseif nPatterns <= 16
    nrows = 4;
    ncols = 4;
elseif nPatterns <= 20
    nrows = 4;
    ncols = 5;
elseif nPatterns <= 25
    nrows = 5;
    ncols = 5;
else
    error('Wow! That is a lot of CAs!');
end


figure;
for i = 1:nPatterns
    subplot(nrows, ncols, i);
    plot(Patterns(:,i), contribution(:,i), 'ko', 'markerfacecolor', [0.75 0.75 0.75]);
    tickpref;
    xlabel('CA weights');
    ylabel('Contribution');
    title(sprintf('Cell contributions to CA #%.0f', i));
end % (for i)




cmb = combnk(1:6,2);
ncmb = size(cmb,1);

if ncmb == 1
    nrows = 1;
    ncols = 1;
elseif ncmb == 2
    nrows = 1;
    ncols = 2;
elseif ncmb <= 4
    nrows = 2;
    ncols = 2;
elseif ncmb <= 6
    nrows = 2;
    ncols = 3;
elseif ncmb <= 9
    nrows = 3;
    ncols = 3;
elseif ncmb <= 12
    nrows = 3;
    ncols = 4;
elseif ncmb <= 16
    nrows = 4;
    ncols = 4;
elseif ncmb <= 20
    nrows = 4;
    ncols = 5;
elseif ncmb <= 25
    nrows = 5;
    ncols = 5;
else
    error('Wow! That is a lot of CAs!');
end

figure;
for i = 1:size(cmb,1)
    subplot(nrows, ncols, i);
    plot(contribution(:,cmb(i,1) ), contribution(:, cmb(i,2) ), 'ko', 'markerfacecolor', [0.75 0.75 0.75]);
    tickpref;
    xlabel(sprintf('Contribution to CA #%.0f', cmb(i,1) ));
    ylabel(sprintf('Contribution to CA #%.0f', cmb(i,2) ));
end % (for i)


return;


