function [cadata] = ca_calc_surrogate_spktrain_from_sta_fio(spk, stimulus, trigger, FsDVD, totalDF)
% ca_calc_cell_assembly_sta_from_spk_stimulus_trigger Estimate STAs from PCA/ICA cell assembly analysis
% 
%     [cadata] = ca_calc_surrogate_sta_fio_spktrain(spk, stimulus, trigger, FsDVD, totalDF)
%     ---------------------------------------------------------------------------------------------
%     spk : struct array of spike times. Obtaines from saveSpikeSort.m
%
%     stimulus : ripple stimulus in matrix format. Each row is one
%     frequency, and each column is one time bin.
%
%     trigger : trigger for ripple stimulus
%
%     FsDVD : sampling of DVD Audio system
%
%     totalDF : The total downsampling factor for the stimulus matrix. 
%     A scalar. This will be the downsampling factor for the original .spr 
%     file, as well as any additional downsampling. The original 
%     downsampling will usually be 48, and additional downsampling will be 
%     [5, 10, 20, 40, or 100]. 
%
%     Therefore, the totalDF would be 48 * [5, 10, 20, 40, or 100], which
%     is equivalent to 5ms, 10ms, etc. for 96 kHz sound stimulus sampling
%     rate.
%
%     cadata : struct holding calculations. Has the form:
%         cadata.spktrain = binned spike train, one row per neuron
%         cadata.fsdvd = sampling rate of sound stimulus
%         cadata.df = total downsampling factor of ripple noise envelope
%         cadata.position = recording depth. position(i) -> spktrain(i,:)
%         cadata.Patterns = cell assemblies
%         cadata.Activities = time course of cell assembly activity; one row per assembly
%         cadata.nf = number of frequencies in STA
%         cadata.nlags = number of time lags in STA
%         cadata.stamat = matrix of spike train STAs.
%         cadata.ca_stamat = matrix of cell assembly STAs.


ca;

% Make sure there are 5 input arguments.
narginchk(5,5);


% Get number of lags for STA to cover 100 ms
if ( (totalDF / FsDVD * 1000) < 5 )
    nlags = ceil( 100 / (totalDF / FsDVD * 1000) ); % 100 ms / time bin size
else
    nlags = 20;
end

% Number of frequencies
nf = size(stimulus, 1);



% Get spike train matrix.
[spktrain_matrix, position] = ...
    ca_create_spktrain_matrix_from_stimulus_spk(spk, stimulus, trigger, FsDVD, totalDF);


% Find cell assemblies
fprintf('\nDetecting Cell Assemblies\n');
[Patterns, Activities] = ca_detect_cell_assemblies_data(spktrain_matrix);

% Plot cell assembly results
ca_plot_cell_assemblies_data(spktrain_matrix, Patterns, Activities, position);


% Calculate STA for each cell assembly. Each row of Activities holds the
% activations for one cell assembly.
ca_stamat = zeros(size(Activities,1),size(stimulus,1)*nlags);
for i = 1:size(Activities, 1)
    sta = ca_calc_sta_from_stimulus_spktrain(stimulus, Activities(i,:), nlags);
    ca_stamat(i,:) = sta(:)';
end







% Calculate STA for each neuron. Each row of spktrain_matrix holds the
% spike train for one neuron.
Nbins = 15; % histogram size of nonlinearity.
stamat = zeros(size(spktrain_matrix,1),size(stimulus,1)*nlags);
stamatPred = zeros(size(spktrain_matrix,1),size(stimulus,1)*nlags);
spktrainPred = zeros(size(spktrain_matrix));

for i = 1:size(spktrain_matrix, 1)
    locator = spktrain_matrix(i,:);

    sta = ca_calc_sta_from_stimulus_spktrain(stimulus, locator, nlags);


    % Calculate projection values
    [xprior, ~] = ca_calc_projection_from_sta_stimulus_spktrain(sta, stimulus, locator);


    % Calculate nonlinearity - p(spk|x)
    [x, px, pxspk, pspk, pspkx] = ca_calc_fio_prob_from_proj_spktrain(xprior, locator, Nbins);

    
    % With nonlinearity, convert projection values to probabilities
    probPred = interp1(x, pspkx, xprior, 'linear', 0);

    % Predict spike train assuming no noise - probably not Poisson
    thr = linspace(0, max(probPred), 500);
    [locatorPred] = noiseFreeSpikePrediction(probPred, locator, thr);

    % Predict spike train assuming no noise - Poisson
    thr = linspace(0, 1, 500);
    [locatorPredNoisy] = noiseFreeSpikePrediction(probPred, locator, thr);

    percent = 0.5;
    [locatorRand, ratio] = noiseFreeAddSpikesSpikePrediction(locator, locatorPred, percent)

    [sum(locator) sum(locatorPred) sum(locatorPredNoisy) sum(locatorRand)]

    


    % For comparison, calculate STA from predicted spikes
    staPred = ca_calc_sta_from_stimulus_spktrain(stimulus, locatorPred, nlags);

    staPredNoisy = ca_calc_sta_from_stimulus_spktrain(stimulus, locatorPredNoisy, nlags);

    staPredRand = ca_calc_sta_from_stimulus_spktrain(stimulus, locatorRand, nlags);


    plotSurrogateSpktrainSTA(sta, staPred, staPredRand);
    plotSurrogateSpktrainISI(locator, locatorPred, locatorRand);

pause


    % Sanity check: plot original STA/nonlinearity and the STA from
    % surrogate spike train
%     hf = figure;
%     subplot(2,2,1);
%     imagesc(sta);
%     title('Original STA');
% 
%     subplot(2,2,3);
%     plot(x,pspkx,'ko-');
%     xlabel('Projection Value');
%     ylabel('Prob(spike|proj)');
%     title('Original Nonlinearity');
% 
%     subplot(2,2,2);
%     imagesc(staPred);
%     title('Predicted STA');
%     pause(1);
%     close(hf);

    stamat(i,:) = sta(:)';

    stamatPred(i,:) = staPred(:)';
    stamatPredNoisy(i,:) = staPredNoisy(:)';
    stamatPredRand(i,:) = staPredRand(:)';

    spktrainPred(i,:) = locatorPred;
    spktrainPredNoisy(i,:) = locatorPredNoisy;
    spktrainPredRand(i,:) = locatorRand;

end



% Find cell assemblies for Predicted Spike Trains
fprintf('\nDetecting Predicted Spike Train Cell Assemblies\n');
[PatternsPred, ActivitiesPred] = ca_detect_cell_assemblies_data(spktrainPred);
ca_plot_cell_assemblies_data(spktrainPred, PatternsPred, ActivitiesPred, position);





% For predicted spike trains, calculate the STA for each cell assembly. 
ca_staMatPred = zeros(size(ActivitiesPred,1),size(stimulus,1)*nlags);
for i = 1:size(ActivitiesPred, 1)
    sta = ca_calc_sta_from_stimulus_spktrain(stimulus, ActivitiesPred(i,:), nlags);
    ca_staMatPred(i,:) = sta(:)';
end



% Assign actual and predicted data to struct
cadata.spktrain = spktrain_matrix;
cadata.fsdvd = FsDVD;
cadata.df = totalDF;
cadata.position = position;
cadata.Patterns = Patterns;
cadata.Activities = Activities;

cadata.nf = nf;
cadata.nlags = nlags;
cadata.stamat = stamat;
cadata.ca_stamat = ca_stamat;

cadata.spktrainPred = spktrainPred;
cadata.PatternsPred = PatternsPred;
cadata.ActivitiesPred = ActivitiesPred;
cadata.stamatPred = stamatPred;
cadata.ca_stamatPred = ca_stamatPred;

return;



function plotSurrogateSpktrainSTA(sta, staPred, staPredNoisy)

% Sanity check: plot original STA and the STAs from surrogate spike trains

figure;
subplot(1,3,1);
imagesc(sta);
title('Original STA');

subplot(1,3,2);
imagesc(staPred);
title('Predicted STA - WithOut Noise');

subplot(1,3,3);
imagesc(staPredNoisy);
title('Predicted STA - With Noise');

return


function plotSurrogateSpktrainISI(locator, locatorPred, locatorPredNoisy)

% Sanity check: plot original ISI and the ISIs from surrogate spike trains

edges = 0:60:3000;
indexOrig = find(locator > 0); % pseudo spike times
isiOrig = diff(indexOrig);
nOrig = histc(isiOrig, edges);

indexPred = find(locatorPred > 0); % pseudo spike times
isiPred = diff(indexPred);
nPred = histc(isiPred, edges);

indexPredNoisy = find(locatorPredNoisy > 0); % pseudo spike times
isiPredNoisy = diff(indexPredNoisy);
nPredNoisy = histc(isiPredNoisy, edges);

maxmax = max([max(nOrig) max(nPred) max(nPredNoisy)]);

figure;
subplot(1,3,1);
bar(edges,nOrig,'histc');
title('ISI Original');
xlim([0 max(edges)]);
ylim([0 maxmax]);

subplot(1,3,2);
bar(edges,nPred,'histc');
title('ISI - WithOut Noise');
xlim([0 max(edges)]);
ylim([0 maxmax]);

subplot(1,3,3);
bar(edges,nPredNoisy,'histc');
title('ISI - With Noise');
xlim([0 max(edges)]);
ylim([0 maxmax]);

set(gcf,'position', [225 448 1404 354]);

[h,p] = kstest2(isiPred,isiPredNoisy);
fprintf('Predicted ISI: with vs without noise: p = %.5f\n', p);


return



function [locatorRand, ratio] = noiseFreeAddSpikesSpikePrediction(locator, locatorPred, percent)


% Find spikes in spike train
index = 1:length(locator);

indexSpikes = find(locatorPred > 0); % index for spikes

[locatorPred(indexSpikes(:)) indexSpikes(:)];

length(indexSpikes);


% Bins without spikes
indexOpen = setdiff(index, indexSpikes);

% Random subset of bins without spikes
ip = randperm(length(indexOpen), ceil(percent*sum(locator)));
indexRandSpk = sort(indexOpen(ip)); % index for random spikes

% Indices for real spikes and random spikes
indexTotal = [indexSpikes(:)' indexRandSpk(:)'];

% Randomly select from total set of indices
ip = randperm(length(indexTotal), sum(locator));

% Assign spikes based on random sampling of indices
locatorRand = zeros(size(locator));
locatorRand(indexTotal(ip)) = 1;

ratio = length(indexSpikes) / length(indexTotal);

return;



function [locatorPred] = noiseFreeSpikePrediction(probPred, locator, thr)

    % Noise-free analysis
    % Try multiple thresholds to find the one that matches the actual spike count
%     thr = linspace(0, max(probPred), 100);
    nPredSpk = zeros(size(thr));
    for j = 1:length(thr)
        predSpk = probPred > thr(j);
        nPredSpk(j) = sum(predSpk);
    end % (for j)

    % Find threshold that matches overall spike count of original train
    d = abs(nPredSpk - sum(locator));
    indMin = find(d == min(d),1);
    locatorPred = probPred > thr(indMin);

return;


function [locatorPred] = noisySpikePrediction(probPred, locator, thr)

    % Add more random-ness
    % Try multiple thresholds to find the one that matches the actual spike count
%     thr = linspace(0, 1, 100);

    nPredSpk = zeros(size(thr));
    rng(1);
    randVals = rand(size(probPred,1),1);
    for j = 1:length(thr)
        predSpk = (randVals+thr(j)) < probPred;
        nPredSpk(j) = sum(predSpk);
    end % (for j)

    % Find threshold that matches overall spike count of original train
    d = abs(nPredSpk - sum(locator));
    indMin = find(d == min(d),1);
    locatorPred = probPred > thr(indMin);

return;









% Extra previous code
%
%     % Noise-free analysis
%     % Try multiple thresholds to find the one that matches the actual spike count
%     thr = linspace(0, max(probPred), 100);
%     nPredSpk = zeros(size(thr));
%     for j = 1:length(thr)
%         predSpk = probPred > thr(j);
%         nPredSpk(j) = sum(predSpk);
%     end % (for j)
% 
%     % Find threshold that matches overall spike count of original train
%     d = abs(nPredSpk - sum(locator));
%     indMin = find(d == min(d),1);
%     locatorPred = probPred > thr(indMin);
% 
% 
%     % Add more random-ness
%     % Try multiple thresholds to find the one that matches the actual spike count
%     thr = linspace(0, 1, 100);
%     nPredSpk = zeros(size(thr));
%     for j = 1:length(thr)
%         predSpk = (rand(size(probPred,1),1)+thr(j)) < probPred;
%         predSpk = probPred > thr(j);
%         nPredSpk(j) = sum(predSpk);
%     end % (for j)
% 
%     % Find threshold that matches overall spike count of original train
%     d = abs(nPredSpk - sum(locator));
%     indMin = find(d == min(d),1);
%     locatorPred = probPred > thr(indMin);

