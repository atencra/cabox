function [] = ca_batch_spet2train(stim, filetype, fsd)
% spet2train - Converts a Spike Event Time Array to a sampled impulse array
%
% [train] = spet2train(spet, fsspk, fsd)
%
% spet  : Spike Event Time Array, in units of sample number
% fsspk : Sampling Rate of spet
% fsd   : Desired Sampling Rate for outpute spike train
%         If fsd = 1000, spet is created using 1 ms bins
%
% train : binned spike train


for i = 1:length(stim)
    
    file = dir( sprintf('*-site*-*um-*db-%s-fs*-*-%s.mat',stim{i},filetype) );
    
    if length(file) > 1
        error('Specify more specific filetype.')
    elseif isempty(file)
        error('No such filetype.')
    end
        
    load(file.name)
    spktimes = {spk.spiketimes};
    lastspk = max(cellfun(@max,spktimes));
    spktrain = cell2mat(cellfun(@(x) histc(x,1:1/fsd*1000:lastspk),spktimes','UniformOutput',0));
    
    dashidx = strfind(file.name,'-');
    outfile = sprintf('%sspktrain-fs%s.mat',file.name(1:dashidx(6)),num2str(fsd));
    save(outfile,'spktrain')

end
    
    