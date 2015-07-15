function ca_plot_arranged_sta(exp_site_cadata)

cadata = exp_site_cadata.cadata;
probetype = exp_site_cadata.probetype;
pos = cell2mat(cadata.position);
stamat = cadata.stamat;
depth = exp_site_cadata.depth;
nf = cadata.nf;
nlags = cadata.nlags;

if strcmp(probetype,'a1x32-poly3')
    %hard code some probe properties
    offset = 150;
    probelength = 500;
    
    %find min and max depth and difference
    mindepth = min(pos(:,2));
    maxdepth = max(pos(:,2));
    difference = maxdepth - mindepth;
    

    %to correct for lack of single units at the extreme ends of the probe
    while difference < probelength
        
        if mindepth > depth - offset - probelength;
            pos = [0 mindepth-50;pos];
            mindepth = min(pos(:,2));
            difference = maxdepth - mindepth;
                
        elseif maxdepth < depth - offset;
            pos = [pos; 0 maxdepth+50];
            maxdepth = max(pos(:,2));
            difference = maxdepth - mindepth;
        else
            %can't imagine what could go wrong yet but just in case
            error('Something''s wrong!!')
        end
        
    end
    
    relpos = pos;
    relpos(:,2) = relpos(:,2) - mindepth;
    
    %create complete channel position matrix
    temp1 = repmat([-50;0;50],11,1);
    temp2 = sort(repmat((0:50:500)',3,1));
    chanpos = [temp1 temp2];
    chanpos = [chanpos(1,:);chanpos(3:end,:)];
    
    numrep = zeros(size(chanpos,1),1);
    
    for i = 1:size(chanpos,1)        
        index = ismember(relpos,chanpos(i,:),'rows');
        numrep(i) = sum(index);
    end
        
    idx = cumsum(numrep);
    idx = [1 ; idx(1:end-1) + 1];
    
    for ii = 1:max(numrep)
        
        figure;
        unitidx = find(numrep >= ii);
                
        for iii = 1:length(unitidx)
            chan = unitidx(iii);            
            neuron = idx(chan);
                       
            if chan == 1
                subplot(11,3,chan)
            else
                subplot(11,3,chan+1)
            end
            
            rfmat = fliplr(reshape(stamat(neuron,:), nf, nlags));
            plot_strf_symmetric_colormap(rfmat);
        end
        idx = idx + 1;

    end
end


    
% elseif strcmp(probetype,'a1x32-poly2')
%     
%     
% elseif strcmp(probetype,'a1x48-poly2')
%     
    
    


end

