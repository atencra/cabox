function tmp_trained_control_plot_fio_params(controldir)
% tmp_plot_fio_params Graph nonlinearity curve parameters
%
% fiofit : holds curve fit parameters. 
%
% If fiofit is not supplied, then the current directory is examined and
% *-spk-strfcmb-loc-sta-fio-fit.mat files are loaded and data extracted.
%
%

narginchk(1,1);


[fiodata_trained] = tmp_read_fiofit_params;

[fiodata_control] = tmp_control_dir_read_fiofit_params(controldir);




% Plot the trained data
theta = fiodata_trained.theta;
sigma = fiodata_trained.sigma;
nmse = fiodata_trained.nmse;
r2 = fiodata_trained.r2;

figure;

subplot(3,1,1);
index = find(theta>-5 & theta < 10 & r2 > 0.8);
hist(theta(index),25);
tickpref;
xlabel('Theta');
ylabel('Count');
title(sprintf('Trained: Nonlinearity Fit; N=%.0f',length(index)));

subplot(3,1,2);
index = find(sigma>=0 & sigma < 20 & r2 > 0.8);
hist(sigma(index),25);
tickpref;
xlabel('Sigma');
ylabel('Count');

subplot(3,1,3);
plot(fiodata_trained.nmse, fiodata_trained.r2, 'ko');
xlabel('NMSE');
ylabel('R2');

tickpref;



% Plot the control data
theta = fiodata_control.theta;
sigma = fiodata_control.sigma;
nmse = fiodata_control.nmse;
r2 = fiodata_control.r2;

figure;

subplot(3,1,1);
index = find(theta>-5 & theta < 10 & r2 > 0.8);
hist(theta(index),25);
tickpref;
xlabel('Theta');
ylabel('Count');
title(sprintf('Control: Nonlinearity Fit; N=%.0f',length(index)));

subplot(3,1,2);
index = find(sigma>=0 & sigma < 20 & r2 > 0.8);
hist(sigma(index),25);
tickpref;
xlabel('Sigma');
ylabel('Count');

subplot(3,1,3);
plot(fiodata_control.nmse, fiodata_control.r2, 'ko');
xlabel('NMSE');
ylabel('R2');

tickpref;



return;







function [fiodata] = tmp_read_fiofit_params

if isempty(spkfiofit) % read from files in current directory

    fiofit_matfiles = dir('*-spk-strfcmb-loc-sta-fio-fit.mat');


    position = [];
    theta = [];
    sigma = [];
    nmse = [];
    r2 = [];

    for i = 1:length(fiofit_matfiles)

        infile = fiofit_matfiles(i).name;

        load(infile,'spkfiofit');


        for ii = 1:length(spkfiofit)

            if ( ~isempty(spkfiofit(ii).fiofit) )
                position = [position spkfiofit(ii).position];
                theta = [theta spkfiofit(ii).fiofit.fitParams(2)];
                sigma = [sigma spkfiofit(ii).fiofit.fitParams(3)];
                nmse = [nmse spkfiofit(ii).fiofit.nmse];
                r2 = [r2 spkfiofit(ii).fiofit.r2];
            end

        end % (for ii)

    end % (for i)

end


fiodata.position = position;
fiodata.theta = theta;
fiodata.sigma = sigma;
fiodata.nmse = nmse;
fiodata.r2 = r2;

return;







function [fiodata] = tmp_control_dir_read_fiofit_params(controldir)

% Remember where we are in directory structure
currentdir = pwd;


% Move to directory holding the control data
cd(controldir)


fiofit_matfiles = dir('*-spk-strfcmb-loc-sta-fio-fit.mat');


position = [];
theta = [];
sigma = [];
nmse = [];
r2 = [];

for i = 1:length(fiofit_matfiles)

    infile = fiofit_matfiles(i).name;

    load(infile,'spkfiofit');

    for ii = 1:length(spkfiofit)

        if ( ~isempty(spkfiofit(ii).fiofit) )
            position = [position spkfiofit(ii).position];
            theta = [theta spkfiofit(ii).fiofit.fitParams(2)];
            sigma = [sigma spkfiofit(ii).fiofit.fitParams(3)];
            nmse = [nmse spkfiofit(ii).fiofit.nmse];
            r2 = [r2 spkfiofit(ii).fiofit.r2];
        end

    end % (for ii)

end % (for i)



fiodata.position = position;
fiodata.theta = theta;
fiodata.sigma = sigma;
fiodata.nmse = nmse;
fiodata.r2 = r2;

% Go back to where we started
cd(currentdif);

return;


















% temp = [rtf_params.position];
% [temp, sort_index] = sort(temp);
% rtf_params = rtf_params(sort_index);
% 
% exp = rtf_params(1).exp;
% site = rtf_params(1).site;
% stim = rtf_params(1).stim;
% 
% layer2 = 200;
% layer3 = 400;
% layer4 = 800;
% layer5 = 1100;
% layer6 = 1500;
% whitematter = 2000;
% 
% layer1label = 100;
% layer2label = 300;
% layer3label = 600;
% layer4label = 950;
% layer5label = 1300;
% layer6label = 1750;
% whitematterlabel = 2200;
% 
% 
% figure;
% 
% suptitle(sprintf('%s site%.0f %s', exp, site, stim));
% 
% % First plot the best temporal modulation frequency vs. depth
% subplot(2,5,1);
% hold on;
% plot(best_tmf, position, 'ko', 'markerfacecolor','k');
% % set(gca,'xscale', 'log');
% box on;
% lax = -5;
% rax = 40;
% axis([lax rax 0 2400]);
% plot([lax rax],[layer2 layer2],'k:');
% plot([lax rax],[layer3 layer3],'k:');
% plot([lax rax],[layer4 layer4],'k:');
% plot([lax rax],[layer5 layer5],'k:');
% plot([lax rax],[layer6 layer6],'k:');
% plot([lax rax],[whitematter whitematter],'k-');
% tax = 30;
% text([tax],[layer1label],'L I');
% text([tax],[layer2label],'L II');
% text([tax],[layer3label],'L III');
% text([tax],[layer4label],'L IV');
% text([tax],[layer5label],'L V');
% text([tax],[layer6label],'L VI');
% text([tax],[whitematterlabel],'WM');
% 
% set(gca,'ytick',[0 layer2 layer3 layer4 layer5 layer6 whitematter]);
% xlabel('bTMF (Hz)');
% ylabel('Depth (um)');
% set(gca,'ydir','rev');
% 
% 
% subplot(2,5,2);
% hold on;
% for i = 1:length(position)
%    plot(tmf, 2400 - position(i) + 150*tmtf(i,:) - mean(150*tmtf(i,:)),'k-');
%    p(i) = 2400 - position(i);
% end
% xlabel('TMTF (cyc/s)');
% set(gca,'ytick',(unique(p)),'yticklabel',fliplr(unique(position)));
% axis([0 40 0 2400]);
% % set(gca,'ydir','rev');
% 
% 
% 
% 
% dt = diff(tmf);
% dt = dt(1);
% dtindex = ceil(1/dt);
% 
% subplot(2,5,3);
% imagesc(tmtf);
% xlabel('TMTF (cyc/s)');
% upos = unique(position);
% for i = 1:length(upos)
%    index(i) = min(find(upos(i) == position));
% end
% xtick = [1 8*dtindex+1 8*2*dtindex+1 8*3*dtindex+1 8*4*dtindex+1 length(tmf)];
% set(gca,'xtick', xtick, 'xticklabel', tmf(xtick));
% set(gca,'ytick',index,'yticklabel',upos);
% set(gca,'fontsize', 8)
% 
% 
% 
% % width of temporal mtf vs. depth
% subplot(2,5,4);
% hold on;
% plot(twidth6db, position, 'ko', 'markerfacecolor','k');
% plot(twidth3db, position, 'ro', 'markerfacecolor','r');
% legend('6db','3db');
% box on;
% lax = -0.1 * max(twidth6db);
% rax = 1.1 * max(twidth6db);
% axis([lax rax 0 2400]);
% plot([lax rax],[layer2 layer2],'k:');
% plot([lax rax],[layer3 layer3],'k:');
% plot([lax rax],[layer4 layer4],'k:');
% plot([lax rax],[layer5 layer5],'k:');
% plot([lax rax],[layer6 layer6],'k:');
% plot([lax rax],[whitematter whitematter],'k-');
% set(gca,'ytick',[0 layer2 layer3 layer4 layer5 layer6 whitematter]);
% xlabel('TMTF Width (cyc/s)');
% set(gca,'ydir','rev');
% 
% 
% subplot(2,5,5);
% tposbins = 0:200:2400;
% tnbp = histc([tbp3db 3000], tposbins);
% tnlp = histc([tlp3db 3000], tposbins);
% bar(tposbins, [tnlp(:) tnbp(:)], 'stacked');
% axis([-100 2500 0 max(max([tnbp+tnlp]))+1]);
% set(gca,'xtick', [0 600 1200 1800 2400], 'xticklabel', [0 600 1200 1800 2400]);
% set(gca,'fontsize', 8);
% children = get(gca,'children');
% set(children(1),'facecolor', [0.75 0.75 0.75]);
% set(children(2),'facecolor', [0 0 0]);
% legend('LP3dB', 'BP3dB');
% 
% 
% % plot best spectral modulation frequency vs. depth
% subplot(2,5,6);
% hold on;
% plot(best_xmf, position, 'ko', 'markerfacecolor','k');
% box on;
% lax = -0.5;
% rax = 4;
% axis([lax rax 0 2400]);
% plot([lax rax],[layer2 layer2],'k:');
% plot([lax rax],[layer3 layer3],'k:');
% plot([lax rax],[layer4 layer4],'k:');
% plot([lax rax],[layer5 layer5],'k:');
% plot([lax rax],[layer6 layer6],'k:');
% plot([lax rax],[whitematter whitematter],'k-');
% set(gca,'ytick',[0 layer2 layer3 layer4 layer5 layer6 whitematter]);
% xlabel('bSMF (cyc/oct)');
% ylabel('Depth (um)');
% set(gca,'ydir','rev');
% 
% 
% subplot(2,5,7);
% hold on;
% for i = 1:length(position)
% %    plot(xmf, 150*xmtf(i,:)+position(i)-75,'k-');
%    plot(xmf, 2400 - position(i) + 150*xmtf(i,:) - mean(150*xmtf(i,:)),'k-');
%    p(i) = 2400 - position(i);
% end
% xlabel('XMTF (cyc/oct)');
% set(gca,'ytick',(unique(p)),'yticklabel',fliplr(unique(position)));
% axis([0 4 0 2400]);
% % set(gca,'ydir','rev');
% 
% 
% 
% dx = diff(xmf);
% dx = dx(1);
% dxindex = ceil(1/dx);
% 
% subplot(2,5,8);
% f = round(xmf * 10)/10;
% imagesc(xmtf);
% xlabel('XMTF (cyc/oct)');
% xtick = [1 dxindex+1 2*dxindex+1 3*dxindex+1 length(xmf)];
% set(gca,'xtick', xtick,'xticklabel', [0 1 2 3 4]); %f(xtick) );
% set(gca, 'ytick', index, 'yticklabel', upos);
% set(gca, 'fontsize', 8)
% 
% 
% subplot(2,5,9);
% hold on;
% plot(xwidth6db, position, 'ko', 'markerfacecolor','k');
% plot(xwidth3db, position, 'ro', 'markerfacecolor','r');
% box on;
% lax = -0.1 * max(xwidth6db);
% rax = 1.1 * max(xwidth6db);
% axis([lax rax 0 2400]);
% plot([lax rax],[layer2 layer2],'k:');
% plot([lax rax],[layer3 layer3],'k:');
% plot([lax rax],[layer4 layer4],'k:');
% plot([lax rax],[layer5 layer5],'k:');
% plot([lax rax],[layer6 layer6],'k:');
% plot([lax rax],[whitematter whitematter],'k-');
% set(gca,'ytick',[0 layer2 layer3 layer4 layer5 layer6 whitematter]);
% xlabel('XMTF Width (cyc/oct)');
% set(gca,'ydir','rev');
% 
% 
% subplot(2,5,10);
% xposbins = 0:200:2400;
% xnbp = histc(xbp3db, xposbins);
% xnlp = histc(xlp3db, xposbins);
% bar(xposbins, [xnlp(:) xnbp(:)], 'stacked');
% axis([-100 2500 0 max(max([xnbp+xnlp]))+1]);
% xlabel('Depth (um)');
% set(gca,'xtick', [0 600 1200 1800 2400], 'xticklabel', [0 600 1200 1800 2400]);
% set(gca,'fontsize', 8)
% children = get(gca,'children');
% set(children(1),'facecolor', [0.75 0.75 0.75]);
% set(children(2),'facecolor', [0 0 0]);
% 
% orient landscape;
% 
% print_mfilename(mfilename);
% 
% 









