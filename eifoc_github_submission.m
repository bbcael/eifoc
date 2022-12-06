clear all; close all; clc; % clear workspace
load HyperRrs.mat; % load data
h = Global_Rrs;
wvlend = 301; % use 400-700nm -- n.b. wavelengths above 648nm are <0.0003 Sr-1 >95% of the time, which is the 1km x 1km pixel absolute error for PACE, but up to 700nm is retained. >700nm data are full of missing values or otherwise mostly zero
h = h(:,1:wvlend);
modeend = 6; % plot 6 modes

seawifs_response_function; % coarsen to SeaWiFS equivalent Rrs
k = seawifs;
k = k(1:wvlend,:);
for i = 1:size(h,1); 
    for j = 1:size(k,2);
        s(i,j) = sum(h(i,:)'.*k(:,j))./sum(k(:,j));
    end
end

modis_response_function; % coarsen to MODIS equivalent Rrs
k = modis;
k = k(1:wvlend,:);
for i = 1:size(h,1); 
    for j = 1:size(k,2);
        m(i,j) = sum(h(i,:)'.*k(:,j))./sum(k(:,j));
    end
end

for i = 1:301; % compute broken stick thresholds
    bsh(i) = 100./301*sum(1./(i:301));
end
for i = 1:6;
    bss(i) = 100./6.*sum(1./(i:6));
end
for i = 1:10;
    bsm(i) = 100./10.*sum(1./(i:10));
end

[~,~,~,~,ph,~] = pca((h-nanmean(h))./nanstd(h),'rows','pairwise'); % calculate percent variance explained by each mode
[~,~,~,~,pm,~] = pca((m-nanmean(m))./nanstd(m),'rows','pairwise');
[~,~,~,~,ps,~] = pca((s-nanmean(s))./nanstd(s),'rows','pairwise');

p1 = plot(1:modeend,ph(1:modeend),'linewidth',3,'color',[0.7373    0.5882    0.0275]); % make scree plot
hold on; box on; 
scatter(1:modeend,ph(1:modeend),100,[0.7373    0.5882    0.0275],'filled')
p2 = plot(1:modeend,pm(1:modeend),'--','linewidth',3,'color',[149 83 62]./256);
p4 = plot(1:modeend,ps(1:modeend),':','linewidth',3,'color',(1/256)*[96 171 149]);
scatter(1:modeend,ps(1:modeend),100,(1/256)*[96 171 149],'filled')
p3 = plot(1:modeend,bsh(1:modeend),'k','linewidth',1.25);
scatter(1:modeend,pm(1:modeend),100,[149 83 62]./256,'filled')
set(gca,'ticklabelinterpreter','latex','fontsize',16,'xtick',1:modeend)
axis([0.5 (modeend+.5) 0 70]);
xlabel('component','fontsize',16,'interpreter','latex');
ylabel('percent of variance explained','fontsize',16,'interpreter','latex');
lgnd = legend([p1 p2 p4 p3],'hyperspectral (54, 33, 8, 2)','MODIS (49, 37, 10, 2)','SeaWiFS (63, 28, 8)','significance threshold');
set(lgnd,'interpreter','latex','fontsize',15,'location','northeast')
title('PCA: 191 global spectra, standardised, 400-700nm','fontsize',16,'interpreter','latex')

DoFh_bs = sum(ph>bsh(1:length(ph))') % DoF, hyperspectral
DoFm_bs = sum(pm>bsh(1:10)') % hyperspectral DoF retained, MODIS
DoFs_bs = sum(ps>bsh(1:6)') % hyperspectral DoF retained, SeaWiFS
clear G* i l lgnd p1 p2 prout w

%%

for i = 1:301; % predict each wavelength from multispectral equivalents via linear regression
    [~,~,~,~,stats] = regress(h(:,i),[m ones(191,1)]);
    rmse_m(i) = sqrt(stats(end));
    [~,~,~,~,stats] = regress(h(:,i),[s ones(191,1)]);
    rmse_s(i) = sqrt(stats(end));
end

figure; % plot RMSE vs. wavelength
p1 = plot(400:700,rmse_m,'--','linewidth',3,'color',[149 83 62]./256);
hold on;
p2 = plot(400:700,rmse_s,':','linewidth',3,'color',(1/256)*[96 171 149]);
p3 = plot(400:700,.05.*nanmean(h),'linewidth',3,'color',[0.7373    0.5882    0.0275]); % 5\% significance level
box on; 
axis([400 700 0 7e-5]);
xlabel('wavelength [nm]','fontsize',16,'interpreter','latex');
ylabel('sr$^{-1}$','fontsize',16,'interpreter','latex');
set(gca,'ticklabelinterpreter','latex','fontsize',15)
lgnd = legend([p1 p2 p3],'MODIS RMSE','SeaWiFS RMSE','5\% of mean $R_{rs}(\lambda)$');
set(lgnd,'interpreter','latex','fontsize',15,'location','northeast')

%%

clear all; close all; clc; % clear workspace and load and reshape climatological data
load clim4cael.mat;
load picpoc.mat;
load PSD_slope.mat;
load rrs_seawifs.mat;
rrs412(rrs412>prctile(rrs412(:),99.9)) = NaN;
rrs443(rrs443>prctile(rrs443(:),99.9)) = NaN;
rrs490(rrs490>prctile(rrs490(:),99)) = NaN;
rrs510(rrs510>prctile(rrs510(:),99)) = NaN;
rrs510(rrs510<prctile(rrs510(:),1)) = NaN;
rrs555(rrs555>prctile(rrs555(:),99)) = NaN;
rrs670(rrs670<prctile(rrs670(:),1)) = NaN;
rrs670(rrs670>prctile(rrs670(:),98)) = NaN;
Rrs = [rrs412(:) rrs443(:) rrs490(:) rrs510(:) rrs555(:) rrs670(:)];
[~,~,~,~,pr,~] = pca((Rrs-nanmean(Rrs))./nanstd(Rrs),'rows','pairwise');

for i = 1:6; % generate SeaWiFS broken stick threshold
    bs6(i) = 100./6.*sum(1./(i:6));
end

% n.b. adding CAFE NPP, CbPMv2 NPP, f_micro, bbp:chl, removing zeros and negative data, removing outliers, or log-transforming do not change the result
Prods = [bbp(:) pic(:) poc(:) Xi(:) exp(logChl(:)) z_eu(:)];
[A,~,C,~,pp,~] = pca((Prods-nanmean(Prods))./nanstd(Prods),'rows','pairwise');

p1 = plot(1:6,pr,'linewidth',3,'color',[87 104 81]./217); % make scree plot
hold on; box on; 
scatter(1:6,pr,100,[87 104 81]./217,'filled')
p2 = plot(1:6,pp(1:6),'--','linewidth',3,'color',[128 40 97]./256);
p3 = plot(1:6,bs6,'k','linewidth',1.25);
scatter(1:6,pp(1:6),100,[128 40 97]./256,'filled')
set(gca,'ticklabelinterpreter','latex','fontsize',16,'xtick',1:6)
axis([0.5 6.5 0 60]);
xlabel('component','fontsize',16,'interpreter','latex');
ylabel('percent of variance explained','fontsize',16,'interpreter','latex');
lgnd = legend([p1 p2 p3],'reflectance (6 SeaWiFS wavelengths)','products (Chl, POC, PIC, C$_{phyto}$, $\xi$, Z$_{eu}$)','significance threshold ($n=6$)');
set(lgnd,'interpreter','latex','fontsize',15,'location','northeast')
title('PCA: reflectance and product climatologies, standardised','fontsize',16,'interpreter','latex')