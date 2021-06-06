%file:   syn1.m

%author:  Dominik A. Kardell
%date:    15 Aug 2020

% Build physical property matrices for 3-layer model with constant
% thickness sediment layer

clear, close all

load phys_props.mat

load /disk/student/dkardell/FWI/CREST/line1A/MCMC_geo/bath_line1A_hydro.txt
wb = bath_line1A_hydro(:,2); 
ocean_thick_vector = wb(1:3:end)*1000;

figure(1)
subplot(4,1,1)
contourf(por_fine, 100, 'LineStyle', 'none')
jetf = flipud(jet);
set(gca, 'colormap', jetf); colorbar
caxis([0 0.3])
title('Porosity from Carlson (2014) [frac]')

subplot(4,1,2)
contourf(-log10(K_fine), 100, 'LineStyle', 'none')
logvec = fliplr(logspace(-10,-20,6));
set(gca, 'colormap', jet) 
caxis([10 20])
colorbar('YTick', -log10(fliplr(logvec)), 'YTickLabel', fliplr(logvec))
title('Permeability [m^2]')

subplot(4,1,3)
contourf(rho_fine, 100, 'LineStyle', 'none')
set(gca, 'colormap', jet); colorbar
caxis([1 3])
title('Matrix Density [g/cm^3]')

subplot(4,1,4)
contourf(CPs_fine, 100, 'LineStyle', 'none')
set(gca, 'colormap', jetf); colorbar
%caxis([1 3])
title('Specific Heat Capacity [J/(kg*K)]')

%%

[rho_syn,por_syn,K_syn,CPs_syn] = deal(zeros(size(K_fine)));

rho_syn(end-9:end,:)   = 2.71;
rho_syn(1:end-10,:)     = 3.33;

por_syn(end-9:end,:)    = 0.6;
por_syn(end-21:end-10,:) = 0.1;
por_syn(1:end-22,:)      = 0.01;

K_syn(end-9:end,:)    = 1e-17;
K_syn(end-21:end-10,:) = 1e-13;
K_syn(1:end-22,:)      = 1e-19;

CPs_syn(end-9:end,:) = 1004;
CPs_syn(1:end-10,:)   = 1160;

figure(2)
subplot(4,1,1)
contourf(por_syn, 100, 'LineStyle', 'none')
jetf = flipud(jet);
set(gca, 'colormap', jetf); colorbar
caxis([0 0.3])
title('Porosity [frac]')

subplot(4,1,2)
contourf(-log10(K_syn), 100, 'LineStyle', 'none')
logvec = fliplr(logspace(-10,-20,6));
set(gca, 'colormap', jet) 
caxis([10 20])
colorbar('YTick', -log10(fliplr(logvec)), 'YTickLabel', fliplr(logvec))
title('Permeability [m^2]')

subplot(4,1,3)
contourf(rho_syn, 100, 'LineStyle', 'none')
set(gca, 'colormap', jet); colorbar
caxis([1 3])
title('Matrix Density [g/cm^3]')

subplot(4,1,4)
contourf(CPs_syn, 100, 'LineStyle', 'none')
set(gca, 'colormap', jetf); colorbar
%caxis([1 3])
title('Specific Heat Capacity [J/(kg*K)]')

%%

% save('phys_props_syn1.mat','rho_syn','por_syn','K_syn','CPs_syn','ocean_thick_vector','n_div');