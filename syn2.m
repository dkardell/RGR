%file:   syn2.m

%author:  Dominik A. Kardell
%date:    15 Aug 2020

% Build physical property matrices for 3-layer model with line1A
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

%Sediment depth
for i = 1:size(K_fine,2)
    k = 0;
    for j = size(K_fine,1):-1:1
        if K_fine(j,i) > 1e-15
            zs(i) = size(K_fine,1)-k;
            break
        else
            k = k+1;
        end
    end
end

%%

[rho_crust,por_crust,K_crust,CPs_crust] = deal(zeros(size(K_fine)));

por_crust(1:12,:)   = 0.1;
por_crust(13:end,:) = 0.01;

K_crust(1:12,:)   = 1e-13;
K_crust(13:end,:) = 1e-19;

K_syn = K_fine;
por_syn = por_fine;
rho_syn = rho_fine;
CPs_syn = CPs_fine;

for i = 1:size(K_fine,2)
    count = 1;
    for j = zs(i):-1:1
        K_syn(j,i) = K_crust(count);
        por_syn(j,i) = por_crust(count);
        count = count + 1;
    end
end

por_syn(por_syn < 0) = 0;

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

% save('phys_props_syn2.mat','rho_syn','por_syn','K_syn','CPs_syn','ocean_thick_vector','n_div');
% layer2A = (K_syn == 1e-13);
% save('layer2A.mat','layer2A');