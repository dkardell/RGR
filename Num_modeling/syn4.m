%file:   syn4.m

%author:  Dominik A. Kardell
%date:    15 Aug 2020

% Build physical property matrices for "averaged" crust with line1A
% thickness sediment layer

clear, close all

load phys_props.mat

dx = 37.5/(2^n_div); dy = dx; %[m] Grid spacing

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

K_crust   = zeros(size(K_fine,1)-max(zs),size(K_fine,2));
por_crust = zeros(size(K_fine,1)-max(zs),size(K_fine,2));

for i = 1:size(K_fine,2)
    count = 1;
    for j = zs(i):-1:1
        K_crust(count,i) = K_fine(j,i);
        por_crust(count,i) = por_fine(j,i);
        count = count + 1;
    end
end

K_sed_ave = K_fine(end:-1:end-10,430);
por_sed_ave = por_fine(end:-1:end-10,430);

K_crust_ave = zeros(size(K_crust,1),1)';
por_crust_ave = zeros(size(por_crust,1),1)';
for j = 1:size(K_crust,1)
    K_crust_ave(j) = mean(nonzeros(K_crust(j,:)));
    por_crust_ave(j) = mean(nonzeros(por_crust(j,:)));
end

K_syn = K_fine;
por_syn = por_fine;
rho_syn = rho_fine;
CPs_syn = CPs_fine;

for i = 1:size(K_fine,2)
    count = 1;
    for j = zs(i):-1:1
        K_syn(j,i) = K_crust_ave(count);
        por_syn(j,i) = por_crust_ave(count);
        count = count + 1;
    end
end

por_syn(por_syn <= 0) = 0.00001;

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

% save('phys_props_syn4.mat','rho_syn','por_syn','K_syn','CPs_syn','ocean_thick_vector','n_div');