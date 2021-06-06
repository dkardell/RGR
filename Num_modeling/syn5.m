%file:   syn5.m

%author:  Dominik A. Kardell
%date:    15 Aug 2020

% Build physical property matrices for real crust with line1A/real
% thickness sediment layer

clear, close all

load phys_props.mat

load /disk/student/dkardell/FWI/CREST/line1A/MCMC_geo/bath_line1A_hydro.txt
wb = bath_line1A_hydro(:,2); 
ocean_thick_vector = wb(1:3:end)*1000;

dx = 37.5/(2^n_div); dy = dx; %[m] Grid spacing

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

K_crust   = repmat(min(K_fine),size(K_fine,1));
por_crust = zeros(size(K_fine,1));

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

K_syn = K_fine;
por_syn = por_fine;
rho_syn = rho_fine;
CPs_syn = CPs_fine;

for i = 1:size(K_syn,2)
    K_syn(end:-1:end-9,i) = K_sed_ave(10:-1:1);
    por_syn(end:-1:end-9,i) = por_sed_ave(10:-1:1);
    rho_syn(end:-1:end-9,i) = rho_fine(end,1);
    CPs_syn(end:-1:end-9,i) = CPs_fine(end,1);
    count = 1;
    for j = size(K_syn,1)-10:-1:1
        K_syn(j,i) = K_crust(count,i);
        por_syn(j,i) = por_crust(count,i);
        rho_syn(j,i) = rho_fine(1,1);
        CPs_syn(j,i) = CPs_fine(1,1);
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

blue1 = [0.71,0.886,0.941];
blue2 = [0.106,0.353,0.655];
c1 = interp1([1 65],[blue1(1) blue2(1)],1:65)';
c2 = interp1([1 65],[blue1(2) blue2(2)],1:65)';
c3 = interp1([1 65],[blue1(3) blue2(3)],1:65)';
cbar = horzcat(c1,c2,c3);

figure 
contourf(-log10(K_syn), 100, 'LineStyle', 'none')
logvec = fliplr(logspace(-10,-20,6));
set(gca, 'colormap', cbar) 
caxis([13 19])
colorbar('YTick', -log10(fliplr(logvec)), 'YTickLabel', fliplr(logvec))
title('Permeability [m^2]')

%%

% save('phys_props_syn5.mat','rho_syn','por_syn','K_syn','CPs_syn','ocean_thick_vector','n_div');