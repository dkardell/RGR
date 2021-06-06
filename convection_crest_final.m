%file:   convection_crest.m

%author:  Dominik A. Kardell, modified from Evan J. Ramos
%date:    07 June 2020

clear, clc, close all

load phys_props.mat %Import physical property grids
dx = 37.5/(2^n_div); dy = dx; %[m] Grid spacing

load /disk/student/dkardell/FWI/CREST/line1A/MCMC_geo/bath_line1A_hydro.txt
load layer2A.mat
wb = bath_line1A_hydro(:,2); 
ocean_thick_vector = wb(1:3:end)*1000;

por_fine(por_fine <= 0) = 0.00001;

temp = 0; %0 for constant bottom T, 1 for variable bottom T accounting for crustal thickness

%% ??

% blue1 = [0.71,0.886,0.941];
% blue2 = [0.106,0.353,0.655];
% c1 = interp1([1 65],[blue1(1) blue2(1)],1:65)';
% c2 = interp1([1 65],[blue1(2) blue2(2)],1:65)';
% c3 = interp1([1 65],[blue1(3) blue2(3)],1:65)';
% cbar = horzcat(c1,c2,c3);
% 
% figure 
% contourf(-log10(K_fine), 100, 'LineStyle', 'none')
% logvec = fliplr(logspace(-10,-20,6));
% set(gca, 'colormap', cbar) 
% caxis([13 19])
% colorbar('YTick', -log10(fliplr(logvec)), 'YTickLabel', fliplr(logvec))
% title('Permeability [m^2]')

%% Dimensional constants

%Space
Depth = size(K_fine,1)*dy;
Length = size(K_fine,2)*dx;

%Time
yrs2s    = 3600*24*365; %conversion factor to seconds from years
endtime  = 50000*yrs2s; %20000
t_totals = 10000;         %5000 Total number of time steps saved
dt       = endtime/t_totals; %[s] time step
times    = (1:t_totals)*endtime/t_totals;

%% Build Grid and operators

Nx = size(K_fine,2); 
Ny = size(K_fine,1);

% Defining the computational grid and conductivities
Grid.xmin = 0; Grid.xmax = Length;  Grid.Nx = Nx;
Grid.ymin = 0; Grid.ymax = Depth;   Grid.Ny = Ny; 
Grid.Nz = 1;
Grid.psi_dir = 'xy';
Grid      = build_grid(Grid);
[Xc,Yc]   = meshgrid(Grid.xc,Grid.yc);
[Xf,Yf]   = meshgrid(Grid.xf,Grid.yf);
[Xfx,Yfx] = meshgrid(Grid.xf,Grid.yc);
[Xfy,Yfy] = meshgrid(Grid.xc,Grid.yf);

%Operators
[D,G,I]           = build_ops(Grid);

ocean_thick_vector = linspace(mean(ocean_thick_vector)+0.5,mean(ocean_thick_vector)-0.5,Nx)';
ocean_thick_matrix = repmat(ocean_thick_vector',Ny,1);

for i = 1:size(K_fine,2)
    for j = size(K_fine,1):-1:1
        if K_fine(j,i) > 1e-15
            cr(i) = j;
            break;
        end
    end
end
% T_bot = (273.15+63+(cr*dy/1000-2)*31)'; %355.15;         %[K] basal temperature. Linearly interpolating between 54 C at 2 km and 82 C at 3 km

% if temp == 0
%     T_bot = linspace(mean(T_bot),mean(T_bot),Nx)';
% end

%% Smooth data and plot

% sy = 7; %Smoothing window size in x-direction (grid points)
% sx = 7; %Smoothing window size in y-direction (grid points)
% 
% K_smooth = smooth2a(K_fine,sy,sx);
% Ks_smooth = smooth2a(Ks_fine,sy,sx);
% por_smooth = smooth2a(por_fine,sy,sx);
% CPs_smooth = smooth2a(CPs_fine,sy,sx);
% rho_smooth = smooth2a(rho_fine,sy,sx);
% % K_smooth = smoothdata(K_fine,2,'movmean',30);
% % por_smooth = smoothdata(por_fine,2,'movmean',30);
% % CPs_smooth = smoothdata(CPs_fine,2,'movmean',30);
% % rho_smooth = smoothdata(rho_fine,2,'movmean',30);
% % K_smooth = smoothdata(K_smooth,1,'gaussian',10);
% % por_smooth = smoothdata(por_smooth,1,'gaussian',10);
% % CPs_smooth = smoothdata(CPs_smooth,1,'movmean',10);
% % rho_smooth = smoothdata(rho_smooth,1,'movmean',10);
% 
% figure(2)
% set(gcf,'units','normalized','outerposition',[0 0 1 0.9])
% subplot(4,2,1)
% contourf(por_fine, 100, 'LineStyle', 'none'); %axis equal;
% jetf = flipud(jet);
% set(gca, 'colormap', jetf); colorbar
% caxis([0 0.3])
% title('Porosity [frac]')
% xlabel('x [m]'); ylabel('z [m]');
% 
% subplot(4,2,3)
% contourf(Xc,Yc,log10(K_fine),100,'LineStyle','none'); %axis equal;
% colormap(gca,jetf); colorbar 
% caxis([-20 -10])
% logvec = fliplr(logspace(-10,-20,6));
% colorbar('YTick', log10(logvec), 'YTickLabel', log10(logvec))
% title('Permeability [m^2]')
% xlabel('x [m]'); ylabel('z [m]');
% 
% subplot(4,2,5)
% contourf(Xc,Yc,CPs_fine,'LineStyle','none'); colorbar; %axis equal 
% xlabel('x (m)'); ylabel('z (m)');
% title('Specific heat');
% 
% subplot(4,2,7)
% contourf(Xc,Yc,rho_fine,'LineStyle','none'); colorbar; %axis equal; 
% xlabel('x (m)'); ylabel('z (m)');
% title('Density');
% 
% subplot(4,2,2)
% contourf(por_smooth, 100, 'LineStyle', 'none'); %axis equal;
% set(gca, 'colormap', jetf); colorbar
% caxis([0 0.3])
% title('Smoothed Porosity [frac]')
% xlabel('x [m]'); ylabel('z [m]');
% 
% subplot(4,2,4)
% contourf(Xc,Yc,log10(K_smooth),100,'LineStyle','none'); %axis equal;
% colormap(gca,jetf); colorbar 
% caxis([-20 -10])
% colorbar('YTick', log10(logvec), 'YTickLabel', log10(logvec))
% title('Smoothed Permeability [m^2]')
% xlabel('x [m]'); ylabel('z [m]');
% 
% subplot(4,2,6)
% contourf(Xc,Yc,CPs_smooth,'LineStyle','none'); colorbar; %axis equal; 
% xlabel('x (m)'); ylabel('z (m)');
% title('Smoothed specific heat');
% 
% subplot(4,2,8)
% contourf(Xc,Yc,rho_smooth,'LineStyle','none'); colorbar; %axis equal; 
% xlabel('x (m)'); ylabel('z (m)');
% title('Density');
% 
% K_fine = K_smooth;
% Ks_fine = Ks_smooth;
% por_fine = por_smooth;
% CPs_fine = CPs_smooth;
% rho_fine = rho_smooth;

%%

%Darcy
Mu       = 5e-5;      %[Pa*s] Dynamic viscosity of water at 20 deg C %%%%%%%%%%
rho_f    = 1030;        %[kg/m^3] Density of water
rho_s    = rho_fine;
grav     = 9.81;        %[m/s^2] Acceleration due to gravity on Earth
n        = 3;        %integer exponent for power law calc
tau      = sqrt(2);     %tortuosity of porosity medium, assuming packed spheres geometry
kvkh     = 1;

%Temperature
T_s      = 275.15;         %[K] surface temp
T_s_vec  = T_s*ones(Grid.Nx,1);
T_bot    = 356.69 * ones(Nx,1);
kf       = 0.65;        %[W/m/K] thermal conductivity of water
ks       = 2.73;        %[W/m/K] thermal conductivity of rock
cp_f     = 4157;        %[J/(kg K)] specific heat of water
% cp_s     = 1171;


%% Create Medium

phi = por_fine; %por_smooth;
k = K_fine; % k_smooth;
Kd = comp_mean(k,1,1,Grid);
cp_s = CPs_fine; %CPs_smooth;

%% Characteristic scaling terms

%Heat transport
ks = Ks_fine;
kbar   = phi*kf + (1-phi).*ks; %phi/tau*kf + (1-phi).*ks;
rhobar = phi*rho_f*cp_f + (1-phi).*rho_s.*cp_s;

%% Boundary conditions pressure

% Parameters: pressure
Param.p.dof_dir   = Grid.dof_ymax;%[Grid.dof_ymin;Grid.dof_ymax]; %[Grid.dof_ymax];
Param.p.dof_f_dir = Grid.dof_f_ymax;%[Grid.dof_f_ymin;Grid.dof_f_ymax]; %[Grid.dof_f_ymax];
Param.p.g         = rho_f*grav*(ocean_thick_vector);%[rho_f*grav*(ocean_thick_vector+Depth);rho_f*grav*(ocean_thick_vector)]; %[rho_f*grav*(ocean_thick_vector)]; 
Param.p.dof_neu   = [];
Param.p.dof_f_neu = [];
Param.p.qb        = [];

[Bp,Np,fn_p] = build_bnd(Param.p,Grid,I);

%% Initialize pressure and temperature

%pressure
Ps                 = flipud(reshape(rho_f*grav*(ocean_thick_matrix+Yc),Grid.N,1));

%fluxes
%Q                = zeros(Grid.Nf,t_totals);

%Heat transport operators
theta_a          = 1;                          %advection solved explicitly
theta_d          = 0;                          %diffusion solved implicitly

Ts = zeros(Ny,Nx);
for i = 1:Nx
    Ts(:,i) = linspace(T_bot(i),T_s,Ny);
end
Ts = reshape(Ts,[Nx*Ny 1]);

%temperature
T_i  = Ts;
T_ii = Ts;

%diagonalized matrices
Rhobar             = spdiags(rhobar(:),0,Grid.N,Grid.N);%%%%%%%%%%
Kbar               = comp_mean(kbar,1,kvkh,Grid); %%%%%%%%%%% 10 Aug 2017: arithmetic mean instead of harmonic mean
Rho_f              = rho_f;
KdMu               = Kd/Mu;

%Initialize loop variables
t = 0;
count = 1;

% mov = matfile('movie_dom_donny.mat','Writable',true);

%% Loop

while t < endtime
    
    t = t + dt;
    %% Solve steady state pressure --> Darcy Flux
    %pressure
    Lp             = -D*KdMu*G;
    fs_grav        = D*(KdMu*Rho_f*grav)*G*(ocean_thick_matrix(:)+Yc(:));
%     fs_grav        = D*(KdMu*Rho_f*grav)*G*(ocean_thick+Yc(:));
    fs_p           = zeros(Grid.N,1);
    fp             = fs_p + fn_p + fs_grav;
    P              = solve_lbvp(Lp,fp,Bp,Param.p.g,Np);
    
    q              = comp_flux_p(D,KdMu,Rho_f*grav*G*Yc(:),G,P,fs_p,Grid,Param.p);
    q(isnan(q)) = 0;
%     Q(:,count+1)   = q;
    
    qtop           = q(Grid.dof_f_ymax);
    qdown          = qtop(qtop<0);
    qup            = qtop(qtop>0);
    qdown_dof      = Grid.dof_ymax(qtop < 0);
    qdown_dof_f    = Grid.dof_f_ymax(qtop < 0);
    
    [qx,qy,qmags]  = center_flux(q,Grid);
    qy_up = ones(size(qy))*1e-20; qy_up(qy>1e-20) = qy(qy>1e-20);
    [psi,~,~]      = comp_streamfun(q,Grid);
    %STORE AT APPROPRIATE TIME STEP
%     if xor(t == times(count), t > times(count))
%         Qx(:,count)    = reshape(qx,Grid.N,1);
%         Qy(:,count)    = reshape(qy,Grid.N,1);
%         Qx_null(:,count) = q(1:Grid.Nfx);
%         Qy_null(:,count) = q(Grid.Nfx+1:end);
%     end
%     Psi(:,count)       = reshape(psi,numel(psi),1);
    % 2nd order upwind of fluid fluxes
    Aq = flux_upwind(q,Grid);
    
     %% Solve for groundwater age
        
%     %PDE
%     L.a = D*Aq;
%     fs.a = phi(:);
% 
%     %Boundry conditions     
%     Param.a.dof_dir   = Grid.dof_ymax(qtop<0);
%     Param.a.dof_f_dir = Grid.dof_f_ymax(qtop<0);
%     Param.a.g         = zeros(length(qdown),1)+0.00001;
%     Param.a.dof_neu   = [];
%     Param.a.dof_f_neu = [];
%     Param.a.qb        = [];
% 
%     [B.a, N.a, fn.a]  = build_bnd(Param.a,Grid,I);
% 
%     %solve for groundwater ages
%     ages = solve_lbvp(L.a,fs.a+fn.a,B.a,Param.a.g,N.a); %[s]
%     for i = 1:length(ages)
%         if ages(i) == Inf || ages(i) <= 0 || isnan(ages(i)) == 1
%             if i == 1
%                 ages(i) = ages(i+1);
%             else
%                 ages(i) = ages(i-1);
%             end
%         end
%     end
%     ages = reshape(ages,Grid.Ny,Grid.Nx);
    
    %% Boundary conditions temperature
    
    % Parameters: temperature
    Param.T.dof_dir   = [Grid.dof_ymin;qdown_dof]; %[Grid.dof_ymin;Grid.dof_ymax];
    Param.T.dof_f_dir = [Grid.dof_f_ymin;qdown_dof_f]; %[Grid.dof_f_ymin;Grid.dof_f_ymax];
    Param.T.g         = [T_bot;T_s_vec(qtop<0)];

    Param.T.dof_neu   = [];%Grid.dof_ymin;
    Param.T.dof_f_neu = [];%Grid.dof_f_ymin;
    Param.T.qb        = [];%qt_b;

    [BT,NT,fn_T] = build_bnd(Param.T,Grid,I);
    
    %% Solve for temperature
    % Build heat operators
    Im_T             = @(theta_a,theta_d) ...
                     Rhobar + dt*(1-theta_a)*D*(rho_f*cp_f*Aq) ...
                            - dt*(1-theta_d)*D*(Kbar*G); %%%%%%
    Ex_T             = @(theta_a,theta_d) ...
                     Rhobar - dt*theta_a*D*(rho_f*cp_f*Aq)...
                            + dt*theta_d*D*(Kbar*G);     %%%%%%
    T_new      = solve_lbvp(Im_T(theta_a,theta_d),...
                            Ex_T(theta_a,theta_d)*T_i + ...
                            dt*(fn_T),BT,Param.T.g,NT);
                        
    T_new(T_new<T_s) = T_s;
    T_new(T_new>max(T_bot)) = max(T_bot);                    
                        
    %store for next time step
    T_ii = T_i;
    T_i  = T_new;
    
%     [T_ave,T_top_ave,power_out] = deal(zeros([1 t_totals]));
    T_ave(count) = mean(T_new-273.15);
    T_new_rect = reshape(T_new-273.15,Grid.Ny,Grid.Nx);
    T_top = T_new_rect(end,:);
    T_top_qup = T_top(qtop > 0); T_top_ave(count) = mean(T_top_qup); 
    sum_qup(count) = sum(qtop(qtop > 0)*dx);
    power_out = 1000 * cp_f * 1000 * (T_top-2).*qtop'; power_out(power_out<0) = 0; %milliWatts per m^2
    mean_power_out(count) = sum(1000 * cp_f * 1000 * ((T_top_qup-2).*qtop(qtop > 0)') / (length(qtop))); %milliWatts per m^2
    
    T_grad_top = (T_new_rect(end-1,:) - T_new_rect(end,:)) / dy;
    cond_power_out = 1000 * T_grad_top(8:end-7) .* kbar(end,8:end-7); %[mW]
    cond_mean_power_out(count) = sum(cond_power_out)/length(cond_power_out); %[mW]
        
    %% Update T-dependent variables
    %update fluid density
    [new_rho,~,new_mu] = steam_table(T_new,P,Grid);
    Rho_f         = comp_mean(reshape(new_rho,Grid.Ny,Grid.Nx),1,kvkh,Grid);
    new_Mu        = comp_mean(reshape(new_mu,Grid.Ny,Grid.Nx),1,kvkh,Grid);
        
    KdMu = comp_mean(k./reshape(new_mu,Grid.Ny,Grid.Nx),1,kvkh,Grid);
    
    if xor(t == times(count),t > times(count))
        fprintf('\n%d of %d iterations\n',count,t_totals)
        count = count + 1;
    end
    
%% Plot while calculating
plotint = 50; %Plotting interval
if mod(count,plotint) == 1

figure(4)
    set(gcf,'units','normalized','outerposition',[0.8 0 0.2 0.5]) %[0 0 1 0.9]
    subplot(3,1,1)
    contourf(Xc,Yc,T_new_rect,100,'LineStyle','none');
%     daspect([3 1 1]);
    colormap(jet);
    xlabel('Distance [m]')
    ylabel('Elevation [m]')
    c = colorbar;
    caxis([T_s-273.15 100]) %max(T_bot_vec)-273.15
    title(c,'T')
%     grid on
%     axis equal
    title(sprintf('Temperature after %5d years',round(times(count-1)/(yrs2s))))
%     set(gca,'xtick',[],'ytick',[]);

    subplot(3,1,2)
    contourf(Xc,Yc,log10(qmags),100,'LineStyle','none'); colorbar; hold on;
%     daspect([3 1 1]);
    logvec = fliplr(logspace(-6,-16,6));
    [qmax,qmaxd] = max(qmags); qmaxd = qmaxd*dy;
    jetf = flipud(jet);
    colormap(gca,jetf)
    colorbar('YTick', log10(logvec), 'YTickLabel', log10(logvec))
    caxis([-16 -6])
    contour(Xf,Yf,psi,20,'LineColor','w','LineWidth',0.5);
    plot(dx:dx:Length,qmaxd,'w.','MarkerSize',8); hold off;
    xlabel('x (m)'); ylabel('z (m)');
    title('Flow magnitude');
%     set(gca,'xtick',[],'ytick',[]);

    subplot(3,1,3)
    contourf(Xc,Yc,log10(K_fine),100,'LineStyle','none'); %axis equal;
    colormap(gca,jetf); colorbar 
    caxis([-20 -10])
    logvec = fliplr(logspace(-10,-20,6));
    colorbar('YTick', log10(logvec), 'YTickLabel', log10(logvec))
    title('Permeability [m^2]')
    xlabel('x [m]'); ylabel('z [m]');
    
%     subplot(3,1,3)
%     contourf(Xc,Yc,log10(ages),100,'LineStyle','none');
%     logvec = fliplr(logspace(18,8,6));
%     colormap(jet);
%     colorbar('YTick', log10(logvec), 'YTickLabel', log10(logvec))
%     caxis([8 18])
%     xlabel('Distance [m]')
%     ylabel('Elevation [m]')
% %     caxis([T_s-273.15 max(T_bot)-273.15])
%     title('Fluid residence time [s]')
    
figure(5)
set(gcf,'units','normalized','outerposition',[0.7 0 0.1 0.5]) %[0 0 1 0.9]
    subplot(5,1,1)
    plot(T_top)
    title('T along top boundary');
    subplot(5,1,2)
    plot(qtop)
    title('Flux across the top boundary [m/s]');
%     plot(T_ave)
%     title('Average System Temperature [C]');
    subplot(5,1,3)
    plot(mean_power_out(2:end))
    title('Average advective power output [mW]');
    subplot(5,1,4)
    plot(cond_power_out)
    title('Conductive power output along top bondary [mW]')
    subplot(5,1,5)
    plot(cond_mean_power_out(2:end))
    title('Average conductive power output [mW]');
   
%     plot(power_out)
%     title('Fluid power output along top bondary [mW]')
%     plot(qtop)
%     title('Flux across the top boundary [m/s]');
%     plot(T_top_ave)
%     title('Average temperature of outflow [C]');
%     plot(sum_qup(2:end))
%     title('Cumulative outflow [m/s]')
 
    drawnow

%     clear F
%     cdata = print('-RGBImage','-r85');
%     F = im2frame(cdata);
%     mov.F(1,count/plotint) = F;

%     figure(2)
%     set(gcf,'units','normalized','outerposition',[0 0 0.5 0.9])
%     subplot(3,1,1)
%     contourf(Xc,Yc,reshape(T_new-273.15,Grid.Ny,Grid.Nx),20,'LineStyle','none');
%     colormap(jet);
% %     set(gcf,'units','normalized','outerposition',[0 0 1 1])
%     ylabel('Elevation [m]')
%     c = colorbar('location','southoutside');
%     caxis([T_s-273.15 T_bot-273.15])
%     ylabel(c,'T')
% %     axis equal
%     title(sprintf('Temperature after t = %.2f yr',times(count-1)/(yrs2s)))
%     
%     subplot(3,1,2)
%     plot(-ocean_thick_vector_orig); hold on;
%     plot(-ocean_thick_vector_smooth); hold off;
%     title('Seafloor');
%     
%     subplot(3,1,3)
%     if t > dt
%         cla(p2a);% cla(p2b);
%     end
%     qtops = q(Grid.dof_f_ymax-1); %grid cells beneath uppermost surface
%     %q at dof_f_ymax seems to be consistently negative, not sure why
%     p2a = area(Grid.xc,qtops,'FaceColor','b'); grid on; hold on;
%     qneg = qtops; qneg(qneg >= 0) = NaN;
%     p2b = area(Grid.xc,qneg,'Facecolor','r');
%     xlabel('Distance [m]')
%     ylabel('Upper boundary q [m/s]')
%     xlim([0 Length]);
%     title('Flow across upper model boundary');

%     cdata = print('-RGBImage','-r85');
%     F = im2frame(cdata);
% %     F = getframe(gcf);
%     mov.F(1,count-1) = F;
%     
%     drawnow

end
    
end

times = [0 times];

%% Plot streamlines

% figure(4)
% set(gcf,'units','normalized','outerposition',[0 0 1 1])
% 
% [PSI,~,~] = comp_streamfun(q,Grid);
% subplot(2,4,1)
% % contour(Xc,Yc,reshape(Qy(:,3),Ny,Nx),20,'LineWidth',1,'LineColor','r'); axis equal;
% % contour(Xc,Yc,qy,20,'LineWidth',1,'LineColor','r'); axis equal;
% plot_flownet(20,20,0,PSI,'b--','r-',Grid);
% title('Streamlines');
% subplot(2,4,2)
% contourf(Xc,Yc,qmags,'LineStyle','none'); colorbar; axis equal;
% xlabel('x (m)'); ylabel('z (m)');
% title('Flux magnitude');
% subplot(2,4,3)
% contourf(Xc,Yc,reshape(Qx(:,end),Ny,Nx),'LineStyle','none'); colorbar; axis equal;
% xlabel('x (m)'); ylabel('z (m)');
% title('Horizontal flux');
% subplot(2,4,4)
% contourf(Xc,Yc,reshape(Qy(:,end),Ny,Nx),'LineStyle','none'); colorbar; axis equal;
% xlabel('x (m)'); ylabel('z (m)');
% title('Vertical flux');
% subplot(2,4,5)
% contourf(Xc,Yc,reshape(Ps,Ny,Nx),'LineStyle','none'); colorbar; axis equal;
% xlabel('x (m)'); ylabel('z (m)');
% title('Updated pressure');
% subplot(2,4,6)
% contourf(Xc,Yc,reshape(new_rho,Ny,Nx),'LineStyle','none'); colorbar; axis equal;
% xlabel('x (m)'); ylabel('z (m)');
% title('Updated fluid density');
% subplot(2,4,7)
% contourf(Xc,Yc,reshape(new_mu,Ny,Nx),'LineStyle','none'); colorbar; axis equal;
% xlabel('x (m)'); ylabel('z (m)');
% title('Updated fluid viscosity');
% % subplot(3,3,7)
% % contourf(Xc,Yc,phi,'LineStyle','none'); colorbar; axis equal;
% % xlabel('x (m)'); ylabel('z (m)');
% % title('Porosity');
% subplot(2,4,8)
% contourf(Xc,Yc,k,'LineStyle','none'); colorbar; axis equal;
% xlabel('x (m)'); ylabel('z (m)');
% title('Permeability');
% % subplot(3,3,9)
% % contourf(Xc,Yc,kbar,'LineStyle','none'); colorbar; axis equal;
% % xlabel('x (m)'); ylabel('z (m)');
% % title('Thermal conductivity');

%% Calculate stats

% layer_upper = zeros(size(layer2A));
% for i = 1:size(layer2A,2)
%     for j = size(layer2A,1):-1:1
%         if layer2A(j,i) == 1
%             layer_upper(j,i) = 1;
%             break;
%         end
%     end
% end
% layer_upper = logical(layer_upper);

flux_sf = sum(qup)/length(qup);
sprintf('Average flux across the seafloor = %d [m/s]',flux_sf)
% flux_sf_stderr = std(qup)/sqrt(length(qup));
% sprintf('Standard error of FASF = %d [m/s]',flux_sf_stderr)
flux_sf_std = std(qup);
sprintf('Standard deviation of FASF = %d [m/s]',flux_sf_std)
% max_flow = max(qmax);
% sprintf('Maximum flow velocity = %d [m/s]', max_flow)
% med_age_2A = median(ages(layer2A)/yrs2s);
% ave_age_2A = mean(ages(layer2A)/yrs2s);
% sprintf('Fluid residence time in extrusive layer: Median = %d; Mean = %d [yrs]',med_age_2A, ave_age_2A)
% med_age_upper = median(ages(layer_upper)/yrs2s);
% ave_age_upper = mean(ages(layer_upper)/yrs2s);
% sprintf('Fluid residence time just below sediment: Median = %d; Mean = %d [yrs]',med_age_upper, ave_age_upper)
% ave_flow_upper = mean(qmags(layer_upper));
ave_flow_2A = mean(qmags(layer2A));
sprintf('Mean flow velocity in the upper crust = %d [m/s]', ave_flow_2A) %, ave_flow_upper)
% flow_2A_stderr = std(qmags(layer2A))/sqrt(length(qmags(layer2A)));
% sprintf('Standard error of F2A = %d [m/s]',flow_2A_stderr)
flow_2A_std = std(qmags(layer2A));
sprintf('Standard deviation of F2A = %d [m/s]',flow_2A_std)

save('T_ave_crest.mat','T_ave');

%% Plots export

% figure
% set(gcf,'units','normalized','outerposition',[0 0 1 1]) 
% contourf(Xc,Yc,T_new_rect,100,'LineStyle','none');
% daspect([2 1 1]);
% colormap(jet); colorbar
% caxis([T_s-273.15 200]) 
% set(gca,'xtick',[],'ytick',[]);
% % print('temp','-dpng','-r300')

% figure
% set(gcf,'units','normalized','outerposition',[0 0 1 1]) 
% contourf(Xc,Yc,log10(qmags),100,'LineStyle','none'); hold on;
% daspect([3 1 1]);
% logvec = fliplr(logspace(-6,-16,6));
% [qmax,qmaxd] = max(qmags); qmaxd = qmaxd*dy;
% jetf = flipud(jet);
% colormap(gca,jetf); colorbar
% caxis([-16 -7])
% contour(Xf,Yf,psi,20,'LineColor','w','LineWidth',1.5);
% % plot(dx:dx:Length,qmaxd,'w.','MarkerSize',8); hold off;
% set(gca,'xtick',[],'ytick',[]);
% % print('flowmag','-dpng','-r300')

% figure
% set(gcf,'units','normalized','outerposition',[0 0 1 1])
% contourf(por_fine,100,'LineStyle','none');
% daspect([3 1 1]);
% set(gca, 'colormap', jetf); %colorbar
% caxis([0 0.3])
% set(gca,'xtick',[],'ytick',[])
% % print('porosity','-dpng','-r300')

% figure
% set(gcf,'units','normalized','outerposition',[0 0 1 1])
% contourf(log10(K_fine),100,'LineStyle','none')
% daspect([3 1 1])
% colormap(gca,jetf); %colorbar 
% caxis([-20 -10])
% set(gca,'xtick',[],'ytick',[])
% % print('permeability','-dpng','-r300')

% figure
% set(gcf,'units','normalized','outerposition',[0 0 0.3 0.2])
% plot(power_out,'b','linewidth',1); hold on %'bo','markerfacecolor','b','markersize',3
% plot(cond_power_out,'r','linewidth',1.5) %'ro','markerfacecolor','r','markersize',3
% % set(gca,'yscale','log')
% xlim([0 701]); ylim([0 1500])
% legend('Conductive','Advective','location','best')
% set(gca,'xtick',[],'ytick',[])
% % print('heatflow','-dpsc')

%% Write out movie

% clear
% load movie_dom_donny.mat
% % F(1) = F(2);
% v = VideoWriter('convection_crest_aug05_dx_75_dy_37.5_sx_1_sy_1_1y_1m_Tbotvar.avi');
% open(v)
% writeVideo(v,F)
% close(v)

%%

% [XF,YF,ZF] = meshgrid(F(5).cdata(1,:,1),F(5).cdata(:,1,1),F(5).cdata(1,1,:));
% for i = 1 %1:length(F)
%     if size(F(i).cdata) ~= size(F(5).cdata)
% %         [Xi,Yi,Zi] = meshgrid(F(i).cdata(:,1,1),F(i).cdata(1,:,1),F(i).cdata(1,1,:));
%         F(i).cdata = interp3(F(i).cdata,XF,YF,ZF);
%     end
% end
