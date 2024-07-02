% This example demonstrates the use of k-Wave for the reconstruction of a
% two-dimensional photoacoustic wave-field recorded  over a linear array of
% sensor elements  The sensor data is simulated using kspaceFirstOrder2D
% and reconstructed using kspaceLineRecon. It builds on the Homogeneous
% Propagation Medium and Heterogeneous Propagation Medium examples. 

%% Clean up matlab
clear; close all; clc; rng shuffle;
addpath('./k-Wave')
%delete(gcp)
%parpool
%% Helper functions
normZero = @(x) x-max(x(:));
rf2Bmode = @(x) 20*log10(abs(hilbert(x)));

%% create the computational grid
elem_pitch=0.30e-3;
%PML_size = [100 350];       % size of the PML in grid points
PML_size = [168 64];       % size of the PML in grid points TONY
%PML_size = [20 20];       % size of the PML in grid points
%PML_size = [40 40];       % size of the PML in grid points

%Ny = 1024 - 2*PML_size(2);   % number of grid points in the y (column) direction
%dx = elem_pitch/5;               % grid point spacing in the x direction [m]
%dx = 0.04*1e-3;               % grid point spacing in the x direction [m]
%Ny = 256;   % number of grid points in the y (column) direction
%Ny = 512;
%Ny = 140;
Ny = 512/2;
ratio = 2;
dy = elem_pitch/ratio;          % grid point spacing in the y direction [m]
dx = dy/4;
%Nx = 0.050/dx;  % number of grid points in the x (row) direction
%Nx = 0.040/dx;  % number of grid points in the x (row) direction
Nx = 600*2;
%Ny = 140
%Nx = 10*floor(Nx/10);
%dy = 4*dx;          % grid point spacing in the y direction [m]

%f1max = (1500/(2*dx))/1e6;
%f2max = (1500/(2*dy))/1e6;

kgrid = makeGrid(Nx, dx, Ny, dy);

%% Define the medium properties
sd = 0.02;
medium.sound_speed = 1500-0*makeDisc(Nx,Ny,0,0,10e-3/dy);
medium.sound_speed = medium.sound_speed+1540*sd*randn(size(medium.sound_speed));
medium.density = 1000-0*makeDisc(Nx,Ny,0,0,10e-3/dy);
medium.density = medium.density + 1000*sd*randn(size(medium.density));

medium.sound_speed(medium.sound_speed>1600) = 1600;
medium.sound_speed(medium.sound_speed<1400) = 1400;

%medium.sound_speed(medium.sound_speed>1600) = 1600;
%medium.sound_speed(medium.sound_speed<1400) = 1400;

% d = Nx*dx;
% medium.sound_speed(round(Nx*0.010/d),round(Ny/2):round(Ny/2)) = 1570;
% medium.sound_speed(round(Nx*0.020/d),round(Ny/2):round(Ny/2)) = 1570;
% medium.sound_speed(round(Nx*0.030/d),round(Ny/2):round(Ny/2)) = 1570;
% 
% medium.density(round(Nx*0.010/d),round(Ny/2):round(Ny/2)) = 1020;
% medium.density(round(Nx*0.020/d),round(Ny/2):round(Ny/2)) = 1020;
% medium.density(round(Nx*0.030/d),round(Ny/2):round(Ny/2)) = 1020;

Nx_tot = Nx;
Ny_tot = Ny;
% Attenuation
rx = ones(Nx_tot,1)*linspace(-Ny_tot*dy/2,Ny_tot*dy/2,Ny_tot);
rz = linspace(0,Nx_tot*dx,Nx)'*ones(1,Ny_tot);

alpha_coeff_map_2D = 0.50*ones([Nx, Ny]);
figure(1000); imagesc(1e2*linspace(-Ny_tot*dy/2,Ny_tot*dy/2,Ny_tot),1e2*linspace(0,Nx_tot*dx,Nx),alpha_coeff_map_2D); colorbar; caxis([0 1.5]);
axis image

radius = sqrt((rx).^2+(rz-0.020).^2);
alpha_coeff_map_2D(radius <= 0.010) = 0.50;
figure(1001); imagesc(1e2*linspace(-Ny_tot*dy/2,Ny_tot*dy/2,Ny_tot),1e2*linspace(0,Nx_tot*dx,Nx),alpha_coeff_map_2D); colorbar; caxis([0 1.5]);
axis image

alpha_coeff_map = alpha_coeff_map_2D;

medium.alpha_coeff = alpha_coeff_map;
%medium.alpha_coeff = 0;
%medium.alpha_power = 1.05;
medium.alpha_power = 1.05;
%figure; imagesc(medium.alpha_coeff)
% Generally want high density sources rather than sound speed sources
% since sound speed effects time step.
% medium.density = medium.density+2000*(makeDisc(Nx,Ny,0.25*Nx,Ny/2,1)+makeDisc(Nx,Ny,0.75*Nx,Ny/2,1));
%Compute time array based on material properties

% create the time array
fs = 60*1e6;
c0 =1540;
t_end = (Nx*dx)*2.2/c0;     % [s]
%t_end = (Nx*dx)*2/c0;     % [s]
kgrid.t_array = 0:1/fs:t_end;
dt = 1/fs;
%[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed,0.40875,1.2*Nx*dx*2/1540);

focal_distance = 0.020;   % center of circle
focal_number = 2;
nAperture = (focal_distance/focal_number)/dy;
nAperture = ratio*floor(nAperture/ratio);
nApertureEle =  nAperture/ratio;

nLines = floor(Ny/ratio);
for ii = 1:nLines, 
  % ii= ceil(nLines/2)
jj = ratio*ii;
axis_x = rz(:,1);
axis_y = rx(1,:);
disp(['Lines: ',num2str(ii),' de ',num2str(nLines)]);
src_ini = max(1,jj-nAperture/2) ;
src_fin = min(Ny,jj+nAperture/2-1) ;

[temp,pos] = min(abs(axis_x-focal_distance));
focus_point = [axis_y(jj) axis_x(pos)];

aperture_point_src = [axis_y(src_ini:src_fin)' axis_x(5)*ones(src_fin-src_ini+1,1)];

figure (6); plot(aperture_point_src(:,1),aperture_point_src(:,2),'sb');
hold on; plot(focus_point(:,1),focus_point(:,1),'sr'); hold off;

%% Define Sensor properties
sensor.mask = zeros(Nx, Ny);
% Need a slight offset here to prevent backwards propagating pulse
sensor.mask(5, src_ini:ratio:src_fin) = 1;
%sensor.mask(5, 2:4:end) = 1;
%sensor.mask(5, 3:4:end) = 1;
%sensor.mask(5, 4:4:end) = 1;
sensor.directivity_size = 0.2698e-3;
sensor.directivity_angle = zeros(size(sensor.mask));


%% Define Source parameters
amp = 5; % [au]
%source.p_mask = zeros(size(sensor.mask));
%mask = source.p_mask(1,:);
%mask(src_ini:src_fin) = 1;
source.p_mask = zeros(Nx, Ny);
%source.p_mask(5,src_ini:src_fin) = sensor.mask(5,:).*mask;
source.p_mask(5,src_ini:src_fin) = 1;
%source.p_mask(5,src_ini:2:src_fin) = 1;
%figure; imagesc(source.p_mask);
%apod = hamming(nnz(source.p_mask));
apod = boxcar(nnz(source.p_mask));
%apod = boxcar(nnz(source.p_mask));
%apod = boxcar(src_fin-src_ini+1);
excitation = amp*toneBurst(1/dt,20/3*1e6,5);
src_exc = apod(:)*excitation(:)';
%%
% 
%  PREFORMATTED
%  TEXT
% 
angle = 0;
%source.p1 = steer_delay(src_exc,angle,dy,dt,1540);
%figure; imagesc(src_exc);
%figure; imagesc(source.p1);
source.p = steer_delay_focus(src_exc,angle,dy,dt,1540,aperture_point_src,focus_point);
%figure; imagesc(source.p2);
%source.p = source.p2;
%% set the input arguments:
% force the PML to be outside the computational grid; 
% switch off p0 smoothing within kspaceFirstOrder2D
PML_alpha = 2;
DATA_CAST       = 'gpuArray-single';     % set to 'single' or 'gpuArray-single' to speed up computations
input_args = {'PMLInside', false, 'PMLAlpha', PML_alpha, 'PMLSize', PML_size, 'PlotPML', false,...
    'Smooth', false,'PlotSim',false,'DataCast',DATA_CAST};

% run the simulation
colormap gray
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

sensor_data=sensor_data';

%% Resample to desired frequency
sample_freq = fs;
[resamp_n,resamp_d] = rat(sample_freq*kgrid.dt,1e-8);
ds_sensor = resample(sensor_data,resamp_n,resamp_d);

% Check resampled spectrum
% figure
% subplot(211)
% plot(((1:size(sensor_data,1))-1)*kgrid.dt,sensor_data(:,end/2),'.-');
% grid on; hold all
% plot(((1:size(ds_sensor,1))-1)/sample_freq,ds_sensor(:,end/2),'.-');
% subplot(212)
% quickSpectrum(sensor_data(:, end/2),1/kgrid.dt);
% hold all
% quickSpectrum(ds_sensor(:,end/2),sample_freq);
% xlim([0,sample_freq]*1e-6)
sensor_data = ds_sensor;
clear ds_sensor

%% Plot pre-BF data
figure(3)
cleanImagesc(normZero(rf2Bmode(sensor_data)),[-60 0])
title('PreBF data')

%% Beamform
%max_apert = 64; 
max_apert = 64;
f_number = 2;
c_bf = 1540;
bf_data = BFangle2(sensor_data,max_apert,sample_freq,c_bf,elem_pitch,'rect',f_number,0);
bf_data2 = BFangle2(sensor_data,max_apert,sample_freq,c_bf,elem_pitch,'bh',f_number,0);

if src_ini <= 1,
    index = size(bf_data,2) - 16;
elseif src_fin>= Ny,
    index = 17;
else
    index = 17;
end
  sensor_data_total{ii} = sensor_data;
  bf_data_total{ii} = bf_data;
  bf_data_total_bh{ii} = bf_data2;
  bf_data_final(:,ii) = bf_data(:,index);
  
  bbb=4;
% index = 
%bf_data_final(:,1 + (jj-1)/4 ) = bf_data(:,jj);
bmode = abs(hilbert(bf_data));
bmode = 20*log10(bmode);
bmode = bmode / max(bmode(:));
figure(10); imagesc(bmode); caxis([-40 0]); colormap gray; colorbar
aaabb = 4;
% temp = bf_data(50:end,:);
% figure; imagesc(abs(temp));
% figure; plot(bf_data(:,round(end/2)));
end
save bf_data_final_homV17_2cm_inc.mat

axAxis = 0:size(bf_data,1)-1; axAxis = axAxis*1/sample_freq*c_bf/2;
latAxis = 0:size(bf_data,2)-1; latAxis = latAxis-mean(latAxis); latAxis = latAxis *elem_pitch;

bmode = abs(hilbert(bf_data_final(200:end,:)));

bmode = 20*log10(bmode);
bmode = bmode - max(bmode(:));
figure; imagesc(latAxis*100,100*axAxis(200:end),bmode); colorbar; colormap gray; caxis([-40 0])
figure; imagesc(bmode); colorbar; colormap gray; caxis([-40 0])

%keyboard
% 
% sam1 = bf_data_final(200:end,:);
% sam1p = sam1(747:975,:);
% sam1p = sam1(221:449,:);
% sam1p = sam1(1273:1501,:);
% 
% bmode = abs(hilbert(sam1p));
% bmode = 20*log10(bmode);
% bmode = bmode - max(bmode(:));
% figure; imagesc(bmode); colorbar; colormap gray; caxis([-40 0])
% 
% xxx = xcorr2(sam1p,sam1p);
% aaa = xxx(round(end/2),:);
% aaa = aaa/max(aaa);
% figure; plot(interp(aaa,8));

bf_data = bf_data_final;
%aaa = 4;
%bf_data_final

%% Plot BF data
axAxis = 0:size(bf_data,1)-1; axAxis = axAxis*1/sample_freq*c_bf/2;
latAxis = 0:size(bf_data,2)-1; latAxis = latAxis-mean(latAxis); latAxis = latAxis *elem_pitch;
%latAxis = 0:size(bf_data,2)-1; latAxis = latAxis-mean(latAxis); latAxis = latAxis *2*dy;

figure
cleanImagesc(latAxis*1e3, axAxis*1e3, normZero(rf2Bmode(bf_data)),[-60 0])
title('BF data')
xlabel('Lat. Dist. [mm]')
ylabel('Axial Dist. [mm]')
axis image


rf_data = bf_data(50:end,:);
z = axAxis(50:end);
x = latAxis;
%x = latAxis(PML_size(2)+1:end-PML_size(2));
figure(99)
cleanImagesc(x*1e3, z*1e3, normZero(rf2Bmode(rf_data)),[-60 0])
title('BF data')
xlabel('Lat. Dist. [mm]')
ylabel('Axial Dist. [mm]')
axis image

save tony_homV17_2cm_inc.mat rf_data x z
return

%figure; plot(bf_data(:,64))