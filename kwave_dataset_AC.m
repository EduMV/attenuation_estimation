clear; close all; clc; rng shuffle;

%gpuInfo = gpuDevice;
%disp(gpuInfo);
rng("default")

%% Iter
for sss = 301:302
% This example demonstrates the use of k-Wave for the reconstruction of a
% two-dimensional photoacoustic wave-field recorded  over a linear array of
% sensor elements  The sensor data is simulated using kspaceFirstOrder2D
% and reconstructed using kspaceLineRecon. It builds on the Homogeneous
% Propagation Medium and Heterogeneous Propagation Medium examples.
    clearvars -except sss; clc; rng shuffle;
    addpath(genpath(pwd))
 	addpath(genpath('/opt/MATLAB Add-Ons'))
    
		% DATA_CAST = 'single';     % set to 'single' or 'gpuArray-single' to speed up computations
    DATA_CAST = 'gpuArray-single';     % set to 'single' or 'gpuArray-single' to speed up computations
    delete(gcp)
    parpool
    %%
    normZero = @(x) x-max(x(:));
    rf2Bmode = @(x) 20*log10(abs(hilbert(x)));
    %%
    elem_pitch = 0.30e-3;
    
    pml_x_size = 2*40;                % [grid points]
    pml_y_size = 2*14;                % [grid points]
    
    % set total number of grid points not including the PML
    Nx = 2*600 - 2 * pml_x_size;      % [grid points]
    Ny = 2*540 - 2 * pml_y_size;      % [grid points]
    
    PML_size = [pml_x_size pml_y_size];       % size of the PML in grid points TONY
    
    ratio = 8;
    dy = elem_pitch/ratio;          % grid point spacing in the y direction [m]
    dx = dy;
    
    %kgrid = makeGrid(Nx, dx, Ny, dy);
    kgrid = kWaveGrid(Nx, dx, Ny, dy);
    
    Nx_tot = Nx;
    Ny_tot = Ny;
    rx = ones(Nx_tot,1)*linspace(-Ny_tot*dy/2,Ny_tot*dy/2,Ny_tot);
    rz = linspace(0,Nx_tot*dx,Nx)'*ones(1,Ny_tot);
    %% Settings
    c0 = 1540;
    t_end = (Nx*dx)*2/c0;     % [s]
    kgrid.makeTime(c0, [], t_end);
    fs = 1/kgrid.dt;
    %[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed,0.40875,1.2*Nx*dx*2/1540);
    
    focal_distance = 20e-3;   % center of circle
    focal_number = 2;
    nAperture = (focal_distance/focal_number)/dy;
    nAperture = ratio*floor(nAperture/ratio);
    nApertureEle =  nAperture/ratio;
    
    nAperture = nApertureEle*ratio;
    
    nLines = floor(Ny/ratio); % vary slightly plm_y to get nLines=128
    
    bf_data_final = nan(kgrid.Nt ,nLines);
    
    %% Parameters
    sos_mean = 1540;
    sos_std = 0;

    rho_mean = [1000 1000]; %randi([18 25],1,2)*100/2; %[1000 1000 1000 1000 1000];
    rho_std = [0.02 0.02]; %randi([2 10],1,2)/200; %[0.01 0.06];

%     alpha_mean = round([(randi([12 24])/20) (randi([12 24])/20)*(randi([10 25])/20)],2); %[0.3 1.5];
%     while abs(alpha_mean(2) - alpha_mean(1)) < 0.2
%         alpha_mean = round([(randi([12 24])/20) (randi([12 24])/20)*(randi([10 25])/20)],2);
%     end
    
    % Inicialización de alpha_mean con un rango y comprobación de la diferencia
    alpha_mean = round([0.3 + (randi([0 5])/10) 1.3 + (randi([0 5])/10)], 2); %[0.3 0.8] y [1.3 1.8]
    while abs(alpha_mean(2) - alpha_mean(1)) < 0.2
        alpha_mean = round([0.3 + (randi([0 5])/10) 1.3 + (randi([0 5])/10)], 2);
    end

    alpha_std = [0 0];

    offset = 5;
    
    % Inclusion
    radius_meters = (randi([8 10])/10)*(randi([16 24],1,2)/20)*1e-2; %[0.6-1.2]
    center_meters = [randi([0 38])/10 randi([16 24])/10]*1e-2;

    radius = radius_meters/dx;
    center = center_meters/dx + [0 offset];
    
    [xGrid, yGrid] = meshgrid(1:Ny, 1:Nx);
    
    distFromCenter = sqrt((xGrid - center(1)).^2 / radius(1)^2 + (yGrid - center(2)).^2 / radius(2)^2);
    circleMatrix = double(distFromCenter <= 1);
    %%
    % SoS
    medium.sound_speed = patterns(sos_mean, sos_std, 'homo', [], Nx, Ny);

    % Density
    density_mean = zeros(size(circleMatrix));
    density_mean(circleMatrix == 0) = rho_mean(1);
    density_mean(circleMatrix == 1) = rho_mean(2);
    
    density_std = zeros(size(circleMatrix));
    density_std(circleMatrix == 0) = rho_std(1);
    density_std(circleMatrix == 1) = rho_std(2);
    
    medium.density = patterns(density_mean, density_std, 'homo', [], Nx, Ny);
    
    % Attenuation
    alpha_coeff_mean = zeros(size(circleMatrix));
    alpha_coeff_mean(circleMatrix == 0) = alpha_mean(1);
    alpha_coeff_mean(circleMatrix == 1) = alpha_mean(2);
    
    alpha_coeff_std = zeros(size(circleMatrix));
    alpha_coeff_std(circleMatrix == 0) = alpha_std(1);
    alpha_coeff_std(circleMatrix == 1) = alpha_std(2);
    
    medium.alpha_coeff = patterns(alpha_coeff_mean, alpha_coeff_std, 'homo', [], Nx, Ny);
    
    medium.alpha_power = 1;
    medium.alpha_mode = 'no_dispersion';
    
    source_strength = 1e6;
    tone_burst_freq = 6e6;        % [Hz]
    tone_burst_cycles = 3.5;
    
    rho0 = 1000;
    
    for ii = 1:nLines
    
        jj = ratio*ii;
        axis_x = rz(:,1);
        axis_y = rx(1,:);
        disp(['Lines: ',num2str(ii),' de ',num2str(nLines)]);
        src_ini = max(1,jj - nAperture/2) ;
        src_fin = min(Ny,jj + nAperture/2-1) ;
    
        [temp,pos] = min(abs(axis_x-focal_distance));
        focus_point = [axis_y(jj) axis_x(pos)];
    
        aperture_point_src = [axis_y(src_ini:src_fin)' axis_x(offset)*ones(src_fin-src_ini+1,1)];
    
        figure (6); plot(aperture_point_src(:,1),aperture_point_src(:,2),'sb');
        hold on; plot(focus_point(:,1),focus_point(:,1),'sr'); hold off;
    
    %%
        sensor.mask = zeros(Nx, Ny);
        % Need a slight offset here to prevent backwards propagating pulse
        sensor.mask(offset, src_ini:ratio:src_fin) = 1;
        %sensor.directivity_size = 0.2698e-3;
        %sensor.directivity_angle = zeros(size(sensor.mask));
        %sensor.directivity_size = 0.4e-3;
        sensor.directivity_size = 10*kgrid.dx;
        sensor.directivity_angle = zeros(size(sensor.mask));
    %%
     %amp = 100000; % [au]
        source.p_mask = zeros(Nx, Ny);
        source.p_mask(offset,src_ini:src_fin) = 1;
        apod = boxcar(nnz(source.p_mask));
    
        input_signal_norm = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);
        input_signal = (source_strength ./ (c0 * rho0)) .* input_signal_norm;
    
        %excitation = amp*toneBurst(1/kgrid.dt,6.66*1e6,5);
        %src_exc = apod(:)*excitation(:)';
        src_exc = apod(:)*input_signal(:)';
    
        angle = 0;
        %source.p1 = steer_delay(src_exc,angle,dy,dt,1540);
        %figure; imagesc(src_exc);
        %figure; imagesc(source.p1);
        %source = steer_delay_focus(src_exc,angle,dy,kgrid.dt,1540,aperture_point_src,focus_point);
        %source3 = steer_delay_focus(src_exc,angle,dy,kgrid.dt,1540,aperture_point_src,[0 Inf]);
        source.p = src_exc;
        %figure; imagesc(source.p2);
        %source.p = source.p2;
    
        medium.sound_speed_ref = 1540;
        
        
        %%
            %PML_alpha = 2;   % Default is 2
        %DATA_CAST       = 'single';     % set to 'single' or 'gpuArray-single' to speed up computations
        %input_args = {'PMLInside', false, 'PMLAlpha', PML_alpha, 'PMLSize', PML_size, 'PlotPML', false,...
        %    'Smooth', false,'PlotSim',false, 'DataCast',DATA_CAST, 'DataRecast', true};
        input_args = {'PMLInside', false, 'PMLSize', PML_size, ... %'PMLAlpha', PML_alpha, ...
            'DataCast',DATA_CAST, 'DataRecast', true, 'PlotSim',false};
    
        % run the simulation
        %colormap gray
        sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
    
        sensor_data=sensor_data';
        %%
        
          max_apert = 64; % elements
        f_number = 2;
        c_bf = 1540;
        bf_data = BFangle(sensor_data,max_apert,fs,c_bf,elem_pitch,'rect',f_number,0);
    
        if src_ini <= 1
            index = size(bf_data,2) - floor(nApertureEle/2);
        elseif src_fin>= Ny
            index = floor(nApertureEle/2)+1;
        else
            index = floor(nApertureEle/2)+1;
        end
    
        bf_data_final(:,ii) = bf_data(:,index);
    
    end
        
    %%
    axAxis = 0:size(bf_data_final,1)-1; axAxis = axAxis*1/fs*c0/2;
    latAxis = 0:size(bf_data_final,2)-1; latAxis = latAxis-mean(latAxis); latAxis = latAxis *elem_pitch;
    
    rf = bf_data_final; %(100:end,:);
    z = axAxis; %(100:end);
    x = latAxis;
    figure; imagesc(x*1e3, z*1e3, normZero(rf2Bmode(rf)),[-60 0])
    title('B-mode'); colormap gray
    xlabel('Lateral distance (mm)')
    ylabel('Depth (mm)')
    axis image

    save("/mnt/nfs/rmarin/attenuation_estimation/out/rf_qus_AC_INCrand_"+num2str(sss)+".mat", "rf", "x", "z", "fs", "medium", "circleMatrix", ...
        "sos_mean", "sos_std", "rho_mean", "rho_std", "alpha_mean", "alpha_std","tone_burst_freq", ...
        "focal_distance", "radius_meters", "center_meters")
%     save("rf_qus_AC_INCrand_ref.mat", "rf", "x", "z", "fs", "medium", "circleMatrix", ...
%         "sos_mean", "sos_std", "rho_mean", "rho_std", "alpha_mean", "alpha_std","tone_burst_freq",...
%         "focal_distance")
    %return
end
function [med, med_mean, med_std] = patterns(var_mean, var_std, patt_type,var_params, Nx, Ny)
    switch patt_type
        case 'homo'
            med_mean = var_mean;
            med_std = var_std;
        case 'circle'
            circle_ind = boolean(makeDisc(Nx, Ny,...
                round(var_params.center_depth/var_params.dx)+var_params.offset, Ny/2, ...
                round(var_params.radius_disk/var_params.dx)));
            med_mean = var_mean(1)*(~circle_ind) + var_mean(2)*(circle_ind);
            med_std = var_std(1)*(~circle_ind) + var_std(2)*(circle_ind);
        case 'layers_vert'
            layer_pos = var_params.layer_pos/var_params.dx;
            [X,~] = meshgrid(1:Ny, 1:Nx);
            med_mean = var_mean(1)*(X<layer_pos) ...
                + var_mean(2)*(X>=layer_pos);
            med_std = var_std(1)*(X<layer_pos) ...
                + var_std(2)*(X>=layer_pos);
        case 'layers_horz'
            layer_pos = var_params.layer_pos/var_params.dy;
            [~,Y] = meshgrid(1:Ny, 1:Nx);
            med_mean = var_mean(1)*(Y<layer_pos) ...
                + var_mean(2)*(Y>=layer_pos);
            
            med_std = var_std(1)*(Y<layer_pos) ...
                + var_std(2)*(Y>=layer_pos);
        case 'layers_vert3'
            layer_pos = var_params.layer_pos/var_params.dx;
            [X,~] = meshgrid(1:Ny, 1:Nx);
            med_mean = var_mean(1)*(X<layer_pos(1)) ...
                + var_mean(2)*((X<layer_pos(2))&(X>=layer_pos(1)))...
                + var_mean(3)*(X>=layer_pos(2));
            med_std = var_std(1)*(X<layer_pos(1)) ...
                + var_std(2)*((X<layer_pos(2))&(X>=layer_pos(1)))...
                + var_std(3)*(X>=layer_pos(2));
        case 'layers_horz3'
            layer_pos = var_params.layer_pos/var_params.dy;
            [~,Y] = meshgrid(1:Ny, 1:Nx);
            med_mean = var_mean(1)*(Y<layer_pos(1)) ...
                + var_mean(2)*((Y<layer_pos(2))&(Y>=layer_pos(1)))...
                + var_mean(3)*(Y>=layer_pos(2));
            med_std = var_std(1)*(Y<layer_pos(1)) ...
                + var_std(2)*((Y<layer_pos(2))&(Y>=layer_pos(1)))...
                + var_std(3)*(Y>=layer_pos(2));
        case 'circleN'
            circle_ind = nan(Nx, Ny, length(var_params.center_depth));
            bg_mask = ones(Nx, Ny);

            med_mean = zeros(Nx, Ny);
            med_std = zeros(Nx, Ny);
            for cc = 1:length(var_params.center_depth)
                circle_ind(:,:,cc) = boolean(makeDisc(Nx, Ny,...
                    round(var_params.center_depth(cc)/var_params.dx)+var_params.offset, var_params.center_lat(cc), ...
                    round(var_params.radius_disk(cc)/var_params.dx)));
                med_mean = med_mean + var_mean(cc)*(circle_ind(:,:,cc));
                med_std = med_std + var_std(cc)*(circle_ind(:,:,cc));
                bg_mask = bg_mask & not(circle_ind(:,:,cc));
            end
            med_mean = med_mean + var_mean(end)*(bg_mask);
            med_std = med_std + var_std(end)*(bg_mask);
        
        otherwise
            med_mean = 0;
            med_std = 0;
            warning('Unexpected type')
    end
    med = med_mean + med_mean.*med_std.*randn(Nx, Ny);
end