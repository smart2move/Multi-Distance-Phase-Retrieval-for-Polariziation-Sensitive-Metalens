clear all;
%% load in data
load("E2_measurements.mat") % Change to right filename

edit = 1; % Do you want to edit the measured data?
show = 0; % Do you want to see the measurement data?
params.HNcheck = 1; % Do you want to check for high-noise measurements?
phase_error_check = 0; % If data is simulated, check phase from algorithm against phase from simulation
perfect_phase_profile = 1; % Do you want a perfect phase profile to compare the algorithm against?
Phase_MDPR = 0; %Phase_esMDPR; % Choose which phase profile you want to compare to

%% set up variables
measured_intensities = E2_esMDPR; % Choose Es or Ep field

% measurements positions
[Nx, Ny, Nz] = size(measured_intensities);



% other variables
params.pixel_size = xspace/(Nx-1);
params.x = xf;
params.y = yf;
params.z = zf;
params.wavelength = 587e-9;
params.max_iterations = 1000;
params.convergence_threshold = 0.01;
params.kalman_start_weight = 1;
params.kalman_end_weight = 0;
params.feedback_weight1 = 0.6;
params.feedback_weight2 = 0.4;
params.intensity_threshold = 1023;
params.noise_level = 1;
params.zeropadding = 0; % in number of pixels, not distance
params.starting_phase = 0; % 0 for zero, 1 for perfect phase, 2 for random
params.lensradius = 171e-6;
params.focal_length = 10e-3;

%% Edit data for testing
% Remove measurements or shorten the data
if edit == 1
    Nm = 17 ;
    temp_E2 = zeros(Nx, Ny, Nm);
    z_temp =zeros(Nm, 1);

    for k = 1:Nm
        temp_E2(:, :, k) = measured_intensities(:, :, k+1);
        z_temp(k) = zf(k+1);
    end

    measured_intensities = temp_E2;
    params.z = z_temp;
    Nz = Nm;
end

%% Show intensity measurements
if show == 1
    for i = 1:Nz
        figure 
        imagesc(xf*1e6, yf*1e6, measured_intensities(:,:,i))
    end
end

%% Run algorithm
[final_phase, reconstructed_field, RMS_error] = multidistance_phase_retrieval_2(measured_intensities, params);

%% Check for high-noise measurements
if params.HNcheck == 1
    figure
    imagesc(RMS_error);
    title("RMS Error");
    xlabel("Intensity Measurement N");
    ylabel("Intensity Measurement N");
    ylabel(colorbar, 'RMS error', 'FontSize', 12)
end

%% plot phase
figure
tiledlayout(2, 2);

[X, Y] = meshgrid(xf, yf);

% Arperature setup
arperature = sqrt(X.^2 + Y.^2) < params.lensradius;
arp_phase = zeros(Nx, Ny);
arp_phase(arperature) = final_phase(arperature);

xi = (params.zeropadding/2 + 1):(Nx - params.zeropadding/2);
yi = (params.zeropadding/2 + 1):(Ny - params.zeropadding/2);


% phase plot from algorithm
nexttile;
imagesc(xf(xi) * 1e6, yf(yi) * 1e6, final_phase(xi, yi));
axis image;
title('Phasemap Algorithm');
ylabel(colorbar, 'Phase', 'FontSize', 12)
colormap('jet'); % You can change to 'hot', 'cool', 'parula', etc.
xlabel('x (um)');
ylabel('y (um)');
set(gca, 'YDir', 'normal');

% phase plot from algorithm with arperature
nexttile;
imagesc(xf(xi) * 1e6, yf(yi) * 1e6, arp_phase(xi, yi));
axis image;
title('Phasemap Algorithm');
ylabel(colorbar, 'Phase', 'FontSize', 12)
colormap('jet'); % You can change to 'hot', 'cool', 'parula', etc.
xlabel('x (um)');
ylabel('y (um)');
set(gca, 'YDir', 'normal');
        
% intensity squared of algorithm
nexttile;
imagesc(xf(xi) * 1e6, yf(yi) * 1e6, abs(reconstructed_field(xi, yi)));
axis image;
ylabel(colorbar, 'Intensity', 'FontSize', 12)
title('Intensity');
colormap('jet'); % You can change to 'hot', 'cool', 'parula', etc.
xlabel('x (um)');
ylabel('y (um)');
set(gca, 'YDir', 'normal');

%phase plot from measurement
nexttile;
field_diff = abs(abs(reconstructed_field) - sqrt(measured_intensities(:,:,1)));
imagesc(xf(xi) * 1e6, yf(yi) * 1e6, field_diff(xi, yi));
axis image;
ylabel(colorbar, 'Intensity', 'FontSize', 12)
title('Intensity Difference With Measurement');
colormap('jet'); % You can change to 'hot', 'cool', 'parula', etc.
xlabel('x (um)');
ylabel('y (um)');
set(gca, 'YDir', 'normal');

%% Hyperbolic phase profile for comparison
if(perfect_phase_profile == 1)
    phase_hyp = zeros(Nx, Ny);
    
    for i = 1:Nx
        for j = 1:Ny
            radius_hyp = sqrt(xf(i)^2 + yf(j)^2);
            if radius_hyp < params.lensradius
                phase_hyp(i, j) = (2*pi / params.wavelength) * (params.focal_length - sqrt(params.focal_length^2  + xf(i)^2 + yf(j)^2));
                while phase_hyp(i, j) < -1*pi
                    phase_hyp(i, j) = 2*pi + phase_hyp(i, j);
                end
            end
        end
    end

    final_phase_centred = final_phase - final_phase(floor(Nx / 2), floor(Ny / 2)) - phase_hyp(floor(Nx / 2), floor(Ny / 2));

    for i = 1:Nx
        for j = 1:Ny
            if(final_phase_centred(i,j) > pi)
                final_phase_centred(i,j) = final_phase_centred(i,j) - 2*pi;
            end
            if(final_phase_centred(i,j) < -1*pi)
                final_phase_centred(i,j) = final_phase_centred(i,j) + 2*pi;
            end
        end
    end
    
     Phase_esMDPR = phase_hyp;

     arp_phase_hyp = zeros(Nx,Ny,1);
     arp_phase_hyp(arperature) = phase_hyp(arperature);

     
   
    figure
    
    nexttile;
    imagesc(xf * 1e6, yf * 1e6, arp_phase_hyp);
    axis image;
    ylabel(colorbar, 'Phase error', 'FontSize', 12)
    title('Perfect phase profile for comparison');
    colormap('jet'); % You can change to 'hot', 'cool', 'parula', etc.
    xlabel('x (um)');
    ylabel('y (um)');
    set(gca, 'YDir', 'normal');

    % phase plot from algorithm with arperature
    nexttile;
    imagesc(xf(xi) * 1e6, yf(yi) * 1e6, final_phase_centred(xi, yi));
    axis image;
    title('Phasemap Algorithm');
    ylabel(colorbar, 'Phase', 'FontSize', 12)
    colormap('jet'); % You can change to 'hot', 'cool', 'parula', etc.
    xlabel('x (um)');
    ylabel('y (um)');
    set(gca, 'YDir', 'normal');
end
%% Phase error check
 if(phase_error_check == 1)
 % the loops are to search for the closest allignment. the lowest error is
 % the best alligned version. it searches around the centre, because its
 % most likely to be the same there.
 zca = 41; % radius of area to check for best zeroing match in pixels
 phase_error_degree = zeros(zca, zca);
 phase_error_degree = phase_error_degree + 360;
    for q = -((zca - 1)/2 - 1):((zca - 1)/2 + 1)
         for r = -((zca - 1)/2 - 1):((zca - 1)/2 + 1)
            % Centering the final phase so the phases are closer to eachother
            final_phase_centred = final_phase - (final_phase(floor(Nx / 2) + q, floor(Ny / 2) + r) - Phase_MDPR(floor(Nx / 2) + q, floor(Ny / 2) + r));
            for i = 1:Nx
                for j = 1:Ny
                    if(final_phase_centred(i,j) > pi)
                        final_phase_centred(i,j) = final_phase_centred(i,j) - 2*pi;
                    end
                    if(final_phase_centred(i,j) < -1*pi)
                        final_phase_centred(i,j) = final_phase_centred(i,j) + 2*pi;
                    end
                end
            end
            
            % Determine phase error
            phase_error = final_phase_centred - Phase_MDPR;
        
            % Use arperature to only capture the field in the lens
            arp_phase_error = zeros(Nx, Ny);
            arp_phase_error(arperature) = phase_error(arperature);
            abs_arp_phase_error = abs(arp_phase_error);
        
            % Adjust for -2pi - 2pi phase difference
            for i = 1:Nx
                for j = 1:Ny
                    if abs_arp_phase_error(i,j) > pi
                        abs_arp_phase_error(i,j) = abs(abs_arp_phase_error(i,j) - 2*pi);
                    end
                end
            end
        
            
            % Calculate mean error from final phase errors
            phase_error_sum = sum(abs_arp_phase_error(:));
            phase_error_mean = phase_error_sum / nnz(arperature);
            phase_error_degree((q + (zca - 1)/2), (r + (zca - 1)/2)) = (phase_error_mean * 360) / (2*pi);

            if (phase_error_degree((q + (zca - 1)/2), (r + (zca - 1)/2)) == min(phase_error_degree(:)))
                final_phase_centred_m = final_phase_centred;
                abs_arp_phase_error_m = abs_arp_phase_error;
            end


         end
    end

    phase_error_degree_final = min(phase_error_degree(:));



    fprintf('Phase error in degrees: %f\n', phase_error_degree_final);
    
    figure
    tiledlayout(1, 3);

    % phase plot from algorithm
    nexttile;
    imagesc(xf(xi) * 1e6, yf(yi) * 1e6, final_phase_centred_m(xi, yi));
    axis image;
    title('Phasemap Algorithm');
    ylabel(colorbar, 'Phase', 'FontSize', 12)
    colormap('jet'); % You can change to 'hot', 'cool', 'parula', etc.
    xlabel('x (um)');
    ylabel('y (um)');
    set(gca, 'YDir', 'normal');
    
    % phase from simulation
    nexttile;
    imagesc(xf(xi) * 1e6, yf(yi) * 1e6, Phase_esMDPR(xi, yi));
    axis image;
    title('Phasemap simulation');
    ylabel(colorbar, 'Phase', 'FontSize', 12)
    colormap('jet'); % You can change to 'hot', 'cool', 'parula', etc.
    xlabel('x (um)');
    ylabel('y (um)');
    set(gca, 'YDir', 'normal');
            
    adjusted phase
    nexttile;
    imagesc(xf(xi) * 1e6, yf(yi) * 1e6, abs_arp_phase_error_m(xi, yi));
    axis image;
    ylabel(colorbar, 'Phase error', 'FontSize', 12)
    title('Phase error');
    colormap('jet'); % You can change to 'hot', 'cool', 'parula', etc.
    xlabel('x (um)');
    ylabel('y (um)');
    set(gca, 'YDir', 'normal');
end
