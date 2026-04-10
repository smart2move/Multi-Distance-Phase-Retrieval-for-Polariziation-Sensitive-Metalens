function [final_phase, reconstructed_field, RMS_error] = multidistance_phase_retrieval_2(measured_intensities, params)
    % Implementation of the multi-distance phase retrieval algorithm
    % 
    % Inputs:
    %   measured_intensities - 3D array (Nx x Ny x Nz) of measured intensities at different z positions
    %   params - structure containing algorithm parameters:
    %       .pixel_size - CCD pixel size in meters
    %       .x - array with pixel locations in x-dimension
    %       .y - array with pixel locations in y-dimension
    %       .z - array of z distances corresponding to measurements
    %       .wavelength - wavelength in meters
    %       .max_iterations - maximum number of iterations
    %       .convergence_threshold - convergence criterion
    %       .kalman_start_weight - initial Kalman filter weight (typically 1)
    %       .kalman_end_weight - final Kalman filter weight (typically 0)
    %       .feedback_weight1 - weight for previous iteration in double feedback
    %       .feedback_weight2 - weight for current iteration in double feedback
    %       .intensity_threshold - for mask function (reject overexposed pixels)
    %       .noise_level - CCD noise level for mask function
    %       .zeropadding - total zero padding added in one dimension
    %       .starting_phase - starting phase of zero or random
    %       .lensradius - known radius of lens for arpature
    %       .HNcheck - 1 for High Noise check, 0 for regular MDPR
    %       .focal_length - focal lenght of the lens for comparison phase
    %    
    %
    % Outputs:
    %   final_phase - retrieved phase at sample plane
    %   reconstructed_field - complex field at sample plane
    %   RMS_error - error between different measurements backpropagated

    %% Initialization
    [Nx, Ny, Nz] = size(measured_intensities);
    k = 2*pi / params.wavelength;
    converg = 1;

    % Setting plot window in case of zero padding of data
    xi = (params.zeropadding/2 + 1):(Nx - params.zeropadding/2);
    yi = (params.zeropadding/2 + 1):(Ny - params.zeropadding/2);

    % Calculate spatial frequency vectors
    dx = params.pixel_size; dy = params.pixel_size;
    kx = 2*pi * (params.x / dx) / (Nx*dx);
    ky = 2*pi * (params.y / dy) / (Ny*dy);
    [Kx, Ky] = ndgrid(kx, ky);

    % Making sure the propegation distances start at zero
    if params.z(1) ~= 0
        params.z = params.z - params.z(1);
    end

    
    % Initialize field with random phase and sqrt of first intensity
    if params.starting_phase == 0
        phase_start = zeros(Nx, Ny);
    elseif params.starting_phase == 1

        phase_hyp = zeros(Nx, Ny);
        for i = 1:Nx
            for j = 1:Ny
            radius_hyp = sqrt(params.x(i)^2 + params.y(j)^2);
                if radius_hyp < params.lensradius
                    phase_hyp(i, j) = (2*pi / params.wavelength) * (params.focal_length - sqrt(params.focal_length^2  + params.x(i)^2 + params.y(j)^2));
                        while phase_hyp(i, j) < -1*pi
                            phase_hyp(i, j) = 2*pi + phase_hyp(i, j);
                        end
                end
            end
        end
        phase_start = phase_hyp/(2*pi);
       
    else
        phase_start = rand(Nx, Ny);
    end

    field_next = sqrt(measured_intensities(:,:,1)) .* exp(1i * 2*pi*phase_start);

    % Allocating space for Double Feedback
    back_propagated_field_prv = zeros(Nx, Ny, (Nz - 1));
    back_propagated_field_prv_prv = zeros(Nx, Ny, (Nz - 1));

    % Allocating space for Error rejection
    back_propagated_field_unaltered = zeros(Nx, Ny, Nz);
    RMS_error = zeros(Nz, Nz);

    % Starting Kalman weight
    kalman_weight = params.kalman_start_weight;

    % Pre-calculate propagation terms for each z position
    propagation_terms = zeros(Nx, Ny, Nz);
    for zi = 2:Nz
        kz = sqrt(k^2 - Kx.^2 - Ky.^2);
        kz = kz - 1i*imag(kz);
        propagation_terms(:,:,zi) = exp(1i * kz * params.z(zi));
    end
    
    %% Main iteration loop
    for iter = 1:params.max_iterations
        field = field_next;

        % Resetting field sum
        field_sum = zeros(Nx, Ny);

        for zi = 2:Nz
            %% Forward propagation to measurement plane zi
            % Decompose into plane waves
            field_f = fftshift(fft2(field));
            
            % Propagate each plane wave
            propagated_field_f = field_f .* propagation_terms(:,:,zi);
            
            % Transform back to spatial domain
            propagated_field = ifft2(ifftshift(propagated_field_f));

            %% Apply constraints at measurement plane
            % Extract amplitude of propegated field
            measured_amplitude = sqrt(measured_intensities(:,:,zi));
            propagated_field_amplitude = abs(propagated_field);
            new_field_amplitude = measured_amplitude;
            %% Apply mask function to reject errorous measurements
            % Check for overexposure
            mask = (measured_intensities(:,:,zi) < params.intensity_threshold);
            
            
            % Check if within noise level
            noise_mask = abs((measured_intensities(:, :, zi) - (abs(propagated_field)).^2)) > params.noise_level;
            mask = mask & noise_mask;

            % Apply mask
            new_field_amplitude(~mask) = propagated_field_amplitude(~mask);

            % Construct new field
            new_field = new_field_amplitude.* exp(1i * angle(propagated_field));
            
            %% Back propagation to sample plane
            % Decompose into plane waves
            new_field_f = fftshift(fft2(new_field));
            
            % Back propagate each plane wave
            back_propagated_field_f = new_field_f .*conj(propagation_terms(:,:,zi));
            
            % Transform back to spatial domain
            back_propagated_field = ifft2(ifftshift(back_propagated_field_f));
            
            % Saving fields for Error check
            back_propagated_field_unaltered(:, :, zi) = back_propagated_field;

            %% Accumulate fields for averaging with Double Feedback
            if iter < 2
                field_sum = field_sum + back_propagated_field;
                back_propagated_field_prv_prv(:, :, (zi - 1)) = back_propagated_field_prv(:, :, (zi - 1));
                back_propagated_field_prv(:, :, (zi - 1)) = back_propagated_field;
            else
                back_propagated_field_DF = (params.feedback_weight1 + params.feedback_weight2 + 1) * back_propagated_field - params.feedback_weight1 * back_propagated_field_prv(:, :, (zi - 1)) - params.feedback_weight2 * back_propagated_field_prv_prv(:, :, (zi - 1));
                field_sum = field_sum + back_propagated_field_DF;
                back_propagated_field_prv_prv(:, :, (zi - 1)) = back_propagated_field_prv(:, :, (zi - 1));
                back_propagated_field_prv(:, :, (zi - 1)) = back_propagated_field_DF;
            end
        end

        %% Update field estimate
        avg_field = field_sum / (Nz - 1); % Average the fields
        phase_avg = angle(avg_field); % Extract the phase
        amplitude_avg = abs(avg_field);
        amplitude_kalman = kalman_weight * sqrt(measured_intensities(:,:,1)) + (1-kalman_weight) * amplitude_avg; % Kalman filtering
        field_next = amplitude_kalman .* exp(1i * phase_avg);

        %% Check convergence and plot results
        if iter > 1
            amp_diff = amplitude_avg(xi, yi) - sqrt(measured_intensities(xi, yi, 1));
            Nmean = size(xi, 1) * size(yi, 1);
            Mamp2 = measured_intensities(xi, yi, 1);
            avgerror = sqrt((sum(amp_diff(:).^2)) / Nmean);
            converg = avgerror / sqrt((sum(Mamp2(:))) / Nmean);
            tiledlayout(1,2);

            % phase plot from algorithm
            nexttile;
            imagesc(params.x(xi) * 1e6, params.y(yi) * 1e6, phase_avg(xi, yi));
            axis image;
            title('Phasemap Algorithm');
            ylabel(colorbar, 'Phase', 'FontSize', 12)
            colormap('jet'); % You can change to 'hot', 'cool', 'parula', etc.
            xlabel('x (um)');
            ylabel('y (um)');
            set(gca, 'YDir', 'normal');
            
            % intensity squared of algorithm
            nexttile;
            imagesc(params.x(xi) * 1e6, params.y(yi) * 1e6, amplitude_avg(xi, yi));
            axis image;
            ylabel(colorbar, 'Intensity', 'FontSize', 12)
            title('Near field intensity');
            colormap('jet'); % You can change to 'hot', 'cool', 'parula', etc.
            xlabel('x (um)');
            ylabel('y (um)');
            set(gca, 'YDir', 'normal');
            drawnow;

            % Check if the convergence thershold is reached
            if converg < params.convergence_threshold
                fprintf('Converged after %d iterations\n', iter);
                break;
            end
        end

        % Printing convergence data
        if iter > 1
            fprintf('Iteration %d, Convergence metric: %f\n', iter, converg);
        else
            fprintf('First iteration, so no convergence\n');
        end

        %% Kalman weight update
        kalman_weight = params.kalman_start_weight + (params.kalman_end_weight - params.kalman_start_weight) * (iter / params.max_iterations); 
    end

    % Max itterations reached
    if converg > params.convergence_threshold
        fprintf("Max iterations reached\n")
    end

    %% Error check
    back_propagated_field_unaltered(:, :, 1) = field_next;
    for k = 1:Nz
        for l = 1:Nz
            temp_error = back_propagated_field_unaltered(:, :, k) - back_propagated_field_unaltered(:, :, l);
            RMS_error(k, l) = sqrt(mean(abs(temp_error(:)).^2));
        end
    end

    %% Final results
    final_phase = angle(field_next);
    reconstructed_field = amplitude_avg .* exp(1i * final_phase);
end