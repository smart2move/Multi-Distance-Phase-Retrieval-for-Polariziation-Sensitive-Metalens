clear all
%% Load the raw phase data from Lumerical
load('Example.mat');  % Your Lumerical export file

%% Do you want plots?
Imageon = 1;

%% Display raw data information
fprintf('Data dimensions: %dx%d\n', size(phase_ss));
fprintf('X-span range: %.3f to %.3f um\n', min(xspan)*1e6, max(xspan)*1e6);
fprintf('Y-span range: %.3f to %.3f um\n', min(yspan)*1e6, max(yspan)*1e6);

%% 2D Phase Unwrapping
fprintf('Performing 2D phase unwrapping...\n');

phase_ss_unwrapped = unwrap(phase_ss, [], 2);  % Unwrap along rows (Y-span)
phase_ss_unwrapped = unwrap(phase_ss_unwrapped, [], 1);  % Then along columns (X-span)

phase_pp_unwrapped = unwrap(phase_pp, [], 1);  % Unwrap along columns (X-span)
phase_pp_unwrapped = unwrap(phase_pp_unwrapped, [], 2);  % Then along rows (Y-span)

% Normalize to start at 0
phase_ss_unwrapped = phase_ss_unwrapped - phase_ss_unwrapped(1,1);
phase_pp_unwrapped = phase_pp_unwrapped - phase_pp_unwrapped(1,1);

% Transpose
phase_ss_unwrapped = transpose(phase_ss_unwrapped);
phase_pp_unwrapped = transpose(phase_pp_unwrapped);

% Check if it worked
if(Imageon == 1)
    figure('Position', [100, 100, 1000, 400]);
    
    subplot(1,2,1)
    imagesc(yspan * 1e6, xspan * 1e6, phase_ss_unwrapped);
    axis xy;
    title('Unwrapped phase Es');
    ylabel(colorbar, 'Phase', 'FontSize', 12)
    colormap('jet'); % You can change to 'hot', 'cool', 'parula', etc.
    xlabel('y (um)');
    ylabel('x (um)');
    
    
    subplot(1,2,2)
    imagesc(yspan * 1e6, xspan * 1e6, phase_pp_unwrapped);
    axis xy;
    title('Unwrapped phase Ep');
    ylabel(colorbar, 'Phase', 'FontSize', 12)
    colormap('jet'); % You can change to 'hot', 'cool', 'parula', etc.
    xlabel('y (um)');
    ylabel('x (um)');
end


%% Create denser dimension grids
xspan_dense = linspace(min(xspan), max(xspan), 51);  % Denser xspan grid
yspan_dense = linspace(min(yspan), max(yspan), 51);  % Denser yspan grid

[XSPAN_DENSE, YSPAN_DENSE] = meshgrid(xspan_dense, yspan_dense);

%% Interpolate phase maps to denser grid
fprintf('Interpolating 2D phase maps...\n');

phase_ss_dense = interp2(xspan, yspan, phase_ss_unwrapped, XSPAN_DENSE, YSPAN_DENSE, 'spline');
phase_pp_dense = interp2(xspan, yspan, phase_pp_unwrapped, XSPAN_DENSE, YSPAN_DENSE, 'spline');

if(Imageon == 1)
    figure('Position', [100, 100, 1000, 400]);
    
    subplot(1,2,1)
    imagesc(xspan * 1e6, yspan * 1e6, phase_ss_dense);
    axis xy;
    title('Interpolated phase Es');
    ylabel(colorbar, 'Phase', 'FontSize', 12)
    colormap('jet'); % You can change to 'hot', 'cool', 'parula', etc.
    xlabel('Length (um)');
    ylabel('Width (um)');
    
    subplot(1,2,2)
    imagesc(xspan * 1e6, yspan * 1e6, phase_pp_dense);
    axis xy;
    title('Interpolated phase Ep');
    ylabel(colorbar, 'Phase', 'FontSize', 12)
    colormap('jet'); % You can change to 'hot', 'cool', 'parula', etc.
    xlabel('Length (um)');
    ylabel('Width (um)');
end
%% Interpolation of E and H fields
fprintf('Interpolating E and H fields...\n');
x = E_ss.x;
y = E_ss.y;

Esx_dense = zeros(50,50,51,51);   % Allocation for field components
Esy_dense = zeros(50,50,51,51);
Esz_dense = zeros(50,50,51,51);

Epx_dense = zeros(50,50,51,51);
Epy_dense = zeros(50,50,51,51);
Epz_dense = zeros(50,50,51,51);

Hsx_dense = zeros(50,50,51,51);
Hsy_dense = zeros(50,50,51,51);
Hsz_dense = zeros(50,50,51,51);

Hpx_dense = zeros(50,50,51,51);
Hpy_dense = zeros(50,50,51,51);
Hpz_dense = zeros(50,50,51,51);


% interpolation over X span and Y span
for i = 1:size(x,2)
    for j = 1:size(y,2)
        Esx_dense(j,i,:,:) = interp2(xspan, yspan, squeeze(Esx(j,i,:,:)), XSPAN_DENSE, YSPAN_DENSE, 'spline');
        Esy_dense(j,i,:,:) = interp2(xspan, yspan, squeeze(Esy(j,i,:,:)), XSPAN_DENSE, YSPAN_DENSE, 'spline');
        Esz_dense(j,i,:,:) = interp2(xspan, yspan, squeeze(Esz(j,i,:,:)), XSPAN_DENSE, YSPAN_DENSE, 'spline');
        
        Epx_dense(j,i,:,:) = interp2(xspan, yspan, squeeze(Epx(j,i,:,:)), XSPAN_DENSE, YSPAN_DENSE, 'spline');
        Epy_dense(j,i,:,:) = interp2(xspan, yspan, squeeze(Epy(j,i,:,:)), XSPAN_DENSE, YSPAN_DENSE, 'spline');
        Epz_dense(j,i,:,:) = interp2(xspan, yspan, squeeze(Epz(j,i,:,:)), XSPAN_DENSE, YSPAN_DENSE, 'spline');

        Hsx_dense(j,i,:,:) = interp2(xspan, yspan, squeeze(Hsx(j,i,:,:)), XSPAN_DENSE, YSPAN_DENSE, 'spline');
        Hsy_dense(j,i,:,:) = interp2(xspan, yspan, squeeze(Hsy(j,i,:,:)), XSPAN_DENSE, YSPAN_DENSE, 'spline');
        Hsz_dense(j,i,:,:) = interp2(xspan, yspan, squeeze(Hsz(j,i,:,:)), XSPAN_DENSE, YSPAN_DENSE, 'spline');

        Hpx_dense(j,i,:,:) = interp2(xspan, yspan, squeeze(Hpx(j,i,:,:)), XSPAN_DENSE, YSPAN_DENSE, 'spline');
        Hpy_dense(j,i,:,:) = interp2(xspan, yspan, squeeze(Hpy(j,i,:,:)), XSPAN_DENSE, YSPAN_DENSE, 'spline');
        Hpz_dense(j,i,:,:) = interp2(xspan, yspan, squeeze(Hpz(j,i,:,:)), XSPAN_DENSE, YSPAN_DENSE, 'spline');
    end
end

%% Downsampling of fields
fprintf('Sampling E and H fields...\n');

ns = 3; % sampling per period
x_less = linspace(min(x), max(x), ns); % Sampled x coordinates
y_less = linspace(min(y), max(y), ns); % Sampled Y coordinates

[X_LESS, Y_LESS] = meshgrid(x_less, y_less);

Esx_sample = zeros(3,3,51,51);   % Allocation for field components
Esy_sample = zeros(3,3,51,51);
Esz_sample = zeros(3,3,51,51);

Epx_sample = zeros(3,3,51,51);
Epy_sample = zeros(3,3,51,51);
Epz_sample = zeros(3,3,51,51);

Hsx_sample = zeros(3,3,51,51);
Hsy_sample = zeros(3,3,51,51);
Hsz_sample = zeros(3,3,51,51);

Hpx_sample = zeros(3,3,51,51);
Hpy_sample = zeros(3,3,51,51);
Hpz_sample = zeros(3,3,51,51);

for i = 1:size(XSPAN_DENSE,1)
    for j = 1:size(YSPAN_DENSE,2)
        Esx_sample(:,:,j,i) = interp2(x, y, squeeze(Esx_dense(:,:,j,i)), X_LESS, Y_LESS, 'spline');
        Esy_sample(:,:,j,i) = interp2(x, y, squeeze(Esy_dense(:,:,j,i)), X_LESS, Y_LESS, 'spline');
        Esz_sample(:,:,j,i) = interp2(x, y, squeeze(Esz_dense(:,:,j,i)), X_LESS, Y_LESS, 'spline');
        
        Epx_sample(:,:,j,i) = interp2(x, y, squeeze(Epx_dense(:,:,j,i)), X_LESS, Y_LESS, 'spline');
        Epy_sample(:,:,j,i) = interp2(x, y, squeeze(Epy_dense(:,:,j,i)), X_LESS, Y_LESS, 'spline');
        Epz_sample(:,:,j,i) = interp2(x, y, squeeze(Epz_dense(:,:,j,i)), X_LESS, Y_LESS, 'spline');

        Hsx_sample(:,:,j,i) = interp2(x, y, squeeze(Hsx_dense(:,:,j,i)), X_LESS, Y_LESS, 'spline');
        Hsy_sample(:,:,j,i) = interp2(x, y, squeeze(Hsy_dense(:,:,j,i)), X_LESS, Y_LESS, 'spline');
        Hsz_sample(:,:,j,i) = interp2(x, y, squeeze(Hsz_dense(:,:,j,i)), X_LESS, Y_LESS, 'spline');

        Hpx_sample(:,:,j,i) = interp2(x, y, squeeze(Hpx_dense(:,:,j,i)), X_LESS, Y_LESS, 'spline');
        Hpy_sample(:,:,j,i) = interp2(x, y, squeeze(Hpy_dense(:,:,j,i)), X_LESS, Y_LESS, 'spline');
        Hpz_sample(:,:,j,i) = interp2(x, y, squeeze(Hpz_dense(:,:,j,i)), X_LESS, Y_LESS, 'spline');
    end
end

%% Create combined lookup table for both polarizations
fprintf('Creating combined polarization lookup table...\n');

% Create target phase grids for both polarizations
num_phase_points = 361;  % 1-degree resolution
target_phases_ep = linspace(0, 2*pi, num_phase_points);  % Ep target phases
target_phases_es = linspace(0, 2*pi, num_phase_points);  % Es target phases

[TARGET_PHASES_EP, TARGET_PHASES_ES] = meshgrid(target_phases_ep, target_phases_es);

% Initialize lookup tables
xspan_lut = zeros(size(TARGET_PHASES_EP));
yspan_lut = zeros(size(TARGET_PHASES_EP));
phase_error_lut = zeros(size(TARGET_PHASES_EP));

% For each target phase combination, find optimal pillar
for i = 1:numel(TARGET_PHASES_EP)
    target_phase_ep = TARGET_PHASES_EP(i);
    target_phase_es = TARGET_PHASES_ES(i);
    
    % Calculate phase errors for both polarizations
    phase_error_ep = abs(phase_pp_dense - target_phase_ep);
    phase_error_es = abs(phase_ss_dense - target_phase_es);
    
    % Combined error metric - finds pillars that work for BOTH polarizations
    total_error = phase_error_ep + phase_error_es;
    
    % Find the pillar with minimum combined error
    [min_error, min_idx] = min(total_error(:));
    [row, col] = ind2sub(size(total_error), min_idx);
    
    % Store results
    xspan_lut(i) = xspan_dense(col);
    yspan_lut(i) = yspan_dense(row);
    phase_error_lut(i) = min_error;
end

%% Create efficient lookup function
fprintf('Creating efficient lookup function...\n');

% Flatten the dense phase data for fast searching
all_xspan = XSPAN_DENSE(:);
all_yspan = YSPAN_DENSE(:);
all_phase_ep = phase_pp_dense(:);
all_phase_es = phase_ss_dense(:);

all_Esx = reshape(Esx_sample, size(Esx_sample, 1), size(Esx_sample, 2), []);
all_Esy = reshape(Esy_sample, size(Esy_sample, 1), size(Esy_sample, 2), []);
all_Esz = reshape(Esz_sample, size(Esz_sample, 1), size(Esz_sample, 2), []);

all_Epx = reshape(Epx_sample, size(Epx_sample, 1), size(Epx_sample, 2), []);
all_Epy = reshape(Epy_sample, size(Epy_sample, 1), size(Epy_sample, 2), []);
all_Epz = reshape(Epz_sample, size(Epz_sample, 1), size(Epz_sample, 2), []);

all_Hsx = reshape(Hsx_sample, size(Hsx_sample, 1), size(Hsx_sample, 2), []);
all_Hsy = reshape(Hsy_sample, size(Hsy_sample, 1), size(Hsy_sample, 2), []);
all_Hsz = reshape(Hsz_sample, size(Hsz_sample, 1), size(Hsz_sample, 2), []);

all_Hpx = reshape(Hpx_sample, size(Hpx_sample, 1), size(Hpx_sample, 2), []);
all_Hpy = reshape(Hpy_sample, size(Hpy_sample, 1), size(Hpy_sample, 2), []);
all_Hpz = reshape(Hpz_sample, size(Hpz_sample, 1), size(Hpz_sample, 2), []);


% Create lookup function
lookup_fun = @(target_phase_ep, target_phase_es) find_best_pillar_simple(target_phase_ep, target_phase_es, all_phase_ep, all_phase_es, all_xspan, all_yspan, all_Epx);

%% Test the lookup function
fprintf('Testing lookup function...\n');

% Example: find pillar for specific phase combination
% test_phase_ep_degree = 360;
% test_phase_es_degree = 360;


test_phase_es = 1*pi; 
test_phase_ep = 1.5*pi;  

[test_xspan, test_yspan, test_error, test_Epx, best_idx] = lookup_fun(test_phase_ep, test_phase_es);

fprintf('Test lookup results:\n');
fprintf('  Target: Ep phase=%.1f°, Es phase=%.1f°\n', test_phase_ep*180/pi, test_phase_es*180/pi);
fprintf('  Found:  xspan=%.1f nm, yspan=%.1f nm\n', test_xspan*1e9, test_yspan*1e9);
fprintf('  Phase error: %.1f°\n ', test_error*180/pi);
fprintf('  index: %.1f\n ', best_idx);

%% Visualization
if(Imageon == 1)
    figure('Position', [100, 100, 1000, 400]);
    
    % Required X-span for phase combinations
    subplot(1,3,1);
    imagesc(target_phases_ep*180/pi, target_phases_es*180/pi, xspan_lut*1e9);
    xlabel('Target Phase Ep (deg)'); ylabel('Target Phase Es (deg)');
    title('Required X-Span (nm)'); colorbar; axis equal;
    axis xy;
    xlim([0 360]);
    ylim([0 360]);
    
    % Required Y-span for phase combinations
    subplot(1,3,2);
    imagesc(target_phases_ep*180/pi, target_phases_es*180/pi, yspan_lut*1e9);
    xlabel('Target Phase Ep (deg)'); ylabel('Target Phase Es (deg)');
    title('Required Y-Span (nm)'); colorbar; axis equal;
    axis xy;
    xlim([0 360]);
    ylim([0 360]);
    
    % Phase error
    subplot(1,3,3);
    imagesc(target_phases_ep*180/pi, target_phases_es*180/pi, phase_error_lut);
    xlabel('Target Phase Ep (deg)'); ylabel('Target Phase Es (deg)');
    title('Phase Error (rad)'); colorbar; axis equal;
    axis xy;
    xlim([0 360]);
    ylim([0 360]);
end
%% Save the lookup table
fprintf('Saving combined polarization lookup table...\n');

z = E_ss.z; % Extracting z value for ansys

save('metalens_phase_and_fields_lut_v7.mat', ...
    'target_phases_ep', 'target_phases_es', ...
    'all_Esx','all_Esy','all_Esz', ...
    'all_Epx','all_Epy','all_Epz', ...
    'all_Hsx','all_Hsy','all_Hsz', ...
    'all_Hpx','all_Hpy','all_Hpz','X_LESS','Y_LESS', ...
    'xspan_lut', 'yspan_lut', 'phase_error_lut', ...
    'all_phase_ep', 'all_phase_es', 'all_xspan', 'all_yspan', ...
    'wavelength', 'period', 'height', 'mat_pillar', 'index_pillar', ...
    'mat_sub', 'index_sub','z', '-v7.3');

fprintf('Combined phase lookup table complete!\n');
fprintf('Use lookup_fun(target_phase_ep, target_phase_es) to find optimal pillars\n');

%% Helper function for pillar lookup
function [xspan, yspan, phase_error, Epx, best_idx ] = find_best_pillar_simple(...
    target_phase_ep, target_phase_es, all_phase_ep, all_phase_es, all_xspan, all_yspan, all_Epx)

    % Calculate phase errors for all pillars
    phase_error_ep = abs(all_phase_ep - target_phase_ep);
    phase_error_es = abs(all_phase_es - target_phase_es);
    
    % Combined error
    total_error = phase_error_ep + phase_error_es;
    
    % Find best pillar
    [phase_error, best_idx] = min(total_error);
    
    % Return dimensions
    xspan = all_xspan(best_idx);
    yspan = all_yspan(best_idx);

    % Return fields
    Epx = all_Epx(:,:,best_idx);


end
