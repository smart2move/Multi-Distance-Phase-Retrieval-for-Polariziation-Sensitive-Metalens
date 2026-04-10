clear all;
% ! Dont forget to set the name of the file you want to save it in !

%% Load data and parameters
load("Example.mat")
Pixel_size = 0.48e-6; % Original pixel size
Pixel_size_new = 2.4e-6; % Wanted pixel size (must be an integer ratio of the original pixel size)
square = 0; % Do you want to square the data?
Square_offset = 0; % + or - offset for squaring

%% Setting size parameters
[Nx, Ny, Nz] = size(E2_esMDPR);

Pixel_increase = Pixel_size_new / Pixel_size; % ratio of pixel size increase

Nx_new = floor(Nx / Pixel_increase);  % New size of data
Ny_new = floor(Ny / Pixel_increase);

New_data = zeros(Nx_new, Ny_new, Nz); % Allocating space


for i = 1:Nz
    for j = 1:Ny_new
        for k = 1:Nx_new

            sum = 0;
            Py = (j - 1) * Pixel_increase  + 1; % Setting starting point for collecting data for averaging
            Px = (k - 1) * Pixel_increase  + 1;

            for l = Py:(Py - 1 + Pixel_increase)
                for m = Px:(Px - 1 + Pixel_increase)
                    sum = sum + E2_esMDPR(m, l, i); % summing data over new pixel size
                end
            end
            New_data(k, j, i) = sum / Pixel_increase^2; % averaging over new pixel size
        end
    end
end

xspace = Pixel_size_new * (Nx_new - 1);
yspace = Pixel_size_new * (Ny_new - 1);

xf = linspace(-xspace/2, xspace/2, Ny_new);
yf = linspace(-yspace/2, yspace/2, Nx_new);

E2_esMDPR = New_data;

%% squaring data
if square == 1
    if Nx_new < Ny_new
        New_data_squared = zeros(Nx_new, Nx_new, Nz);
        N_diff = Ny_new - Nx_new;

        if mod(N_diff, 2) == 0
            for i = 1:Nz
                New_data_squared(:, :, i) = New_data(:, (Square_offset + N_diff/2):(Square_offset + Ny_new - N_diff/2 - 1), i);
            end
        else
            for i = 1:Nz
                New_data_squared(:, :, i) = New_data(:, ((Square_offset + ceil(N_diff/2)):(Square_offset + Ny_new - ceil(N_diff/2))), i);
            end
        end

        Ny_new = Nx_new;
    end
    if Nx_new > Ny_new
        New_data_squared = zeros(Ny_new, Ny_new, Nz);
        N_diff = Nx_new - Ny_new;

        if mod(N_diff, 2) == 0
            for i = 1:Nz
                New_data_squared(:, :, i) = New_data((Square_offset + N_diff/2):(Square_offset + Nx_new - N_diff/2 - 1), :, i);
            end
        else
            for i = 1:Nz
                New_data_squared(:, :, i) = New_data(((Square_offset + ceil(N_diff/2)):(Square_offset + Nx_new - ceil(N_diff/2))), :, i);
            end
        end

        Nx_new = Ny_new;
    end

    xspace = Pixel_size_new * (Nx_new - 1);
    yspace = Pixel_size_new * (Ny_new - 1);
    
    xf = linspace(-xspace/2, xspace/2, Nx_new);
    yf = linspace(-yspace/2, yspace/2, Ny_new);


    E2_esMDPR = New_data_squared;
end

%% Saving data
save("Example.mat", "E2_esMDPR", 'xf', 'yf', 'zf', 'xspace', 'yspace');
