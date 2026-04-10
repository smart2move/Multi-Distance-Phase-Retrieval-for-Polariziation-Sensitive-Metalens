clear all;
%% image convertion to matlab matrix
addpath 'Example_Location/Fiji/scripts' % Update for your ImageJ installation as appropriate
ImageJ;
%% parameters
Distance_z = 0.06e-3; % Distance between measurements
Pixel_size = 0.48e-6;
%% getting data from image
% Before this step, open the image in imagej
IJM.getDataset;

%% Declare matrix with z dimension
E2_epMDPR = zeros(1280, 1024, 20); 

%% Put all data into the matrix
E2_epMDPR(:, :, 20) = Es_650nm_d1_14_tif;

%% All other variables in the dataset
zf = zeros(20, 1);
for i = 2:20
    zf(i) = zf(i - 1) + Distance_z;
end

[Nx, Ny] = size(Es_650nm_d0_00_tif);

xspace = Pixel_size * (Nx - 1);
yspace = Pixel_size * (Ny - 1);

xf = linspace(-xspace/2, xspace/2, Nx);
yf = linspace(-yspace/2, yspace/2, Ny);

%% Save variables in new matlab file
save("Example.mat"); 
