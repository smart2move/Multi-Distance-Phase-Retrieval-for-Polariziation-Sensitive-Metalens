clear all;
clc;
%% Parameters
lamda = 587e-9; % Wavelenght
NA = 0; % NA (put 0 if unkown)
f = 10e-3; % Focal length (put 0 if unkown)
R = 100e-6; % Lens radius
a = 4.2e-6; % Pixel size
N = 308; % Amount of pixels in one dimention

%% Calculations
if(f == 0) % If focal lenght is unknown
    f = R / (tan(asin(NA))); % Focal lenght needed to achieve NA
end
if(NA == 0) % If NA is unkown
    NA = sin(atan(R/f));
end

amax = lamda/(2*NA); % Maximal pixelsize
Nmin = 4*R*NA/lamda; % Minimal amouth of pixels in 1 dimension
Zmax = N*a^2/lamda; % Maximal total scanning lenght


fprintf("Maximal pixel size: %.0f um\n", amax*1e6);
fprintf("Minimal amouth of pixels in 1 dimension: %.0f\n", ceil(Nmin));
fprintf("Maxial scanning distance: %.0f um\n", Zmax*1e6);
