clear all;
%% Unit cell max size 
f = 1200e-6;            % focal lenght
r = 110e-6;            % lens radius
lamda = 800e-9;        % wavelength
n = 1.4533;              % refractive index

umax(1) = lamda/n;       % Maximal unit cell size 0th order diffraction
NA = sin(atan(r/f));   % Calculate NA
umax(2) = lamda/(2*NA);  % Maximal unit cell size Nyquist

Umax = min(umax);
