% blade dipole

close all
clear all

L = 16.34e-3;         % Total half-length of blade dipole
W = 25.63e-3;         % Width of blade dipole
Ld = 12.88e-3;        % Half-length excluding impedance match taper
fw = 0.1e-3;          % Feed width
g = 0.1e-3;           % Feed gap



bladeDipole = dipoleBlade('Length', L, 'width', W, 'TaperLength', Ld,...
                           'FeedWidth', fw,'FeedGap', g);

figure(1);
show(bladeDipole);

fmin = 2e9;
fmax = 6e9;
Nfreq = 21;
freq = linspace(fmin,fmax,Nfreq);
s_blade = sparameters(bladeDipole,freq);

s11 = s_blade.Parameters

figure(2)

rfplot(s_blade,1,1);

