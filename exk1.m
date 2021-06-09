% Klopfenstein taper example1
%
% Using very short (lambda/100) transformers,
% a tapered match can be approximated.
% 
% The total length is given by : (Number of sections)*(Transformer Length)
% In this example N*Tlen=0.50 wavelength at Fo.
%
% Load=100 Ohm
% Line=50 Ohm
% Number of sections N=50
% Operating band ripple -35dB
%
% Transformer lengths Tlen=lambda/100
% Fo=1000 Mhz
% Plot 1 to 10000 MHz
%


clc;
close all;
help exk1

Zload=100;  % Load impedance, to matched (Ohms)
Zo=50;      % Characteristic impedance to match to (Ohms)
Fo=1000;    % Lower cut-off frequency (MHz)
F1=1;       % Start frequency for response plot (MHz)
F2=10000;   % Stop frequency for response plot (MHz)
Tlen=0.01;  % Transformer length as a fraction of wavelength
N=50;       % Number of transformer sections
RdB=-35;    % Operating band ripple (dB)
Er=3.48;    % Dielectric constant for microstrip
d=1.52;     % Thickness of substrate (mm)

Zlist=bklop(Zo,Zload,N,RdB);
bplot(Zlist,Tlen,Fo,F1,F2);
axis([0 F2 -50 0]);
filename=setfname;
bphysical(Zlist,Tlen,Fo,Er,d,filename);