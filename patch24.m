% single 2.45GHz horizontal patch

% from http://www.emtalk.com/mpacalc.php

clear all
close all

c = 299792458;
mu0 = 1.256e-6;
eps0 = 8.854e-12;

f0 = 2.45; %GHz
flow = 2.4; %GHz
fhigh = 2.5; %GHz

%dielectric: FR4
epsr = 4.4;
h = 1.56; %mm

targetZ = 50;



f = flow*1e9:10e6:fhigh*1e9;

% Do some stripline calcs
h_m = h/1000
f0_hz = f0*1e9;
W_m = (c/(2*f0_hz))*(sqrt(2/(epsr+1)))
epseff = (epsr+1)/2 + ((epsr-1)/2)*(1+12*h_m/W_m)^(-0.5);
deltaL = h_m*0.412*((epseff+0.3)*((W_m/h_m)+0.264))/((epseff-0.258)*((W_m/h_m)+0.8));
L = c/(2*f0_hz*sqrt(epseff));
L*1000/2
L_m = L - 2*deltaL
y0 = 1e-4*(0.001699*(epsr^7)+0.13761*(epsr^6)-6.1783*(epsr^5)+93.187*(epsr^4)-682.69*(epsr^3)+2561.9*(epsr^2)-4043*(epsr^1)+6697*(epsr^0))*L_m/2;
syms eff tz hh ww p
eff = epseff;
tz = targetZ;
hh = h_m;
p = pi;
eqn = 120*p/(sqrt(eff)*(1.393+(ww/hh)+(2/3)*log((ww/hh)+1.444))) == tz

w_m = eval(solve(eqn));
if (w_m<0)
    disp('desired impedance too high');
    return;
end

Z0 = 120*pi/(sqrt(epseff)*(1.393+(w_m/h_m)+(2/3)*log((w_m/h_m)+1.444)))

% do array spacing
vertdist = L*2;
hordist = L*1.5;


%stripline
% w = 1; %mm
% 
% L = 29; %mm
% W = 37.4; %mm
% l2 = 8; %Notch length in mm


w = w_m*1000;
L = L_m*1000;
W = W_m*1000;
l2 = y0*1000;
w2 = 2*w; %Notch width in mm




d = dielectric('FR4');
m = metal();
le = lumpedElement('Impedance',50,'Frequency',f0*1e9,'Location',[-W/2000 0 0]);

% insetpatch = patchMicrostripInsetfed('Length',L/1000,...
%                                      'Width',W/1000,...
%                                      'Height',h/1000,...
%                                      'GroundPlaneLength',L*2/1000,...
%                                      'GroundPlaneWidth',W*2/1000,...
%                                      'Substrate',d,...
%                                      'PatchCenterOffset',[0 0],...
%                                      'FeedOffset',[-L/2000 0],...
%                                      'StripLineWidth',w/1000,...
%                                      'NotchLength',l2/1000,...
%                                      'NotchWidth',w2/1000,...
%                                      'Load',le);

insetpatch = patchMicrostripInsetfed('Length',L/1000,...
                                     'Width',W/1000,...
                                     'Height',h/1000,...
                                     'GroundPlaneLength',L*2/1000,...
                                     'GroundPlaneWidth',W*2/1000,...
                                     'Substrate',d,...
                                     'PatchCenterOffset',[0 0],...
                                     'FeedOffset',[-L/1500 0],...
                                     'StripLineWidth',w/1000,...
                                     'NotchLength',l2/1000,...
                                     'NotchWidth',w2/1000);
                                     
                                 
insetpatch = patchMicrostripInsetfed('Length',L/1000,...
                                     'Width',W/1000,...
                                     'Height',h/1000,...
                                     'Substrate',d,...
                                     'PatchCenterOffset',[0 0],...
                                     'FeedOffset',[-L/1000 0],...
                                     'StripLineWidth',w/1000,...
                                     'NotchLength',l2/1000,...
                                     'NotchWidth',w2/1000);                         
                                 
show(insetpatch)


% array stuff
arrayObject = rectangularArray;

arrayObject.Size = [2  4];
arrayObject.RowSpacing = vertdist;
arrayObject.ColumnSpacing = hordist;

%Define Array Elements
Element1 = patchMicrostrip;

% approximate the patch info to the inset fed info
Element1.Length = insetpatch.Width;
Element1.Width = insetpatch.Length;
%Element1.FeedLocation = [y0 - insetpatch.Length/2 0 0 ];
Element1.Height = insetpatch.Height;
Element1.GroundPlaneLength = .06;
Element1.GroundPlaneWidth = .07;
Element1.FeedOffset = [0 y0 - insetpatch.Length/2];
Element1.Substrate = d;
Element1.Substrate.EpsilonR = epsr;

arrayObject.Element = [Element1 ];





figure;
show(arrayObject) 
axis ( [-200 200 -200 200 ]);

figure
patternAzimuth(arrayObject,f0_hz,89,'Azimuth',[-180:30:180]);
return

figure
patternElevation(arrayObject,f0_hz,0);

return;

figure;
pattern(arrayObject,f0_hz);

%optAnt = optimize(insetpatch, f0*1e9, 'maximizeGain',{'Length', 'Width'}, {35e-3 20e-3; 40e-3 55e-3})


%vswr(insetpatch,f,50);