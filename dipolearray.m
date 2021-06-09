%optimise patch for vswr

clear all
close all


DO_5GHZ = 1
DO_VERTICAL = 1


c = 299792458;



if DO_5GHZ
    flow = 5.5; %GHz
    fhigh = 5.7; %GHz
    f0 = (fhigh-flow)/2 + flow; %GHz
    f0_hz = f0*1e9; %Hz
    fstep = 20; %MHz

else
    flow = 2.4; %GHz
    fhigh = 2.5; %GHz
    f0 = (fhigh-flow)/2 + flow; %GHz
    f0_hz = f0*1e9; %Hz    
    fstep = 10; %MHz    


end

nocols = 1;
norows = 4;

if (nocols==1)
    DOLINARRAY = 1;
else
    DOLINARRAY = 0;
end

L_des = 24.10;
W_des = 22.43;
Ld_des = 16.82;
w_des = 2.04;

% Transform to 0mm gap feed
L = (L_des + w_des/2)/1000;
W = W_des/1000;
Ld = (Ld_des + w_des/2)/1000;
g = 0.1e-3;
fw = 0.1e-3


% L_mm = 24.95; 
% L_mm = 19.16
% W_mm = 23.53; 
% W_mm = 18.87;
% Ld_mm = 23.07; 
% Ld_mm = 14.19;
% g_mm = 0.1; 
% % g_mm = 2.53;
% fw_mm = 0.1; %2.53;
s_mm = 15;

% L = L_mm/1000;
% W = W_mm/1000;
% Ld = Ld_mm/1000;
% g = g_mm/1000;
% fw = fw_mm/1000;
s = s_mm/1000;

%dielectric: FR4
epsr = 4.2;
dissipation_factor = 0.008;
h = 1.6; %mm
h_m = h/1000; %m

d = dielectric('FR4');
d.EpsilonR = epsr;
d.LossTangent = dissipation_factor;
d.Thickness = h_m;



lambda_0 = physconst('lightspeed')/1e9;
lambda_d = lambda_0/sqrt(epsr);

f = flow*1e9:fstep*1e6:fhigh*1e9;

if DO_VERTICAL
    bladeDipole = dipoleBlade('Length', L, 'width', W, 'TaperLength', Ld,...
                           'FeedWidth', fw,'FeedGap', g,'Tilt', 90,'TiltAxis',[0 1 0]);
else
    bladeDipole = dipoleBlade('Length', L, 'width', W, 'TaperLength', Ld,...
                           'FeedWidth', fw,'FeedGap', g,'Tilt',[90 90],'TiltAxis',[0 1 0; 0 0 1]);
end
                       
% b = bowtieTriangular('Length',0.05)
if DO_VERTICAL                       
    vertdist = L*2.1;
    hordist = W*1.3; 
else
    vertdist = W*1.3;
    hordist = L*2.1;    
end
                       

% ant = reflector('Exciter',rectArr)

% d = dipole('Length',0.15,'Width',0.015, 'Tilt',90,'TiltAxis',[0 1 0]);
ww = 0.05*nocols;
if DO_VERTICAL 
    rf = reflector('Exciter',bladeDipole,'GroundPlaneLength',vertdist, 'GroundPlaneWidth',.1,...
              'Spacing',s);
else
    rf = reflector('Exciter',bladeDipole,'GroundPlaneLength',vertdist, 'GroundPlaneWidth',.1,...
              'Spacing',s);
end
figure(1);
if DOLINARRAY
    linArr = linearArray('Element',rf,'elementSpacing',vertdist);
    linArr.NumElements = norows;
    show(linArr)
else
    rectArr = rectangularArray('Element',rf,'RowSpacing',hordist,'ColumnSpacing',vertdist);   
    rectArr.Size = [nocols  norows];
    show(rectArr);
end




% rectArr.GroundPlaneLength = lambda_d;
% rectArr.GroundPlaneWidth = lambda_d/4;


% 
 h = figure(2)
 if DOLINARRAY
     pattern(linArr,f0_hz)
 else
     pattern(rectArr,f0_hz)
 end
% 
% figure(3);
% dira = patternAzimuth(linArr,f0_hz,90);


% figure(3);
% dirb = patternElevation(linArr,f0_hz,90,'Elevation',angs);
if DOLINARRAY
    dirb = patternElevation(linArr,f0_hz,90);
else
    dirb = patternElevation(rectArr,f0_hz,90);
end
[maxer maxi] = max(dirb);
three = maxer-3;
new = (dirb>=three)
azbeam = sum(new)

% figure(4);
% patternElevation(linArr,f0_hz,[-180 -135 -90 -45 0 45 90 135 180]);
if DOLINARRAY
    dirb = patternElevation(linArr,f0_hz,0);
else
    dirb = patternElevation(rectArr,f0_hz,0);
end
[maxer maxi] = max(dirb);
three = maxer-3;
new = (dirb>=three)
elbeam = sum(new)

txt1 = ['Center Freq: ' num2str(f0) ' GHz'];
txt2 = ['azimuth beamwidth: ' num2str(azbeam) ' deg'];
txt3 = ['elevation beamwidth: ' num2str(elbeam) ' deg'];

txt = {txt1,txt2,txt3};
annotation(h,'textbox', [0.3, 0.1, 0.1, 0.1], 'String',txt,'FitBoxToText','on');  



return;



