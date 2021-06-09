%optimise patch for vswr

clear all
close all


DO_5GHZ = 1
USE_CALC = 1


c = 299792458;



if DO_5GHZ
    flow = 5.5; %GHz
    fhigh = 5.7; %GHz
    f0 = (fhigh-flow)/2 + flow; %GHz
    f0_hz = f0*1e9; %Hz
    fstep = 20; %MHz
    nocols = 2;
    norows = 4;
    L_mm = 12.11; %mm
    W_mm = 14.38; %mm
    y0_mm = 0.323; %mm
else
    flow = 2.4; %GHz
    fhigh = 2.5; %GHz
    f0 = (fhigh-flow)/2 + flow; %GHz
    f0_hz = f0*1e9; %Hz    
    fstep = 10; %MHz    
    nocols = 2;
    norows = 4;
    L_mm = 27.67; %mm
    W_mm = 35.37; %mm
    y0_mm = 2.58; %mm
end

L_m = L_mm/1000;
W_m = W_mm/1000;
y0_m = y0_mm/1000;

f = flow*1e9:fstep*1e6:fhigh*1e9;



%dielectric: FR4
epsr = 4.2;
dissipation_factor = 0.008;
h = 1.6; %mm
h_m = h/1000; %m

d = dielectric('FR4');
d.EpsilonR = epsr;
d.LossTangent = dissipation_factor;
d.Thickness = h_m;

targetZ = 100;



% do array spacing
vertdist = 0.0512; %L*2;
hordist = 0.0512; %L*1.5;
vertdist = L_m*2;
hordist = L_m*1.5;

% optimse patch overide
if (USE_CALC == 0)
    
end


% array stuff
arrayObject = rectangularArray;

arrayObject.Size = [nocols  norows];
arrayObject.RowSpacing = vertdist;
arrayObject.ColumnSpacing = hordist;

%Define Array Elements
Element1 = patchMicrostrip;

% approximate the patch info to the inset fed info
Element1.Length = W_m;
Element1.Width = L_m;
%Element1.FeedLocation = [y0 - insetpatch.Length/2 0 0 ];
Element1.Height = h_m;
Element1.GroundPlaneLength =  hordist;
Element1.GroundPlaneWidth = vertdist;
Element1.FeedOffset = [0 L_m/2 - y0_m ];
Element1.Substrate = d;


arrayObject.Element = [Element1 ];





figure;
show(arrayObject) 
%axis ( [-200 200 -200 200 ]);

% figur
% patternAzimuth(arrayObject,f0_hz,90);
%  return
% 
% figure
% patternElevation(arrayObject,f0_hz,0);
% 

% 
figure;
pattern(arrayObject,f0_hz,'Termination',targetZ);

