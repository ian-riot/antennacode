%optimise patch for vswr

clear all
close all


DO_5GHZ = 0
DO_GENETIC_ALGORITHM = 0
USE_GRADIENTS = 0


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

f = flow*1e9:fstep*1e6:fhigh*1e9;



%dielectric: FR4
epsr = 4.2;
dissipation_factor = 0.008;
h = 1.56; %mm
h_m = h/1000; %m

targetZ = 100;


h = 1.6; %mm
L2 = 26.13; %mm
W2 = 33.15; %mm
y02 = 2.74; %mm


% 5 GHz patch
L5 = 11.7; %mm
W5 = 13.2; %mm
y05 = 0.335; %mm


w2 = 2; % stripline multiplier

if DO_5GHZ
    f0_hz = 5.6e9; %GHz
    L_m = L5/1000;
    W_m = W5/1000;
    y0 = y05/1000;
else
    f0_hz = 2.44e9; %GHz
    L_m = L2/1000;
    W_m = W2/1000;
    y0 = y02/1000;
end

h_m = h/1000;

Z0 = 100;
dd = h;

d = dielectric('FR4');
d.EpsilonR = epsr;
d.LossTangent = 0.008;
d.Thickness = h_m;

 Zx=Z0;
 Er = epsr;

 Aa=Zx/60;
 Ab=sqrt((Er+1)/2);
 Ac=(Er-1)/(Er+1);
 Ad=(0.23+(0.11/Er));
 A=Aa*Ab+Ac*Ad;
 Wdr1=(8*exp(A))/(exp(2*A)-2);      % W/d ratio < 2


 B=(377*pi)/(2*Zx*sqrt(Er));
 Wdr2a=(2/pi);
 Wdr2b=(B-1-log(2*B-1));
 Wdr2c=(Er-1)/(2*Er);
 Wdr2d=(log(B-1)+0.39-(0.61/Er));
 Wdr2=Wdr2a*(Wdr2b+Wdr2c*Wdr2d);    % W/d ratio > 2
 
 if Wdr1<2  
   Wdr=Wdr1;
 else
   Wdr=Wdr2;
 end
 
 % Calculate line width (mm)
 w=Wdr*dd; 

 w_m = w/1000;
 
w2_m = w2*w_m;




d = dielectric('FR4');
d.EpsilonR = epsr;
d.LossTangent = dissipation_factor;
d.Thickness = h_m;
m = metal();


                                 
insetpatch = patchMicrostripInsetfed('Length',L_m,...
                                     'Width',W_m,...
                                     'Height',h_m,...
                                     'Substrate',d,...
                                     'PatchCenterOffset',[0 0],...
                                     'FeedOffset',[-L_m 0],...
                                     'StripLineWidth',w_m,...
                                     'NotchLength',y0,...
                                     'GroundPlaneLength',L_m*2,...
                                     'GroundPlaneWidth',W_m*2,...
                                     'NotchWidth',2*w_m);                         

                                 
show(insetpatch);
output = vswr(insetpatch,f0_hz,100);

vals_opt(1) = L_m;
vals_opt(2) = W_m;
vals_opt(3) = y0;
%vals_opt(4) = 2*w_m;
% vals_opt(5) = w_m;

vals_opt_min(1) = L_m-L_m/5;
vals_opt_min(2) = W_m-W_m/5;
vals_opt_min(3) = 1e-4;
% vals_opt_min(5) = w_m-1e-4;
%vals_opt_min(4) = (w_m-1e-4)*2;

vals_opt_max(1) = L_m+L_m/5;
vals_opt_max(2) = W_m+W_m/5;
vals_opt_max(3) = L_m-L_m/4;
% vals_opt_max(5) = w_m+1e-4;
%vals_opt_max(4) = (w_m+1e-4)*2;

% vals_opt = [9.9660   17.6272    5.6265    4.8679    2.5989]./1e3;
% vals_opt = [9.9660   17.6272    y0*1e3    4.8679    2.5989]./1e3;
% vals_opt = [9.9660   17.6272    y0*1e3     ]./1e3;
%vals_opt = [11.3083   18.5259    3.5170     ]./1e3;
% vals_opt = [28.2 36.6 8.8]*1e-3;
% take parameters out
FitnessFunction = @(vals_opt) optvswr(insetpatch,vals_opt,f0_hz,Z0);

options = optimset('Display','iter','TolFun',0.1);
options3 = optimset('Display','iter','TolFun',0.1,'UseParallel', true);
options2 = optimoptions(@ga,'Display','iter','UseParallel', true);
if DO_GENETIC_ALGORITHM
    vals_opt = ga(@(vals_opt) optvswr(insetpatch,vals_opt,f0_hz,Z0),length(vals_opt),[],[],[],[],vals_opt_min,vals_opt_max,[],options2);
    vals_opt = fminsearch(@(vals_opt) FitnessFunction(vals_opt),vals_opt,options);
else
    if USE_GRADIENTS
        vals_opt = fminunc(@(vals_opt) FitnessFunction(vals_opt),vals_opt,options3);
        vals_opt = fminsearch(@(vals_opt) FitnessFunction(vals_opt),vals_opt,options);
    else
        vals_opt = fminsearch(@(vals_opt) FitnessFunction(vals_opt),vals_opt,options);
    end
end


 insetpatch.Length = vals_opt(1);
    insetpatch.Width = vals_opt(2);
    insetpatch.NotchLength = vals_opt(3);
%     insetpatch.NotchWidth = vals_opt(4);

h = figure(1);
show(insetpatch);
txt1 = ['L: ' num2str(vals_opt(1)*1000) ' mm'];
txt2 = ['W: ' num2str(vals_opt(2)*1000) ' mm'];
txt3 = ['Notch Depth: ' num2str(vals_opt(3)*1000) ' mm'];
%txt4 = ['Notch Width: ' num2str(vals_opt(4)*1000) ' mm'];
% txt5 = ['50 ohm width: ' w_m*1000) ' mm'];
txt = {txt1,txt2,txt3};
annotation(h,'textbox', [0.05, 0.9, 0.1, 0.1], 'String',txt,'FitBoxToText','on');                                 
figure(2)
vswr(insetpatch,f,Z0)


function output = optvswr(patch,vals_opt,freq,Z0)
    patch.Length = max(1e-5,vals_opt(1));
    patch.Width = max(1e-5,vals_opt(2));
    patch.NotchLength = max(1e-5,vals_opt(3));
    patch.NotchLength = min(patch.Length,patch.NotchLength);
    
    patch.Length = vals_opt(1);
    patch.Width = vals_opt(2);
    patch.NotchLength = vals_opt(3);
   
    
    %patch.NotchWidth = vals_opt(4);
%     patch.StripLineWidth = vals_opt(5);
    h = figure(1);
    show(patch);
    txt1 = ['L: ' num2str(vals_opt(1)*1000) ' mm'];
    txt2 = ['W: ' num2str(vals_opt(2)*1000) ' mm'];
    txt3 = ['Notch Depth: ' num2str(vals_opt(3)*1000) ' mm'];
    %txt4 = ['Notch Width: ' num2str(vals_opt(4)*1000) ' mm'];
%     txt5 = ['50 ohm width: ' num2str(vals_opt(5)*1000) ' mm'];
    txt = {txt1,txt2,txt3};
    annotation(h,'textbox', [0.05, 0.9, 0.1, 0.1], 'String',txt,'FitBoxToText','on');
    pause(0.1);
    
    output = vswr(patch,freq,Z0);
    disp(['L:' num2str(vals_opt(1)*1000) ' W: ' num2str(vals_opt(2)*1000) ' y0: ' num2str(vals_opt(3)*1000) ' = ' num2str(output) ]);
end


