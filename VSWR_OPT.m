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
epsr = 4.6;
dissipation_factor = 0.018;
h = 1.56; %mm
h_m = h/1000; %m

targetZ = 100;


% temp = f0_hz;
% f0_hz = fhigh*1e9;
% INITIAL CONDITIONS
% Do some stripline calcs
%Width of patch
% f_offset = 0;
W_m = (c/(2*(f0_hz)))*(sqrt(2/(epsr+1)))
% Length of patch
epseff = (epsr+1)/2 + ((epsr-1)/2)*(1+12*h_m/W_m)^(-0.5);
deltaL = h_m*0.412*((epseff+0.3)*((W_m/h_m)+0.264))/((epseff-0.258)*((W_m/h_m)+0.8));
L = c/(2*f0_hz*sqrt(epseff));
% L*1000/2
L_m = L - 2*deltaL

% depth of notch
y0 = 1e-4*(0.001699*(epsr^7)+0.13761*(epsr^6)-6.1783*(epsr^5)+93.187*(epsr^4)-682.69*(epsr^3)+...
      2561.9*(epsr^2)-4043*(epsr^1)+6697)*L_m/2;

% width of desired connect stripline
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

% lets solve for Z0= 70.7 ohms for quarte wave transformer
syms eff tz hh ww p
eff = epseff;
tz = 70.7;
hh = h_m;
p = pi;
eqn = 120*p/(sqrt(eff)*(1.393+(ww/hh)+(2/3)*log((ww/hh)+1.444))) == tz

w70_m = eval(solve(eqn));
if (w70_m<0)
    disp('desired impedance too high');
    return;
else
    disp(['70.7 ohm track width is ' num2str(w70_m*1000) ' mm']);
    disp(['quarter wavelength is ' num2str(L*1000/4) ' mm']);
end

% lets solve for Z0= 100 ohms for quarte wave transformer
syms eff tz hh ww p
eff = epseff;
tz = 100;
hh = h_m;
p = pi;
eqn = 120*p/(sqrt(eff)*(1.393+(ww/hh)+(2/3)*log((ww/hh)+1.444))) == tz

w100_m = eval(solve(eqn));
if (w100_m<0)
    disp('desired impedance too high');
    return;
else
    disp(['100 ohm track width is ' num2str(w100_m*1000) ' mm']);
end

% lets solve for Z0= 50 ohms for quarte wave transformer
syms eff tz hh ww p
eff = epseff;
tz = 50;
hh = h_m;
p = pi;
eqn = 120*p/(sqrt(eff)*(1.393+(ww/hh)+(2/3)*log((ww/hh)+1.444))) == tz

w50_m = eval(solve(eqn));
if (w50_m<0)
    disp('desired impedance too high');
    return;
else
    disp(['50 ohm track width is ' num2str(w50_m*1000) ' mm']);
end



% f0_hz = temp;


%stripline
% w = 1; %mm
% 
% L = 29; %mm
% W = 37.4; %mm
% l2 = 8; %Notch length in mm



w2_m = 2*w_m;




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
                                     'GroundPlaneLength',L_m*2 +1e-3,...
                                     'GroundPlaneWidth',L_m*2 + 1e-3,...
                                     'NotchWidth',2*w_m);                         

                                 
show(insetpatch);

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
vals_opt = [28.2 36.6 8.8]*1e-3;
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


