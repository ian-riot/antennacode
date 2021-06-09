% looking at microstrip widths and lengths

clear all

targetZ = 50

DO_5GHZ = 0



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
Er = epsr;
dissipation_factor = 0.018;
h = 1.6*2; %mm
d = h;
h_m = h/1000; %m
cu_h = 35e-6;


% basic start points:
H = h_m-2*cu_h;
W = 1e-3;  %2 mm guess

% work out both W/H<1 and W/H>=1

% W/H<1
% e_eff = ((epsr+1)/2) + ((epsr-2)/2)*((1/sqrt(1+12*H/W))+0.04*(1-(W/H)*(W/H)));
% Z0 = (60/sqrt(e_eff))*log((8*H/W)+(W/(4*H)));

% W/H>=1
% e_eff = ((epsr+1)/2) + (epsr-1)/(2*sqrt(1+12*H/W));
% Z0 = (120*pi)/(sqrt(e_eff)*((W/H) + 1.393+(2/3)*log((W/H)+1.444)));

% solve for both cases
syms eff tz hh ww p
oldW = 0;

while (abs(oldW - W) < 0.00001)
    oldW = W;
    if (W/H < 1)
        e_eff = ((epsr+1)/2) + ((epsr-1)/2)*((1/sqrt(1+12*(H/W)))+0.04*(1-((W/H)*(W/H))));
        e_eff = ((epsr+1)/2) + (epsr-1)/(2*sqrt(1+12*(H/W)));
%         Z0 = (60/sqrt(e_eff))*log((8*H/W)+(W/(4*H)));
    else
        e_eff = ((epsr+1)/2) + (epsr-1)/(2*sqrt(1+12*(H/W)));
%         Z0 = (120*pi)/(sqrt(e_eff)*((W/H) + 1.393+(2/3)*log((W/H)+1.444)));
    end
    eff = e_eff;
    tz = targetZ;
    hh = H;
    p = pi;
    if (W/H < 1)
        eqn = (60/sqrt(eff))*log((8*hh/ww)+(ww/(4*hh))) == tz;
    else
        eqn = 120*p/(sqrt(eff)*(1.393+(ww/hh)+(2/3)*log((ww/hh)+1.444))) == tz;
    end    
    
%     eqn = 120*p/(sqrt(eff)*(1.393+(ww/hh)+(2/3)*log((ww/hh)+1.444))) == tz;

    w_m = eval(solve(eqn));

    W = w_m(1);
    
    if (W/H < 1)
        e_eff = ((epsr+1)/2) + ((epsr-2)/2)*((1/sqrt(1+12*H/W))+0.04*(1-((W/H)*(W/H))));
        e_eff = ((epsr+1)/2) + (epsr-1)/(2*sqrt(1+12*(H/W)));
    else
        e_eff = ((epsr+1)/2) + (epsr-1)/(2*sqrt(1+12*(H/W)));
    end

end


% 
% % guess a start w_m
% w_m = 2e-3;  %2 mm
% 
% % Length of patch
% if ((w_m/h_m) < 1)
%     epseff = ((epsr+1)/2) + ((epsr-1)/2)*((1/(sqrt(1+12*h_m/w_m))+0.04*(1-(w_m/h_m)^2)));
% else
%     epseff = (epsr+1)/2 + ((epsr-1)/2)*(1+12*h_m/w_m)^(-0.5);
% end






% % width of desired connect stripline
syms eff tz hh ww p
oldw_m = 0;
% 
% while (oldw_m ~= w_m)
%     oldw_m = w_m;
%     eff = e_eff;
%     tz = targetZ;
%     hh = h_m;
%     p = pi;
%     eqn = 120*p/(sqrt(eff)*(1.393+(ww/hh)+(2/3)*log((ww/hh)+1.444))) == tz;
% 
%     w_m = eval(solve(eqn));
% 
% 
% 
% 
%     if ((w_m/h_m) < 1)
%         e_eff = ((epsr+1)/2) + ((epsr-1)/2)*((1/(sqrt(1+12*h_m/w_m))+0.04*(1-(w_m/h_m)^2)));
%     else
%         e_eff = (epsr+1)/2 + ((epsr-1)/2)*(1+12*h_m/w_m)^(-0.5);
%     end
% end
% 
% 
% 
% if (W<0)
%     disp('desired impedance too high');
%     return;
% else
%     L = c/(2*f0_hz*sqrt(e_eff));
%     disp([num2str(targetZ) ' ohm track width is ' num2str(W*1000) ' mm (WRONG)']);
%     disp(['quarter wavelength is ' num2str(L*1000/4) ' mm']);
% end

 Zx=targetZ;

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
 W=Wdr*d; 
 
if (W<0)
    disp('desired impedance too high');
    return;
else
    
    disp([num2str(targetZ) ' ohm track width is ' num2str(W) ' mm (RIGHT)']);
    
end
