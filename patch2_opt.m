close all
clear all

DO_GENETIC_ALGORITHM = 0;  % A from scratch optimisation (takes time)
DO_OPT = 0; % A normal minimum seek optimisation
DO_5GHZ = 1; %5 or 2.4GHz design
DO_PATTERNS = 1; % display antenna patterns (time consuming)
DO_ER43 = 1; % Er of 4.3 or 4.8


LOG_FILE = "c:/Users/iwg/test.txt";

[fid,errmsg] = fopen(LOG_FILE,'w');
fprintf(fid,"L, D, Ds, gap, w_100\n");
fclose(fid);
if DO_5GHZ
    f = [5.5e9 5.6e9 5.7e9];
else
    f = [2.4e9 2.44e9 2.48e9]; % 5.5e9 5.6e9 5.7e9];
end
f1 = f(2);
fmin = min(f);
fmax = max(f);
% er = 4.2;
if DO_ER43
    er = 4.3;
else
    er = 4.8;
end

f1_wavelength = 3e8/(sqrt(er)*f1);












if (DO_5GHZ)
    if DO_ER43
        L = 13.318e-3;
        W = 12.46e-3;
        y0 = 5.11e-3; 
%         vertdist = W*1.3;
    else
        L = 12.1297e-3;
        W = 11.5248e-3;
        y0 = 1.6424e-3;
    end
else
    if DO_ER43
        L = 32.086e-3;
        W = 29.96e-3;
        y0 = 4.67e-3; 
    else
        L = 29.7e-3;
        W = 29.6252e-3;
        y0 = 1.6243e-3;   
    end
end


% fixed parameters
h = 1.465e-3;  % pcb thickness
Z0 = 100;
feedD = .5e-3;  % 1mm
%  vertdist = L*3.0;
 norows = 1;

edge = 5e-3;
Nfreq = 21;
freq = linspace(fmin,fmax,Nfreq);

% fixed parameters
% h = 1.55e-3;  % pcb thickness

vertdist = L*2.0;
hordist = L*1.5;
norows = 1;
nocols = 1;

% Z0 = 100;
dd = h*1e3;

d = dielectric('FR4');
d.EpsilonR = er;
d.LossTangent = 0.008;
d.Thickness = h;

 Zx=Z0;
 Er = er;

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
 y = 2*w_m;

 feedD = w_m/2;  % 1mm

Nfreq = 21;
freq = linspace(fmin,fmax,Nfreq);

vals_opt(1) = L;
vals_opt(2) = W;
vals_opt(3) = y0;

% limits for ga
vals_opt_min(1) = 10e-3;
vals_opt_min(2) = 10e-3;
vals_opt_min(3) = 1e-4;



vals_opt_max(1) = 40e-3;
vals_opt_max(2) = 40e-3;
vals_opt_max(3) = 5e-3;

% linear constraints according to ga
% Ds < (L*sqrt(2)- feedD*3 (Ds-L*sqrt(2) < -feedD*3)
% D > L*sqrt(2) or L*sqrt(2) - D < 0
% gap > 5e-3 or -gap < -5e-3
% gap < 15e-3
% Ds > feedD*3 or -Ds < -feedD*3
% A > 6e-3 or -A < -6e-3
% D < lam/4
% A < lam/8
A = [ 1, 0, 0 ;...
      0, 1 , 0 ;...
      -.45 , 0 , 1 ;...
      -1 , 0 , 0 ;...
      0 , -1 , 0 ;...
      0 , 0 , -1 ;...      
     ]; 
%      0, 0, 0, 0, 0; ...
%      0, 0, 0, 0, 0];
 b = [f1_wavelength/1.5 ; f1_wavelength/1.5; -1e-4; -5e-3; -5e-3; -1e-4];

% mesh(blade_dipole, 'MaxEdgeLength',f2_wavelength/5,'MinEdgeLength',f2_wavelength/20)


FitnessFunction = @(vals_opt) optvswr(vals_opt,f,edge,h,er,Z0,LOG_FILE,w_m,feedD,vertdist,norows,hordist,nocols);

options = optimset('Display','iter','TolFun',0.1);
options2 = optimoptions(@ga,'Display','iter','UseParallel', true,'FunctionTolerance',1e-2,'MaxTime',60*2);

if DO_OPT
    if DO_GENETIC_ALGORITHM
        vals_opt = ga(@(vals_opt) FitnessFunction(vals_opt),3,A,b,[],[],[],[],[],options2);
        vals_opt = fminsearch(@(vals_opt) FitnessFunction(vals_opt),vals_opt,options);
    else
        vals_opt = fminsearch(@(vals_opt) FitnessFunction(vals_opt),vals_opt,options);
    end
end
L = vals_opt(1);
W = vals_opt(2);
y0 = vals_opt(3);

vertdist = L*1.5;
hordist = L*2;

patch = arrgeom(L,W,y0,edge,h,er,w_m,feedD,vertdist,norows,hordist,nocols);

figure(3);
impedance(patch,freq);

figure(4);
s_blade1 = sparameters(patch,freq,Z0);
   
imp = s_blade1.Parameters;
imp2 = reshape(imp(1,1,:),[length(freq) 1]);
vv = vswr(imp2);
plot(freq,vv);

if DO_5GHZ
    norows = 8;
else
    norows = 4;
end
nocols = 2;

vertspace = vertdist-W
horspace = hordist - L
patch = arrgeom(L,W,y0,edge,h,er,w_m,feedD,vertdist,norows,hordist,nocols);
if DO_PATTERNS
    figure(5);
    % pattern(patch,f1);
    figure(6)
    patternElevation(patch,f1,90);
    dirb = patternElevation(patch,f1,90);
    [maxer maxi] = max(dirb);
    three = maxer-3;
    new = (dirb>=three);
    elbeam = sum(new)
    figure(7)
    patternElevation(patch,f1);
    dirb = patternElevation(patch,f1);
    [maxer maxi] = max(dirb);
    three = maxer-3;
    new = (dirb>=three);
    azbeam = sum(new)
end


function output=optvswr(vals_opt,f,edge,h,er,Z0,LOG_FILE,w_m,feedD,vertdist,norows,hordist,nocols)
   
    L = vals_opt(1);
    W = vals_opt(2);
    y0 = vals_opt(3);


    fid = fopen(LOG_FILE,'a');
    fprintf(fid,"L=%f mm, W=%f mm, y0=%f mm\n",L*1000,W*1000,y0*1000);
    fclose(fid);

    patch = arrgeom(L,W,y0,edge,h,er,w_m,feedD,vertdist,norows,hordist,nocols);


    
    s_blade1 = sparameters(patch,f,Z0);
   
    imp = s_blade1.Parameters;
    ss = size(imp);
    lngth = prod(ss);
    imp = reshape(imp,[lngth 1]);
    output = mean(abs(sum(vswr(imp))));
    


    fid = fopen(LOG_FILE,'a');
    fprintf(fid,"VSWR=%f\n",output);
    fclose(fid);

end



function anten=arrgeom(L,W,y0,edge,h,er,w_m,feedD,vertdist,norows,hordist,nocols)
    
   
    Wgnd = (W+edge)*2 + vertdist*(norows-1);;       % pcb size
    Lgnd = (L + edge)*2 + hordist*(nocols-1);;
    
%     Lgnd = (L*sqrt(2)+D2+edge)*2;       % pcb size
%     Wgnd = (L + edge)*2 + vertdist*(norows-1);
%     
%     
    
    % make the top pcb structure
    % patch
    patchtemp = antenna.Rectangle('Length',L,'Width',W,...
                               'Center',[0 , 0 ],...
                               'NumPoints', [2,2,2,2]);
                           
    feedpoint = -L/2-L/3;
                           
    patchslot = antenna.Rectangle('Length',L/3 + y0,'Width',w_m*2,...
                               'Center',[-L/2+y0-((L/3+y0)/2) , 0 ],...
                               'NumPoints', [2,2,2,2]);
                           
    patchfeed = antenna.Rectangle('Length',L/3 + y0,'Width',w_m,...
                               'Center',[-L/2+y0-((L/3+y0)/2) , 0 ],...
                               'NumPoints', [2,2,2,2]);                       
                           

    patch = patchtemp-patchslot+patchfeed;
    patchtemp2 = patchtemp-patchslot+patchfeed;
    
    topLayer = patch;
    for (i=2:nocols)
        patchtemp3 = translate(patchtemp2,[hordist,0,0]);
        topLayer = topLayer+patchtemp3;
        patchtemp2 = patchtemp-patchslot+patchfeed;
    end 
    
   
    for (i=2:norows)
        patchtemp3 = translate(patchtemp2,[0,vertdist*(i-1),0]);
        topLayer = topLayer+patchtemp3;
        patchtemp2 = patchtemp-patchslot+patchfeed;
        for (k=2:nocols)
           patchtemp4 = translate(patchtemp2,[hordist*(k-1),vertdist*(i-1),0]); 
           topLayer = topLayer+patchtemp4;
           patchtemp2 = patchtemp-patchslot+patchfeed;
        end
    end
  
    

    boardShape = antenna.Rectangle('Length',Lgnd,'Width',Wgnd,'Center',[hordist*(nocols-1)/2 , vertdist*(norows-1)/2 ]);
                                 

   
    substrate1 = dielectric('Name','FR4','EpsilonR', er, 'Thickness', h);
    
    
    metal = antenna.Rectangle('Length',Lgnd,'Width',Wgnd,...
                               'Center',[hordist*(nocols-1)/2 , vertdist*(norows-1)/2 ],...
                               'NumPoints', [2,2,2,2]);
    % stackup
    anten = pcbStack;
    anten.Name = 'patch';

        anten.BoardThickness = h;
        anten.Layers = {topLayer,substrate1,metal};
        gndlayer = 3;

    anten.BoardShape = boardShape;

    anten.FeedDiameter = feedD;
    anten.FeedLocations = [feedpoint, 0, 1, gndlayer ];
    for k=2:nocols
            anten.FeedLocations = [anten.FeedLocations;feedpoint + (k-1)*hordist, 0, 1, gndlayer ];
    end 
   
    for (i=2:norows)
        anten.FeedLocations = [anten.FeedLocations;feedpoint, (i-1)*vertdist, 1, gndlayer ];
        for k=2:nocols
            anten.FeedLocations = [anten.FeedLocations;feedpoint + (k-1)*hordist, (i-1)*vertdist, 1, gndlayer ];
        end 
    end

    anten.FeedPhase = zeros([1 norows*nocols]);
    h = figure(2);
    show(anten);
    txt1 = ['L: ' num2str(L*1000) ' mm'];
    txt2 = ['W: ' num2str(W*1000) ' mm'];
    txt3 = ['y0: ' num2str(y0*1000) ' mm'];

    txt = {txt1,txt2,txt3};
    annotation(h,'textbox', [0.05, 0.9, 0.1, 0.1], 'String',txt,'FitBoxToText','on'); 

end
