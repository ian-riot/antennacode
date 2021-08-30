% cavity backed diagonal micropatch

close all
clear all

DO_GENETIC_ALGORITHM = 0;
DO_OPT = 1;


LOG_FILE = "c:/Users/iwg/test.txt";

[fid,errmsg] = fopen(LOG_FILE,'w');
fprintf(fid,"L, D, Ds, gap, w_100\n");
fclose(fid);

f = [2.4e9 2.44e9 2.48e9]; % 5.5e9 5.6e9 5.7e9];
f = [5.1e9 5.6e9 5.7e9];
f1 = f(2);
fmin = min(f);
fmax = max(f);
% er = 4.2;
er = 4.0;
f1_wavelength = 3e8/(sqrt(er)*f1);


% patch
L = 30e-3
W = 30e-3
y0 = 5e-3

feedD = 1e-4;

% fixed parameters
h = 1.55e-3;  % pcb thickness

vertdist = L*2.0;
hordist = L*1.5;
norows = 1;
nocols = 1;

Z0 = 100;
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
% A = [-sqrt(2), 0, 1, 0;...
%       sqrt(2), -1 , 0 , 0;...
%       0 , 0 , 0 , -1;...
%       0 , 0 , 0 , 1;...
%       0, 0, -1, 0;...
%       -1, 0, 0, 0;...
%       0, 1, 0, 0;...
%       1, 0, 0, 0;...
%      ]; 
%      0, 0, 0, 0, 0; ...
%      0, 0, 0, 0, 0];
%  b = [-feedD*6 ; -1e-3; -5e-3; 15e-3; -feedD*6; -2e-3; f1_wavelength/2 ; f1_wavelength/4];

FitnessFunction = @(vals_opt) optvswr(vals_opt,f,edge,h,er,Z0,LOG_FILE,w_m,vertdist,norows,hordist,nocols);

options = optimset('Display','iter','TolFun',0.5);
options2 = optimoptions(@ga,'Display','iter','UseParallel', false,'FunctionTolerance',1e-2,'MaxTime',60*60);
if DO_OPT
    if DO_GENETIC_ALGORITHM
        vals_opt = ga(@(vals_opt) optvswr(vals_opt,f,edge,h,er,Z0,LOG_FILE,w_m,vertdist,norows,hordist,nocols),3,[],[],[],[],vals_opt_min,vals_opt_max,[],options2);
%         vals_opt = fminsearch(@(vals_opt) FitnessFunction(vals_opt),vals_opt,options);
    else
        vals_opt = fminsearch(@(vals_opt) FitnessFunction(vals_opt),vals_opt,options);
    end
end
L = vals_opt(1);
W = vals_opt(2);
y0 = vals_opt(3);



patch = arrgeom(L,W,y0,edge,h,er,w_m,vertdist,norows,hordist,nocols);

figure(3);
impedance(patch,freq);

figure(4);
s_blade1 = sparameters(patch,freq,Z0);
   
imp = s_blade1.Parameters;
imp2 = reshape(imp(1,1,:),[length(freq) 1]);
vv = vswr(imp2);
plot(freq,vv);

norows = 8;
nocols = 2;
patch = arrgeom(L,W,y0,edge,h,er,w_m,vertdist,norows,hordist,nocols);
figure(5);
pattern(micropatch,f1);
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




function output = optvswr(vals_opt,f,edge,h,er,Z0,LOG_FILE,w_m,vertdist,norows,hordist,nocols)
   
    L = vals_opt(1);
    W = vals_opt(2);
    y0 = vals_opt(3);


    fid = fopen(LOG_FILE,'a');
    fprintf(fid,"L=%f mm, W=%f mm, y0=%f mm\n",L*1000,W*1000,y0*1000);
    fclose(fid);

    patch = arrgeom(L,W,y0,edge,h,er,w_m,vertdist,norows,hordist,nocols);


    
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



function anten=arrgeom(L,W,y0,gap,edge,h,er,w_m,vertdist,norows,hordist,nocols)
    
    D2 = D/2;
    Lgnd = (W+edge)*2;       % pcb size
    Wgnd = (L + edge)*2;
    

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
    
    topLayer = patch;

    boardShape = antenna.Rectangle('Length',Lgnd,'Width',Wgnd,'Center',[0 , 0 ]);
                                 

   
    substrate1 = dielectric('Name','FR4','EpsilonR', er, 'Thickness', h);
    
    
    metal = antenna.Rectangle('Length',Lgnd,'Width',Wgnd,...
                               'Center',[0 , 0  ],...
                               'NumPoints', [2,2,2,2]);
    % stackup
    anten = pcbStack;
    anten.Name = 'patch';

        anten.BoardThickness = h+h+gap;
        anten.Layers = {topLayer,substrate1,metal};
        gndlayer = 3;

    anten.BoardShape = boardShape;

    anten.FeedDiameter = 1e-4;
    anten.FeedLocations = [feedpoint, 0, 1, gndlayer ];
   
    anten.ViaDiameter = feedD;


    anten.FeedPhase = 0;
    h = figure(2);
    show(anten);
    txt1 = ['L: ' num2str(L*1000) ' mm'];
    txt2 = ['W: ' num2str(W*1000) ' mm'];
    txt3 = ['y0: ' num2str(y0*1000) ' mm'];

    txt = {txt1,txt2,txt3};
    annotation(h,'textbox', [0.05, 0.9, 0.1, 0.1], 'String',txt,'FitBoxToText','on'); 

end
