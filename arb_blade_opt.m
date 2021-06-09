% blade dipole

close all
clear all

DO_GENETIC_ALGORITHM = 0

%GeometryFunction = @(anten) geom(L,W,Ld,g,w_100,edge,h,er);

L = 16.34e-3;         % Total half-length of blade dipole
W = 25.63e-3;         % Width of blade dipole
Ld = 12.88e-3;        % Half-length excluding impedance match taper
g = 1e-3;          % Feed gap (apex point)
w_100 = 0.749e-3;    % transmission line width



L = 18.06e-3;         % Total half-length of blade dipole
W = 18.34e-3;         % Width of blade dipole
Ld = 13.47e-3;        % Half-length excluding impedance match taper
g = 3.89e-3;          % Feed gap (apex point)
w_100 = 0.749e-3*2;    % transmission line width

L = 24.1e-3;         % Total half-length of blade dipole
W = 22.43e-3;         % Width of blade dipole
Ld = 16.82e-3;        % Half-length excluding impedance match taper
g = 2.04e-3;          % Feed gap (apex point)
w_100 = 1.58e-3;    % transmission line width

% fixed parameters
h = 1.6e-3;
er = 4.2;
edge = 5e-3;
f1 = 2.44e9;
f2 = 5.6e9;
Z0 = 100;
f2_wavelength = 3e8/(sqrt(er)*f2);
fmin = f1-1e9;
fmax = f2+1e9;
Nfreq = 21;
freq = linspace(fmin,fmax,Nfreq);

vals_opt(1) = L;
vals_opt(2) = W;
vals_opt(3) = Ld;
vals_opt(4) = g;
vals_opt(5) = w_100;

% limits for ga
vals_opt_min(1) = L*.6;
vals_opt_min(2) = W*.6;
vals_opt_min(3) = Ld*.6;
vals_opt_min(4) = w_100/2;
vals_opt_min(5) = w_100*.9;

vals_opt_max(1) = L*1.5;
vals_opt_max(2) = W*1.5;
vals_opt_max(3) = Ld*1.5;
vals_opt_max(4) = 10e-3;
vals_opt_max(5) = w_100*1.1;

% linear constraints according to ga
% L>Ld (L-Ld>=0 or Ld-L<=0) 
% g>w_100 (g-w_100>=0 or w_100-g<=0)
A = [-1, 0, 1, 0, 0;...
     0,  0, 0, -1, 1; ...
     -1 ,0, 1, -0.5, 0.5]; 
%      0, 0, 0, 0, 0; ...
%      0, 0, 0, 0, 0];
 b = [-1e-4; -1e-4; -1e-4];


 blade_dipole = geom(L,W,Ld,g,w_100,edge,h,er);



% mesh(blade_dipole, 'MaxEdgeLength',f2_wavelength/5,'MinEdgeLength',f2_wavelength/20)

f = [2.4e9 2.44e9 2.48e9 5.5e9 5.6e9 5.7e9];

FitnessFunction = @(vals_opt) optvswr(vals_opt,f,edge,h,er,Z0);

options = optimset('Display','iter','TolFun',0.1);
options2 = optimoptions(@ga,'Display','iter','UseParallel', true,'FunctionTolerance',1e-2,'MaxTime',60*20);

if DO_GENETIC_ALGORITHM
    vals_opt = ga(@(vals_opt) optvswr(vals_opt,f,edge,h,er,Z0),5,A,b,[],[],vals_opt_min,vals_opt_max,[],options2);
    vals_opt = fminsearch(@(vals_opt) FitnessFunction(vals_opt),vals_opt,options);
else
    vals_opt = fminsearch(@(vals_opt) FitnessFunction(vals_opt),vals_opt,options);
end


L = vals_opt(1);
W = vals_opt(2);
Ld = vals_opt(3);
g = vals_opt(4);
w_100 = vals_opt(5);

blade_dipole = geom(L,W,Ld,g,w_100,edge,h,er);

figure(3);
impedance(blade_dipole,freq);

figure(4);
s_blade1 = sparameters(blade_dipole,freq,Z0);
   
imp = s_blade1.Parameters;
imp = reshape(imp,[length(imp) 1]);
vv = vswr(imp);
plot(freq,vv);

function output = optimp(vals_opt,f,edge,h,er)
    outofbounds = 0;
    L = vals_opt(1);
    W = vals_opt(2);
    Ld = vals_opt(3);
    if (Ld > L-0.1e-4)
        Ld = L-0.1e-4;
        outofbounds = 1;
    end
    g = vals_opt(4);

    w_100 = vals_opt(5);
    if (g<w_100)
        g = w_100;
        outofbounds=1;
    end    
    disp(['L:' num2str(L*1000) ' W: ' num2str(W*1000) ' taper: ' num2str(Ld*1000)  ' g: ' num2str(g*1000)  ' w: ' num2str(w_100*1000)  ]);   
   
    blade_dipole = geom(L,W,Ld,g,w_100,edge,h,er);
%     f2_wavelength = 3e8/(sqrt(er)*f2);
%     mesh(blade_dipole, 'MaxEdgeLength',f2_wavelength/5,'MinEdgeLength',f2_wavelength/20);
%     twof = [f1 f2];
    imp = impedance(blade_dipole,f);
    output = sum(abs(real(imp)-Z0));
    output = output + sum(abs(imag(imp)-0));

    output = output + outofbounds*output;

end

function output = optvswr(vals_opt,f,edge,h,er,Z0)
    outofbounds = 0;
    L = vals_opt(1);
    W = vals_opt(2);
    Ld = vals_opt(3);
%     if (Ld > L-0.1e-4)
%         Ld = L-0.1e-4;
%         outofbounds = 1;
%     end
    g = vals_opt(4);

    w_100 = vals_opt(5);
%     if (g<w_100)
%         g = w_100;
%         outofbounds=1;
%     end    
%      disp(['L:' num2str(L*1000) ' W: ' num2str(W*1000) ' taper: ' num2str(Ld*1000)  ' g: ' num2str(g*1000)  ' w: ' num2str(w_100*1000)  ]);   
   
    blade_dipole = geom(L,W,Ld,g,w_100,edge,h,er);
%     f2_wavelength = 3e8/(sqrt(er)*f2);
%     mesh(blade_dipole, 'MaxEdgeLength',f2_wavelength/5,'MinEdgeLength',f2_wavelength/20);
%     twof = [f1 f2];
%     s_blade = sparameters(dipole,freq);
% 
%     f = s_blade.Frequencies;
%     imp = s_blade.Parameters;
%     imp = reshape(imp,[length(imp) 1]);
%     
%     h = figure(2)
% 
%     vswrout = vswr(imp);
%     plot(f,vswrout);

    
    s_blade1 = sparameters(blade_dipole,f,Z0);
   
    imp = s_blade1.Parameters;
    imp = reshape(imp,[length(imp) 1]);
    output = abs(sum(vswr(imp)));

    output = output + outofbounds*output;

end



function anten=geom(L,W,Ld,g,w_100,edge,h,er)
    Lgnd = L+L+g+edge;       % pcb size
    Wgnd = W + edge;
    g2 = g/2;

    % make the top pcb structure
    % patch
    patch1 = antenna.Rectangle('Length',Ld,'Width',W,...
                               'Center',[-g2 - (L-Ld) - Ld/2 , 0 ],...
                               'NumPoints', [2,2,2,2]);
                           
    patch2 = antenna.Rectangle('Length',Ld,'Width',W,...
                               'Center',[+g2 + (L-Ld) + Ld/2 , 0 ],...
                               'NumPoints', [2,2,2,2]);

    % taper
    taper1 = antenna.Polygon('Vertices',[ -g2-(L-Ld) , W/2 , 0 ; ...
                                          -g2-(L-Ld) , -W/2 , 0 ; ...
                                          -g2 , 0 , 0]);
                                      
    taper2 = antenna.Polygon('Vertices',[ +g2+(L-Ld) , W/2 , 0 ; ...
                                          +g2+(L-Ld) , -W/2 , 0 ; ...
                                          +g2 , 0 , 0]);                                  

    % feed 1
    feed1 = antenna.Rectangle('Length',g2+(L-Ld)-w_100/2,'Width',w_100,...
                               'Center',[(-g2 - (L-Ld) - w_100/2)/2 , 0 ],...
                               'NumPoints', [2,2,2,2]);
    
    feed2 = antenna.Rectangle('Length',g2+(L-Ld)-w_100/2,'Width',w_100,...
                               'Center',[(+g2 + (L-Ld) + w_100/2)/2 , 0 ],...
                               'NumPoints', [2,2,2,2]);

                           
    % chamfer                       
    chamfer1 = antenna.Polygon('Vertices',[-w_100/2 , w_100/2 , 0 ; ...
                                           -w_100/2 , -w_100/2 , 0 ; ...
                                           w_100/2 , -w_100/2 , 0]);
                                       
                                       
    chamfer2 = antenna.Polygon('Vertices',[+w_100/2 , w_100/2 , 0 ; ...
                                           +w_100/2 , -w_100/2 , 0 ; ...
                                           -w_100/2 , -w_100/2 , 0]);

    feed12 = antenna.Rectangle('Length',w_100,'Width',Wgnd/2 - w_100/2,...
                               'Center',[0 , (-Wgnd/2 - w_100/2)/2  ],...
                               'NumPoints', [2,2,2,2]);   
                           
    feed22 = antenna.Rectangle('Length',w_100,'Width',Wgnd/2 - w_100/2,...
                               'Center',[0 , (-Wgnd/2 - w_100/2)/2  ],...
                               'NumPoints', [2,2,2,2]);  


     bottomLayer = patch2+taper2+feed2+chamfer2+feed22;
    
    % bottomlayer is mirror around y=0
%     rotateY(bottomLayer,180);
    topLayer = patch1+taper1+feed1+chamfer1+feed12;
    boardShape = antenna.Rectangle('Length',Lgnd,'Width',Wgnd);
    figure(1);
    hold on;
    plot(topLayer)
    plot(bottomLayer)
    plot(boardShape)
    grid on
    hold off;

    % place a dielectric
    substrate = dielectric('Name','FR4','EpsilonR', er, 'Thickness', h);

    % stackup
    anten = pcbStack;
    anten.Name = 'blade dipole';
    anten.BoardThickness = h;
    anten.BoardShape = boardShape;
    anten.Layers = {topLayer,substrate,bottomLayer};
    anten.FeedLocations = [0, -(Wgnd/2), 1, 3];
    anten.FeedDiameter = w_100/2;
    figure(2);
    show(anten);

end


