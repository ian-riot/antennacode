% blade dipole

close all
clear all

DO_GENETIC_ALGORITHM = 1

L = 24.34e-3;         % Total half-length of blade dipole
W = 22.63e-3;         % Width of blade dipole
Ld = 10.88e-3;        % Half-length excluding impedance match taper
fw = .1e-3;          % Feed width
g = .1e-3;           % Feed gap

f1 = 2.44e9;
f2 = 5.6e9;

bladeDipole = dipoleBlade('Length', L, 'width', W, 'TaperLength', Ld,...
                           'FeedWidth', fw,'FeedGap', g);

figure(1);
show(bladeDipole);

fmin = 2e9;
fmax = 6e9;
Nfreq = 21;
freq = linspace(fmin,fmax,Nfreq);
s_blade = sparameters(bladeDipole,freq);

f = s_blade.Frequencies;
imp = s_blade.Parameters;
imp = reshape(imp,[length(imp) 1]);


vals_opt(1) = L;
vals_opt(2) = W;
vals_opt(3) = Ld;

% limits for ga
vals_opt_min(1) = L*.8;
vals_opt_min(2) = W*.8;
vals_opt_min(3) = Ld*.8;


vals_opt_max(1) = L*1.5;
vals_opt_max(2) = W*1.5;
vals_opt_max(3) = Ld*1.5;


figure(2)
vswrout = vswr(imp);

plot(f,vswrout);

FitnessFunction = @(vals_opt) optvswr(bladeDipole,vals_opt,freq,f1,f2);


A = [-1, 0, 1;...
      1, 0, -1 ]; 
%      0, 0, 0, 0, 0; ...
%      0, 0, 0, 0, 0];
 b = [-1e-4; 1e-2];

options = optimset('Display','iter','TolFun',0.1);
options2 = optimoptions(@ga,'Display','iter','UseParallel', true,'FunctionTolerance',1e-2,'MaxTime',60*20);

if DO_GENETIC_ALGORITHM
    vals_opt = ga(@(vals_opt) optvswr(bladeDipole,vals_opt,freq,f1,f2),3,A,b,[],[],vals_opt_min,vals_opt_max,[],options2);
    vals_opt = fminsearch(@(vals_opt) FitnessFunction(vals_opt),vals_opt,options);
else
    vals_opt = fminsearch(@(vals_opt) FitnessFunction(vals_opt),vals_opt,options);
end

vals_opt = fminsearch(@(vals_opt) FitnessFunction(vals_opt),vals_opt,options);


figure(3);
show(bladeDipole);

vals_opt*1000

function output = optvswr(dipole,vals_opt,freq,f1,f2)
    dipole.Length = vals_opt(1);
    dipole.Width = vals_opt(2);
    dipole.TaperLength = min([vals_opt(3) vals_opt(1)-vals_opt(1)/100]);
    disp(['L:' num2str(dipole.Length*1000) ' W: ' num2str(dipole.Width*1000) ' taper: ' num2str(dipole.TaperLength*1000) ]);   
    s_blade = sparameters(dipole,freq,100);

    f = s_blade.Frequencies;
    imp = s_blade.Parameters;
    imp = reshape(imp,[length(imp) 1]);
    
    h = figure(2)

    vswrout = vswr(imp);
    plot(f,vswrout);

    
    s_blade1 = sparameters(dipole,[f1 f2],100);
   
    imp = s_blade1.Parameters;
    imp = reshape(imp,[length(imp) 1]);
    output = sum(vswr(imp));
end


