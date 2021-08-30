% cavity backed diagonal micropatch

close all
clear all

DO_GENETIC_ALGORITHM = 0;
DO_OPT = 0;
DO_5G = 1;
DO_ANT = 0;
DO_FULL_POP = 1;

LOG_FILE = "c:/Users/iwg/test.txt";

[fid,errmsg] = fopen(LOG_FILE,'w');
fprintf(fid,"L, D, Ds, gap, w_100\n");
fclose(fid);


if DO_5G
    f = [5.5e9 5.6e9 5.7e9];
else
    f = [2.4e9 2.44e9 2.48e9]; 
end
f1 = f(2);
fmin = min(f);
fmax = max(f);
% er = 4.2;
er = 4.0;
f1_wavelength = 3e8/(sqrt(er)*f1);
f1_wavelength = 3e8/(f1);

% micropatch

L = f1_wavelength/10;
D = f1_wavelength/4;
Ds = L/10;
gap =10e-3;  % height above reflector


% fixed parameters
h = 1.55e-3;  % pcb thickness
Z0 = 100;
feedD = .5e-3;  % 1mm
 vertdist = L*2.2;
 norows = 1;
 
    Lb = 3.7845e-3;
%     Db = 21.7111e-3;    
    Db = 25.7111e-3; % made larger for better azimuth beamwidth
    Dsb = 4.7854e-3;
    gapb = 9.3632e-3; 
    
    La = 11.1223e-3;
    Da = 32.5215e-3;
    Dsa = 7.7287e-3;
    gapa = 15.5208e-3;

if DO_5G
    L = Lb;
    D = Db;
    Ds = Dsb;
    gap = gapb;
    vertdist = 20.e-3;
else
    L = La;
    D = Da;
    Ds = Dsa;
    gap = gapa;
    vertdist = 20.e-3;
end




edge = 5e-3;
Nfreq = 21;
freq = linspace(fmin,fmax,Nfreq);

vals_opt(1) = L;
vals_opt(2) = D;
vals_opt(3) = Ds;
vals_opt(4) = gap;




% linear constraints according to ga
% Ds < (L*sqrt(2)- feedD*3 (Ds-L*sqrt(2) < -feedD*3)
% D > L*sqrt(2) or L*sqrt(2) - D < 0
% gap > 5e-3 or -gap < -5e-3
% gap < 15e-3
% Ds > feedD*3 or -Ds < -feedD*3
% A > 6e-3 or -A < -6e-3
% D < lam/4
% A < lam/8
A = [-sqrt(2), 0, 1, 0;...
      sqrt(2), -1 , 0 , 0;...
      0 , 0 , 0 , -1;...
      0 , 0 , 0 , 1;...
      0, 0, -1, 0;...
      -1, 0, 0, 0;...
      0, 1, 0, 0;...
      1, 0, 0, 0;...
     ]; 

 b = [-feedD*6 ; -1e-3; -1e-3; 15e-3; -feedD*8; -2e-3; f1_wavelength/2 ; f1_wavelength/4];




FitnessFunction = @(vals_opt) optvswr(vals_opt,f,edge,h,er,Z0,LOG_FILE,feedD,vertdist,norows);

options = optimset('Display','iter','TolFun',0.5);
options2 = optimoptions(@ga,'Display','iter','UseParallel', true,'FunctionTolerance',1e-2,'MaxTime',60*60);
if DO_OPT
    if DO_GENETIC_ALGORITHM
        vals_opt = ga(@(vals_opt) optvswr(vals_opt,f,edge,h,er,Z0,LOG_FILE,feedD,vertdist,norows),4,A,b,[],[],[],[],[],options2);
%         vals_opt = fminsearch(@(vals_opt) FitnessFunction(vals_opt),vals_opt,options);
    else
        vals_opt = fminsearch(@(vals_opt) FitnessFunction(vals_opt),vals_opt,options);
    end
end
L = vals_opt(1);
D = vals_opt(2);
Ds = vals_opt(3);
gap = vals_opt(4);


micropatch = arrgeom(L,D,Ds,gap,edge,h,er,feedD,vertdist,norows);

figure(3);
impedance(micropatch,freq);

figure(4);
s_blade1 = sparameters(micropatch,freq,Z0);
   
imp = s_blade1.Parameters;
imp2 = reshape(imp(1,1,:),[length(freq) 1]);
vv = vswr(imp2);
plot(freq,vv);

if DO_ANT
    norows = 8;
    micropatch = arrgeom(L,D,Ds,gap,edge,h,er,feedD,vertdist,norows);
    figure(5);
    pattern(micropatch,f1);
    figure(6)
    patternElevation(micropatch,f1,90);
    dirb = patternElevation(micropatch,f1,90);
    [maxer maxi] = max(dirb);
    three = maxer-3;
    new = (dirb>=three);
    elbeam = sum(new)
    figure(7)
    patternElevation(micropatch,f1);
    dirb = patternElevation(micropatch,f1);
    [maxer maxi] = max(dirb);
    three = maxer-3;
    new = (dirb>=three);
    azbeam = sum(new)
end
if DO_FULL_POP
    norows = 8;
%     micropatch = arrgeomdualside(La,Da,Dsa,gapa,Lb,Db,Dsb,gapb,edge,h,er,feedD,vertdist,norows,75e-3);
    micropatch = arrgeomdual(La,Da,Dsa,gapa,Lb,Db,Dsb,gapb,edge,h,er,feedD,vertdist,norows);
    figure(5);
    pattern(micropatch,f1);
    figure(6)
    patternElevation(micropatch,f1,90);
    dirb = patternElevation(micropatch,f1,90);
    [maxer maxi] = max(dirb);
    three = maxer-3;
    new = (dirb>=three);
    elbeam = sum(new)
    figure(7)
    patternElevation(micropatch,f1);
    dirb = patternElevation(micropatch,f1);
    [maxer maxi] = max(dirb);
    three = maxer-3;
    new = (dirb>=three);
    azbeam = sum(new)

end

% make a linear array

% linArr = array(micropatch,'linear','NumElements',norows,'ElementSpacing',vertdist);
% figure(6);
% show(linArr)
% vertdist = L*2.1;
%  norows = 10;
% arr = arrgeom(L,D,Ds,gap,edge,h,er,feedD,vertdist,norows);




function output = optvswr(vals_opt,f,edge,h,er,Z0,LOG_FILE,feedD,vertdist,norows)
    outofbounds = 0;
    L = vals_opt(1);
    D = vals_opt(2);
    Ds = vals_opt(3);
%     if (Ld > L-0.1e-4)
%         Ld = L-0.1e-4;
%         outofbounds = 1;
%     end
    gap = vals_opt(4);

    fid = fopen(LOG_FILE,'a');
    fprintf(fid,"L=%f mm, D=%f mm, Ds=%f mm, gap=%f mm\n",L*1000,D*1000,Ds*1000,gap*1000);
    fclose(fid);

    micropatch = arrgeom(L,D,Ds,gap,edge,h,er,feedD,vertdist,norows);


    
    s_blade1 = sparameters(micropatch,f,Z0);
   
    imp = s_blade1.Parameters;
    ss = size(imp);
    lngth = prod(ss);
    imp = reshape(imp,[lngth 1]);
    output1 = mean(abs(sum(vswr(imp))));
    
    dirb = patternElevation(micropatch,f(2));
    [maxer maxi] = max(dirb);
    three = maxer-3;
    new = (dirb>=three);
    % fill in the center if need be
    sw = 1;
    for i=1:length(new)
        if sw == 1
            if new(i) == 0
                new(i) = 1;
            else
                sw = 0;
            end

        end
    end
    sw = 1;
    for i=length(new):-1:1
        if sw == 1
            if new(i) == 0
                new(i) = 1;
            else
                sw = 0;
            end

        end
    end    
    azbeam = sum(new);
    
%     output = output1 + abs(azbeam-90)/2;
%     output = output + (20-maxer);
    output = output1;

    fid = fopen(LOG_FILE,'a');
    fprintf(fid,"VSWR=%f\n",output);
    fclose(fid);

end



function anten=arrgeom(L,D,Ds,gap,edge,h,er,feedD,vertdist,norows)
    DO_TOP_FR4 = 0;
    DO_BOTTOM_FR4 = 1;
    D2 = D/2;
    Lgnd = (L*sqrt(2)+D2+edge)*2;       % pcb size
    Wgnd = (L + edge)*2 + vertdist*(norows-1);
    

    % make the top pcb structure
    % patch
    patchtemp1 = antenna.Rectangle('Length',L,'Width',L,...
                               'Center',[0 , 0 ],...
                               'NumPoints', [2,2,2,2]);
                           
    patchtemp2 = antenna.Rectangle('Length',L,'Width',L,...
                               'Center',[0 , 0 ],...
                               'NumPoints', [2,2,2,2]);
                           
                           
    patch1 = rotateZ(patchtemp1,45);
    patch2 = rotateZ(patchtemp2,45);
    patch1 = translate(patch1,[-D2,0,0]);
    patch2 = translate(patch2,[D2,0,0]);
    
    topLayer = patch1+patch2;
    for (i=2:norows)
        patch1 = translate(patch1,[0,vertdist,0]);
        patch2 = translate(patch2,[0,vertdist,0]);
        topLayer = topLayer+patch1+patch2;
    end

    boardShape = antenna.Rectangle('Length',Lgnd,'Width',Wgnd,'Center',[0 , vertdist*(norows-1)/2  ]);
                                 

    % place a dielectric
    substrate1 = dielectric('Name','FR4','EpsilonR', er, 'Thickness', h);
    substrate2 = dielectric('Name','FR4','EpsilonR', er, 'Thickness', h);
    air =  dielectric('Name','Air','Thickness', gap);
    metal = antenna.Rectangle('Length',Lgnd,'Width',Wgnd,...
                               'Center',[0 , vertdist*(norows-1)/2  ],...
                               'NumPoints', [2,2,2,2]);
    % stackup
    anten = pcbStack;
    anten.Name = 'micropatch';
    if (DO_BOTTOM_FR4)
        anten.BoardThickness = h+gap;
        anten.Layers = {topLayer,air,substrate2,metal};
        gndlayer = 4;
    elseif (DO_TOP_FR4)
        anten.BoardThickness = h+gap;
        anten.Layers = {topLayer,substrate1,air,metal};
        gndlayer = 4;       
    else
        anten.BoardThickness = gap;
        anten.Layers = {topLayer,air,metal};
        gndlayer = 3;
    end
    anten.BoardShape = boardShape;
%     anten.FeedViaModel = 'square';
    anten.FeedDiameter = feedD;
    anten.FeedLocations = [-D2-(L/sqrt(2))+feedD*2.5, 0, 1, gndlayer ];
    anten.FeedLocations = [anten.FeedLocations;D2+(L/sqrt(2))-feedD*2.5, 0, 1, gndlayer ];
    for (i=2:norows)
        anten.FeedLocations = [anten.FeedLocations;-D2-(L/sqrt(2))+feedD*2, (i-1)*vertdist, 1, gndlayer ];
        anten.FeedLocations = [anten.FeedLocations;D2+(L/sqrt(2))-feedD*2, (i-1)*vertdist, 1, gndlayer ];
    end
    anten.ViaDiameter = feedD;

    anten.ViaLocations =[ -D2-(L/sqrt(2))+Ds , 0 , 1 , gndlayer ];
    anten.ViaLocations =[anten.ViaLocations; D2+(L/sqrt(2))-Ds , 0 , 1 , gndlayer ];
    for (i=2:norows)
        anten.ViaLocations =[anten.ViaLocations;  -D2-(L/sqrt(2))+Ds , (i-1)*vertdist , 1 , gndlayer ];
        anten.ViaLocations =[anten.ViaLocations; D2+(L/sqrt(2))-Ds , (i-1)*vertdist , 1 , gndlayer ];
    end
   aa = [0 180];
    for (i=2:norows)
        aa =[aa 0 ];
        aa =[aa 180 ];
    end
    anten.FeedPhase = aa;
    h = figure(2);
    show(anten);
    txt1 = ['L: ' num2str(L*1000) ' mm'];
    txt2 = ['D: ' num2str(D*1000) ' mm'];
    txt3 = ['Ds: ' num2str(Ds*1000) ' mm'];
    txt4 = ['gap: ' num2str(gap*1000) ' mm'];
    % txt5 = ['50 ohm width: ' w_m*1000) ' mm'];
    txt = {txt1,txt2,txt3,txt4};
    annotation(h,'textbox', [0.05, 0.9, 0.1, 0.1], 'String',txt,'FitBoxToText','on'); 

end


function anten=arrgeomdual(La,Da,Dsa,gapa,Lb,Db,Dsb,gapb,edge,h,er,feedD,vertdist,norows)
    DO_FR4 = 0;
    D2a = Da/2;
    D2b = Db/2;
    Lgnd = (La*sqrt(2)+D2a+edge)*2;       % pcb size
    Wgnd = (La + edge)*2 + vertdist*(norows-1);
    

    % make the top pcb structure
    % patch
    % toplayer
    patchtemp1 = antenna.Rectangle('Length',La,'Width',La,...
                               'Center',[0 , 0 ],...
                               'NumPoints', [2,2,2,2]);
                           
    patchtemp2 = antenna.Rectangle('Length',La,'Width',La,...
                               'Center',[0 , 0 ],...
                               'NumPoints', [2,2,2,2]);
                           
                           
    patch1 = rotateZ(patchtemp1,45);
    patch2 = rotateZ(patchtemp2,45);
    patch1 = translate(patch1,[-D2a,0,0]);
    patch2 = translate(patch2,[D2a,0,0]);
    
    topLayer = patch1+patch2;
    for (i=2:norows)
        patch1 = translate(patch1,[0,vertdist,0]);
        patch2 = translate(patch2,[0,vertdist,0]);
        topLayer = topLayer+patch1+patch2;
    end
    
    % midlayer
    patchtemp1 = antenna.Rectangle('Length',Lb,'Width',Lb,...
                               'Center',[0 , 0 ],...
                               'NumPoints', [2,2,2,2]);
                           
    patchtemp2 = antenna.Rectangle('Length',Lb,'Width',Lb,...
                               'Center',[0 , 0 ],...
                               'NumPoints', [2,2,2,2]);
                           
                           
    patch1 = rotateZ(patchtemp1,45);
    patch2 = rotateZ(patchtemp2,45);
    patch1 = translate(patch1,[-D2b,vertdist/2,0]);
    patch2 = translate(patch2,[D2b,vertdist/2,0]);
    
    midLayer = patch1+patch2;
    for (i=2:norows)
        patch1 = translate(patch1,[0,vertdist,0]);
        patch2 = translate(patch2,[0,vertdist,0]);
        midLayer = midLayer+patch1+patch2;
    end    

    boardShape = antenna.Rectangle('Length',Lgnd,'Width',Wgnd,'Center',[0 , vertdist*(norows-1)/2  ]);
                                 

    % place a dielectric
    substrate1 = dielectric('Name','FR4','EpsilonR', er, 'Thickness', h);
    substrate2 = dielectric('Name','FR4','EpsilonR', er, 'Thickness', h);
    aira =  dielectric('Name','Air','Thickness', gapa-gapb);
    airb =  dielectric('Name','Air','Thickness', gapb);
    metal = antenna.Rectangle('Length',Lgnd,'Width',Wgnd,...
                               'Center',[0 , vertdist*(norows-1)/2  ],...
                               'NumPoints', [2,2,2,2]);
    % stackup
    anten = pcbStack;
    anten.Name = 'micropatch';
    if (DO_FR4)
        anten.BoardThickness = h+h+gap;
        anten.Layers = {topLayer,substrate1,air,substrate2,metal};
        gndlayer = 5;
    else
        anten.BoardThickness = gapa;
        anten.Layers = {topLayer,aira,midLayer,airb,metal};
        gndlayer = 5;
    end
    anten.BoardShape = boardShape;
%     anten.FeedViaModel = 'square';
    anten.FeedDiameter = feedD;
    anten.FeedLocations = [-D2a-(La/sqrt(2))+feedD*2.5, 0, 1, gndlayer ];
    anten.FeedLocations = [anten.FeedLocations;D2a+(La/sqrt(2))-feedD*2.5, 0, 1, gndlayer ];
    for (i=2:norows)
        anten.FeedLocations = [anten.FeedLocations;-D2a-(La/sqrt(2))+feedD*2, (i-1)*vertdist, 1, gndlayer ];
        anten.FeedLocations = [anten.FeedLocations;D2a+(La/sqrt(2))-feedD*2, (i-1)*vertdist, 1, gndlayer ];
    end
    
    anten.FeedLocations = [anten.FeedLocations;-D2b-(Lb/sqrt(2))+feedD*2.5, vertdist/2, 3, gndlayer ];
    anten.FeedLocations = [anten.FeedLocations;D2b+(Lb/sqrt(2))-feedD*2.5, vertdist/2, 3, gndlayer ];
    for (i=2:norows)
        anten.FeedLocations = [anten.FeedLocations;-D2b-(Lb/sqrt(2))+feedD*2, (i-1)*vertdist + vertdist/2, 3, gndlayer ];
        anten.FeedLocations = [anten.FeedLocations;D2b+(Lb/sqrt(2))-feedD*2, (i-1)*vertdist + vertdist/2, 3, gndlayer ];
    end
    anten.ViaDiameter = feedD;

    anten.ViaLocations =[ -D2a-(La/sqrt(2))+Dsa , 0 , 1 , gndlayer ];
    anten.ViaLocations =[anten.ViaLocations; D2a+(La/sqrt(2))-Dsa , 0 , 1 , gndlayer ];
    for (i=2:norows)
        anten.ViaLocations =[anten.ViaLocations;  -D2a-(La/sqrt(2))+Dsa , (i-1)*vertdist , 1 , gndlayer ];
        anten.ViaLocations =[anten.ViaLocations; D2a+(La/sqrt(2))-Dsa , (i-1)*vertdist , 1 , gndlayer ];
    end
    anten.ViaLocations =[anten.ViaLocations; -D2b-(Lb/sqrt(2))+Dsb , vertdist/2 , 3 , gndlayer ];
    anten.ViaLocations =[anten.ViaLocations; D2b+(Lb/sqrt(2))-Dsb , vertdist/2 , 3 , gndlayer ];
    for (i=2:norows)
        anten.ViaLocations =[anten.ViaLocations;  -D2b-(Lb/sqrt(2))+Dsb , (i-1)*vertdist + vertdist/2 , 3 , gndlayer ];
        anten.ViaLocations =[anten.ViaLocations; D2b+(Lb/sqrt(2))-Dsb, (i-1)*vertdist + vertdist/2 , 3 , gndlayer ];
    end
   aa = [0 180];
    for (i=2:norows)
        aa =[aa 0 ];
        aa =[aa 180 ];
    end
    aa = [aa 0 180];
    for (i=2:norows)
        aa =[aa 0 ];
        aa =[aa 180 ];
    end
    anten.FeedPhase = aa;
    h = figure(2);
    show(anten);
%     txt1 = ['L: ' num2str(L*1000) ' mm'];
%     txt2 = ['D: ' num2str(D*1000) ' mm'];
%     txt3 = ['Ds: ' num2str(Ds*1000) ' mm'];
%     txt4 = ['gap: ' num2str(gap*1000) ' mm'];
%     % txt5 = ['50 ohm width: ' w_m*1000) ' mm'];
%     txt = {txt1,txt2,txt3,txt4};
%     annotation(h,'textbox', [0.05, 0.9, 0.1, 0.1], 'String',txt,'FitBoxToText','on'); 

end

function anten=arrgeomdualside(La,Da,Dsa,gapa,Lb,Db,Dsb,gapb,edge,h,er,feedD,vertdist,norows,hordist)
    DO_FR4 = 0;
    D2a = Da/2;
    D2b = Db/2;
    Lgnd = (La*sqrt(2)+D2a+edge)*2 + hordist;       % pcb size
    Wgnd = (La + edge)*2 + vertdist*(norows-1);
    

    % make the top pcb structure
    % patch
    % toplayer
    patchtemp1 = antenna.Rectangle('Length',La,'Width',La,...
                               'Center',[0 , 0 ],...
                               'NumPoints', [2,2,2,2]);
                           
    patchtemp2 = antenna.Rectangle('Length',La,'Width',La,...
                               'Center',[0 , 0 ],...
                               'NumPoints', [2,2,2,2]);
                           
                           
    patch1 = rotateZ(patchtemp1,45);
    patch2 = rotateZ(patchtemp2,45);
    patch1 = translate(patch1,[-D2a,0,0]);
    patch2 = translate(patch2,[D2a,0,0]);
    
    topLayer = patch1+patch2;
    for (i=2:norows)
        patch1 = translate(patch1,[0,vertdist,0]);
        patch2 = translate(patch2,[0,vertdist,0]);
        topLayer = topLayer+patch1+patch2;
    end
    
    % midlayer
    patchtemp1 = antenna.Rectangle('Length',Lb,'Width',Lb,...
                               'Center',[0 , 0 ],...
                               'NumPoints', [2,2,2,2]);
                           
    patchtemp2 = antenna.Rectangle('Length',Lb,'Width',Lb,...
                               'Center',[0 , 0 ],...
                               'NumPoints', [2,2,2,2]);
                           
                           
    patch1 = rotateZ(patchtemp1,45);
    patch2 = rotateZ(patchtemp2,45);
    patch1 = translate(patch1,[-D2b+hordist,0,0]);
    patch2 = translate(patch2,[D2b+hordist,0,0]);
    
    midLayer = patch1+patch2;
    for (i=2:norows)
        patch1 = translate(patch1,[0,vertdist,0]);
        patch2 = translate(patch2,[0,vertdist,0]);
        midLayer = midLayer+patch1+patch2;
    end    

    boardShape = antenna.Rectangle('Length',Lgnd,'Width',Wgnd,'Center',[hordist/2 , vertdist*(norows-1)/2  ]);
                                 

    % place a dielectric
    substrate1 = dielectric('Name','FR4','EpsilonR', er, 'Thickness', h);
    substrate2 = dielectric('Name','FR4','EpsilonR', er, 'Thickness', h);
    aira =  dielectric('Name','Air','Thickness', gapa-gapb);
    airb =  dielectric('Name','Air','Thickness', gapb);
    metal = antenna.Rectangle('Length',Lgnd,'Width',Wgnd,...
                               'Center',[hordist/2 , vertdist*(norows-1)/2  ],...
                               'NumPoints', [2,2,2,2]);
    % stackup
    anten = pcbStack;
    anten.Name = 'micropatch';
    if (DO_FR4)
        anten.BoardThickness = h+h+gap;
        anten.Layers = {topLayer,substrate1,air,substrate2,metal};
        gndlayer = 5;
    else
        anten.BoardThickness = gapa;
        anten.Layers = {topLayer,aira,midLayer,airb,metal};
        gndlayer = 5;
    end
    anten.BoardShape = boardShape;
%     anten.FeedViaModel = 'square';
    anten.FeedDiameter = feedD;
    anten.FeedLocations = [-D2a-(La/sqrt(2))+feedD*2.5, 0, 1, gndlayer ];
    anten.FeedLocations = [anten.FeedLocations;D2a+(La/sqrt(2))-feedD*2.5, 0, 1, gndlayer ];
    for (i=2:norows)
        anten.FeedLocations = [anten.FeedLocations;-D2a-(La/sqrt(2))+feedD*2, (i-1)*vertdist, 1, gndlayer ];
        anten.FeedLocations = [anten.FeedLocations;D2a+(La/sqrt(2))-feedD*2, (i-1)*vertdist, 1, gndlayer ];
    end
    
    anten.FeedLocations = [anten.FeedLocations;-D2b-(Lb/sqrt(2))+feedD*2.5 + hordist, 0, 3, gndlayer ];
    anten.FeedLocations = [anten.FeedLocations;D2b+(Lb/sqrt(2))-feedD*2.5 + hordist, 0, 3, gndlayer ];
    for (i=2:norows)
        anten.FeedLocations = [anten.FeedLocations;-D2b-(Lb/sqrt(2))+feedD*2 + hordist, (i-1)*vertdist, 3, gndlayer ];
        anten.FeedLocations = [anten.FeedLocations;D2b+(Lb/sqrt(2))-feedD*2 + hordist, (i-1)*vertdist, 3, gndlayer ];
    end
    anten.ViaDiameter = feedD;

    anten.ViaLocations =[ -D2a-(La/sqrt(2))+Dsa , 0 , 1 , gndlayer ];
    anten.ViaLocations =[anten.ViaLocations; D2a+(La/sqrt(2))-Dsa , 0 , 1 , gndlayer ];
    for (i=2:norows)
        anten.ViaLocations =[anten.ViaLocations;  -D2a-(La/sqrt(2))+Dsa , (i-1)*vertdist , 1 , gndlayer ];
        anten.ViaLocations =[anten.ViaLocations; D2a+(La/sqrt(2))-Dsa , (i-1)*vertdist , 1 , gndlayer ];
    end
    anten.ViaLocations =[anten.ViaLocations; -D2b-(Lb/sqrt(2))+Dsb+ hordist , 0 , 3 , gndlayer ];
    anten.ViaLocations =[anten.ViaLocations; D2b+(Lb/sqrt(2))-Dsb+ hordist , 0 , 3 , gndlayer ];
    for (i=2:norows)
        anten.ViaLocations =[anten.ViaLocations;  -D2b-(Lb/sqrt(2))+Dsb + hordist , (i-1)*vertdist  , 3 , gndlayer ];
        anten.ViaLocations =[anten.ViaLocations; D2b+(Lb/sqrt(2))-Dsb + hordist, (i-1)*vertdist  , 3 , gndlayer ];
    end
   aa = [0 180];
    for (i=2:norows)
        aa =[aa 0 ];
        aa =[aa 180 ];
    end
    aa = [aa 0 180];
    for (i=2:norows)
        aa =[aa 0 ];
        aa =[aa 180 ];
    end
    anten.FeedPhase = aa;
    h = figure(2);
    show(anten);
%     txt1 = ['L: ' num2str(L*1000) ' mm'];
%     txt2 = ['D: ' num2str(D*1000) ' mm'];
%     txt3 = ['Ds: ' num2str(Ds*1000) ' mm'];
%     txt4 = ['gap: ' num2str(gap*1000) ' mm'];
%     % txt5 = ['50 ohm width: ' w_m*1000) ' mm'];
%     txt = {txt1,txt2,txt3,txt4};
%     annotation(h,'textbox', [0.05, 0.9, 0.1, 0.1], 'String',txt,'FitBoxToText','on'); 

end

