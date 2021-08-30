% simple plotter for array
close all;

% General board
L = 200; %mm
W = 200; %mm
epsr = 4.6;
diss = 0.086;
h = 1.6; %mm

% 2.4 GHz patch
L2 = 26.13; %mm
W2 = 33.15; %mm
y02 = 2.74; %mm
w2 = 2; % stripline multiplier
f0_hz = 2.45; %GHz


% 5 GHz patch
L5 = 11.7; %mm
W5 = 13.2; %mm
y05 = 0.335; %mm
w5 = 2; % stripline multiplier
f0_hz = 5.6; %GHz

% stripline widths
w_100 = 0.531; %mm
w_50 = 2.3891; %mm

% 2.4GHz array
rows2 = 4;
cols2 = 2;
rowspacing2 = 44.4; %mm
colspacing2 = 59.3; %mm
% rowspacing2 = L2*1.7; %mm
% colspacing2 = L2*2.27; %mm

% 5GHz array
rows5 = 8;
cols5 = 2;
rowspacing5 = 20; %mm
colspacing5 = 26.7; %mm
% rowspacing5 = L5*1.7; %mm
% colspacing5 = L5*2.27; %mm

% Klopfenstein taper
L_klop = 20; %mm

% Miter
% calculated

% radius of 100 ohm lines
r_100 = 5; % x the width

% stripline clearance from patch
clear2 = L2/5;
clear5 = L5/5;

% array centers
X2 = W/2;
Y2 = L/4 + clear2*2;
X5 = W/2;
Y5 = L*3/4 + clear5*5;

% create board
xboard = [ 0 , 0 , W , W , 0];
yboard = [ 0 , L , L , 0 , 0 ];

% create primitive 2.4 GHz patch
slotwidth = w2*w_100;
ypatch2 = [-W2/2 , -W2/2 , W2/2 , W2/2 , slotwidth/2 , slotwidth/2 , w_100/2, w_100/2 ,- w_100/2 ,-w_100/2 , -slotwidth/2 , -slotwidth/2 , -W2/2];
xpatch2 = [-L2/2 , L2/2 , L2/2 , -L2/2 , -L2/2 , -L2/2 + y02 ,  -L2/2 + y02 , -L2/2 , -L2/2 , -L2/2 + y02 , -L2/2 + y02 , -L2/2 , -L2/2];

% create primitive 5 GHz patch
ypatch5 = [-W5/2 , -W5/2 , W5/2 , W5/2 , slotwidth/2 , slotwidth/2 , w_100/2, w_100/2 ,- w_100/2 ,-w_100/2 , -slotwidth/2 , -slotwidth/2 , -W5/2];
xpatch5 = [-L5/2 , L5/2 , L5/2 , -L5/2 , -L5/2 , -L5/2 + y05 ,  -L5/2 + y05 , -L5/2 , -L5/2 , -L5/2 + y05 , -L5/2 + y05 , -L5/2 , -L5/2];

% make 2.4G array
xarray2 = [];
yarray2 = [];
startpoint2x = -(cols2/2)*colspacing2/2;
startpoint2y = -(rows2/2)*rowspacing2/2 - rowspacing2/2;
for x=1:cols2
    xoffset = startpoint2x + (x-1)*colspacing2 + Y2;
    for y=1:rows2
        yoffset = startpoint2y + (y-1)*rowspacing2 +X2;
        xarray2 = [ xarray2 NaN xpatch2+xoffset ];
        yarray2 = [ yarray2 NaN ypatch2+yoffset ];
    end
end

% make 5G array
xarray5 = [];
yarray5 = [];
startpoint5x = -(cols5/2)*colspacing5/2;
startpoint5y = -(rows5/2)*rowspacing5/2 - rowspacing5*3/2;
for x=1:cols5
    xoffset = startpoint5x + (x-1)*colspacing5 + Y5;
    for y=1:rows5
        yoffset = startpoint5y + (y-1)*rowspacing5 +X5;
        xarray5 = [ xarray5 NaN xpatch5+xoffset ];
        yarray5 = [ yarray5 NaN ypatch5+yoffset ];
    end
end

% make total array
xarray = [xarray2 NaN xarray5];
yarray = [yarray2 NaN yarray5];

plot(xarray,yarray);
hold on;
plot(xboard,yboard);
plot(Y2,X2,'rx');
plot(Y5,X5,'bx');
