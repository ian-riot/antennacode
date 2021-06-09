% DXFtool example 2: Making a nice graphical legend 
% using dxf-pictures for two data series

clc; clear vars; close all

% some data
x1 = [25 38 50 75 100];
y1 = [136 108 95 82 80];
x2 = [9 20 40 80 160];
y2 = [141 125 104 89 74];

% open figure window & plot data
figure

plot(x1,y1,'ko','markerfacecolor','r'); hold on
plot(x2,y2,'ks','markerfacecolor','b');    

xlabel('Thickness, t [mm]')
ylabel('Fatigue strength, \Delta\sigma_{R} [MPa]')



% make small axes on top of the main one
axes('position',[0.67 0.67 .2 .2]);

% plot the dxf in the small axes
DXFtool('ex2.dxf');
axis off % hide ticks, labels, etc. in small axes


% show relevant marker next dxf picture
plot(40,63,'ko','markerfacecolor','r')
plot(40,33,'ks','markerfacecolor','b')

% add text to dxf picture
text(88,33,'t','fontsize',11)

