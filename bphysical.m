function bphysical(Zlist,L,Fo,Er,d,filename);
% Generates the physical realisation, in microstrip 
% of the impedance list as calculated by bmatch.m
%
% The results are output as polygon .DXF file that can
% be imported directly into an EM simulator.
%
% Usage : bphysical(Zlist,L,Fo,Er,d,filename)
%
% Zlist.....Impedance list returned by bmatch (Ohms)
% L.........Length of transformer sections in wavelengths
% Fo........Centre frequency in (MHz)
% Er........Dielectric constant
% d.........Dielectric thickness (mm)
% filename..Full pathname for the .DXF file (string)
%                                               
%
% E.g. To calculate the dimensions for a transformer with Fo=1000MHz 
% using 1/100 wave transformer sections on 0.76mm board Er=3.48 use 
% the following :-
%
% filename='c:\matlab\toolbox\RFutils_M\dxf\match.dxf';
% bphysical(Zlist,0.01,1000,3.48,0.76,filename) 


% Reference Microwave Engineering 2nd Ed Page162   D.M. Pozar
% N. Tucker ActiveFrance.com 2010


vo=3e8;
lambda=vo/(Fo*1e6);
ko=2*pi/lambda;  
Lo=L*lambda*1e3;     % Free space length of transformer section (mm)

[Row,Col]=size(Zlist);

Ztran=Zlist;                         
N=Col;                      

Wx=zeros(1,N);
for x=1:N
  
 Zx=Ztran(1,x);

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
 Wx(1,x)=W;

 % Calculate effictive dielectric constant for microstrip
 % line of width W on dielectric material of constant Er
 Ereff=((Er+1)/2)+((Er-1)/2)*(1+12*(d/(W/1e3))).^-0.5;
 Lx(1,x)=Lo/sqrt(Ereff);
  
end


Lx=Lx;         % Convert length to (mm)
X=cumsum(Lx);  % Cumulative length of transformer (x-axis for plotting)
Ltot=X(1,N);   % Total length of transformer (mm)
Wx=Wx;         % Convert line width to (mm)


% Add extra points at each end to allow plotting of Zo / Load sections

Xp1=[X,(X(1,N)+Lx(1,N))]-Lx(1,1); % Concatenate array of X-coords (dist along line)
Wxp1=[Wx,Wx(1,N)];                % Concatenate array of Y-coords (width profile)


% Construct the plotting point pairs if there are less that 10 transformers. 
% Remembering that N includes the values of Zo and Zload, hence N<12
% This is so the steps in line width are plotted not just lines between the points,
% However for transformers with many sections it is advantageous to interpolate
% to reduce the edge complexities for the EM solver or board fabricator.

if N<12        % Interpolate option 

  x1=1;
  Xp(1,1)=Xp1(1,1);
  for x=1:(N+0)
     Xp(1,x1+1)=Xp1(1,x+1);
     Xp(1,x1+2)=Xp1(1,x+1);
   
     Wxp(1,x1)=Wxp1(1,x);
     Wxp(1,x1+1)=Wxp1(1,x);
     x1=x1+2;
  end
  Wxp(1,x1)=Wxp1(1,N+1);
else
   Xp=Xp1;    % Plot individual steps option
   Wxp=Wxp1;
end   


figure(11);
plot(Xp,Wxp/2,Xp,-Wxp/2);
xlabel('Zo                   Distance along transformer (mm)                   Zload');
ylabel('Line profile (mm)')
T1=sprintf('Microstrip Dimensions Er=%g  d=%g',Er,d);
title(T1);

chartname=sprintf(' Microstrip Physical Profile ');
set(11,'name',chartname);

% **************** Write DXF file **************

% Concatenate top and bottom line profile data to produce
% data to plot the polygon.

Xpoly=[0,Xp,fliplr(Xp),0];
Ypoly=[0,Wxp/2,-fliplr(Wxp/2),0];
Zpoly=zeros(size(Xpoly));

pcolour=5; % 1=red,2=yellow,3=green,4=cyan,5=blue,6=magenta
layername='MATCH';

mat2dxfp(Xpoly,Ypoly,Zpoly,5,layername,filename);

