close all
clear all


f1 = 5.6e9;


Z0 = 100;
er = 4.1;
h = 1.54e-3;


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
 w=Wdr*h; 

 w_m = w


% calculate wavelength on microstrip
A = (er+1)/2;
B = (er-1)/2;
C = sqrt(1+(12*h/w_m));

if ((w_m/h)<1)
   eeff = A+B*((1/C) + 0.04*((1-(w_m/h))^2)) ;  
else
    eeff = A+(B/C);
end

lam = 3e11/(f1*sqrt(eeff));

lam2 = lam/2



