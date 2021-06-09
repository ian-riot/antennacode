function Zlist=bklop(Zo,Zload,N,RdB)
% Calculate impedance list for a Klopfenstein taper 
% of unit length. 
% 
% Zlist=bklop(Zo,Zload,N,RdB)
%
% Zo......Characteristic impedance (Ohms)
% Zload...Load impedance to match to (Ohms)
% N.......Number of sections used to approximate taper (integer)
% RdB.....Operating band ripple (dB)
%
% e.g.  Zlist=bklop(50,100,60,-25) % Match a 100ohm load to a 50ohm line
%                                  % with operating band ripple at -25dB.
%                                  % Taper defined as list of 60 sections  
%
% Note : Only valid for Zload>Zo
%
% Matches a load impedance Zload to a standard line impedance Zo
% using Klopfenstein taper. 
% Taper profile is returned as a list of impedances, the optimum length
% for this design of taper is 0.565 lambda.
%
%            Impedance Values
% Zo --->    [ Z1 ] [ Z2 ] ....    [ ZN ]   <-- Zload
%
% Ref D.M Pozar Microwave Engineering 2nd Ed Page 291

% N.Tucker www.activefrance.com 2010



   
Tld=log(Zload/Zo)*0.5;            % Reflection coefficient of load
Trip=10.^(RdB/20);                % Lin value of ripple in operating band
 

A=acosh(Tld/Trip);                % Intermediate variable in calculation


z=0;                              % Fractional distance along taper
dz=1/(N-1);                       % Incremental distance


Zx=zeros(1,N);
for c=1:N                         % Loop for impedance values along taper
 
 M=round(z*100+25);               % Number of steps for the numerical integration          
 PsiXA=0;
 y=0;
 dy=((2*z-1)/(M-2));              % Increment for numerical integration
 
 for d=1:M                        % Loop for PsiXA numerical integration
  PsiXA=PsiXA+besseli(1,A*sqrt(1-y.^2))/(A*sqrt(1-y.^2))*dy;
  y=y+dy;
 end

 % Calculate impedance as a function of distance along the 
 % unit length transformer
 LNZx=0.5*log(Zo*Zload)+(Tld/cosh(A))*(A.^2)*real(PsiXA);
 Zx(1,c)=exp(LNZx);
 z=z+dz;
end


Zlist=[Zo,Zx,Zload];          % Assemble the list of impedances for output
X=1:1:N;                      % X-axis vector for plotting

figure(10);
plot(X,(Zx),'b-',X,(Zx),'+');
xlabel('Zo     Matching Section Number    Zload');
ylabel('Impedance (Ohms)');
title('Transformer Impedances')
grid on;

chartname=sprintf(' Transformer Impedances ');
set(10,'name',chartname);