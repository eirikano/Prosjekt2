clear all
close all
clc

%-------------Constants------------------
c.D0 = 3.46*10^-5; % [m^2/s]
c.Q = 123.8; %[kJ/mol]
c.R=8.3145; %[J/K*mol]
c.Cstar=2.17*10^3; %[wt%]
c.dH_0=50.8; % [kJ/mol]
c.B0=0.001*10^-6; %[m]
%----------------------------------------

%-------------Variables------------------
v.T_iso=400+273; %[K]
v.C_p=100;
v.C_0=0;
%----------------------------------------

%Diffusion
D_eq = @(T) c.D*exp(-c.Q/(c.R*v.T));  % [mm^2/s]
D_T=D_eq(v.T_iso);

%Ci
Ci_eq = @(T) c.Cstar*exp(-c.dH_0/(c.R*v.T));
C_i=Ci_eq(v.T_iso);



%Analytic solution
C_an = @(x,t) C_i-(C_i - C_0)*erf((x-c.R)/(2*sqrt(D_T*t)));



%------------------------------------------
% %k
% k_eq = @(C_i) 2*(C_i-v.C_0)/(v.C_p-v.C_0);
% 
% %B
% B = @(k) c.B0- (k/sqrt(pi))*sqrt(D_T*t);
% 
% t1star=(pi/D_T(T))*(
% scaledvolf=1-sqrt(t/t1star);7

%eirik er kulere enn tobias
