clear all
close all
clc

%-------------Constants------------------
c.D0 = 3.46*10^-5; % [m^2/s]
c.Q = 123800; %[J/mol]
c.R=8.3145; %[J/K*mol]
c.Cstar=2.17*10^3; %[wt%]
c.dH_0=50800; % [J/mol]
c.B0=0.001*10^-6; %[m]
%----------------------------------------

%-------------Variables------------------
v.T_iso=400+273; %[K]
v.C_p=100;
v.C_0=0;
%----------------------------------------

%Diffusion
D_eq = @(T) c.D0*exp(-c.Q/(c.R*T));  % [mm^2/s]
D_T=D_eq(v.T_iso);

%Ci
Ci_eq = @(T) c.Cstar*exp(-c.dH_0/(c.R*T));
C_i=Ci_eq(v.T_iso);

t=[0.5,1,2,3]*60*60; %hours
x=linspace(10^-6,10^-5,100); 

%Analytic solution
C_an_eq = @(x,t) C_i-(C_i - v.C_0).*erf((x-c.B0)./(2*sqrt(D_T*t)));
for i=1:length(t)
    C_an(:,i)= C_an_eq(x,t(i));
end
figure
hold on
plot(x,C_an)
grid
title('Concentration profile')
xlabel('Position [µm]')
ylabel('Composition B')
leg = strtrim(cellstr(num2str((t./(60^2))'))');
legend(strcat(leg,'  hours'))
%Numerisk del Eirik
%Numerical iii)a)----------------------------

%Initial values

dt=0.01; %time step
dx=10^-6; %distance step

%Boundary conditions (NEEDS TO BE CONFIRMED!)
C_num(1:length(x),1)=v.C_0;
C_num(c.B0,1:length(t))=v.C_i;


%Stability
r=(D_T*dt)/(dx^2);
if r>0.5
    display('stability error')
end

for i=1:length(x)
    for j=1:length(t)
        C_num(i+1,j)=C(i,j) + r*(C(i+1,j)-2*C(i,j)+C(i-1,j));
    end
end


%Isokinetic solution iii)b)----------------------------

%k
%k_eq = @(C_i) 2*(C_i-v.C_0)/(v.C_p-v.C_0);

%B
%B = @(k) c.B0- (k/sqrt(pi))*sqrt(D_T*t);

%t1star=(pi/D_T(T))*(
%scaledvolf=1-sqrt(t/t1star);7

%lukta sjukt myggsprøy
