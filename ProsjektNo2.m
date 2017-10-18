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

x=logspace(log(c.B0),-5); 

xgrid=300;
xend=10^-3;
x=linspace(c.B0,xend,xgrid); 
dt=100; %time step
dx=x(2)-x(1); %distance step


%Analytic solution
C_an_eq = @(x,t) C_i-(C_i - v.C_0).*erf((x-c.B0)./(2*sqrt(D_T*t)));
for i=1:length(t)
    C_an(:,i)= C_an_eq(x,t(i));
end
figure
hold on
plot(x(1:15),C_an(1:15,:)*100)
grid
title('Concentration profile, analytic 2D')
xlabel('Position [µm]')
ylabel('% B')
leg = strtrim(cellstr(num2str((t./(60^2))'))');
legend(strcat(leg,'  hours'))
%Numerisk del Eirik
%Numerical iii)a)----------------------------

%Initial values


t=0;

%Boundary conditions (NEEDS TO BE CONFIRMED!)
C_num(1:length(x),1)=v.C_0;
C_num(1,1)=C_i;
%


%Stability
r=(D_T*dt)/(dx^2);
while r>0.5
    r=(D_T*dt)/(dx^2);
    dt=dt/2;
end

j=1;
while t<1*60*60 %hours
    C_num(1,j+1)=C_i; %boundary cond.
    for i=2:length(x)-1
        C_num(i,j+1)=C_num(i,j) + r*(C_num(i+1,j)-2*C_num(i,j)+C_num(i-1,j));
    end
    if i==length(x)-1
        C_num(i+1,j+1)=v.C_0; %boundary cond.
    end
j=j+1;
t=t+dt;
end


figure
plot(x(1:15),C_an(1:15,2))
hold on
plot(x(1:15),C_num(1:15,j-1),'.')
grid
legend('Analytic','Numeric')
title('Concentration profile, 2D')
xlabel('Position [µm]')
ylabel('Composition B')
leg = strtrim(cellstr(num2str((t./(60^2))'))');


%Isokinetic solution iii)b)----------------------------

%k_eq = @(C_i) 2*(C_i-v.C_0)/(v.C_p-v.C_0);

%B
%B = @(k) c.B0- (k/sqrt(pi))*sqrt(D_T*t);

%tr1=(pi/D_T(T))*(c.B0/B0r).^2 
%t1star=trl*()
%scaledvolf=1-sqrt(t/t1star);


