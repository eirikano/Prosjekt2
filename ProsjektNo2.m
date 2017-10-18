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
plot(x(1:15),C_an(1:15,:)*100)
grid
title('Concentration profile, analytic 1D')
xlabel('Position [µm]')
ylabel('% B')
leg = strtrim(cellstr(num2str((t./(60^2))'))');
legend(strcat(leg,'  hours'))

%iii)a)----------------------------

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
plot(x(1:15),C_num(1:15,j-1))
grid
legend('Analytic','Numeric')
title('Concentration profile, 1D')
xlabel('Position [µm]')
ylabel('Composition B')
leg = strtrim(cellstr(num2str((t./(60^2))'))');



%Isothermal annealing, analytical solution iii)b)--------------------


k_eq = @(C_i) 2*(C_i-v.C_0)/(v.C_p-v.C_0);

k=k_eq(C_i);

B_eq = @(k,t) c.B0 - (k/sqrt(pi))*sqrt(D_T*t);
B_norm(1)=1;
i=1;
t(1)=0;
dt=0.1;
while B_norm>0
B_norm(i)= B_eq(k,t(i))/c.B0;
t(i+1)=t(i)+dt;
i=i+1;
end

figure
subplot(2,2,1)
plot(t(1:length(t)-1),B_norm)
axis([0 20 0 1])
title('Plate dissolution, 1D')
ylabel('scaled volume fraction')
xlabel('time[s]')
grid



%iii)c)----------------------------
clear t
dt=0.001;
B_num(1)=c.B0;
B_num_norm(1)=1;

t(1)=0;
j=1;
while B_num_norm(j)>0
t(j+1)=t(j)+dt;
%Backwards euler
    B_num(j+1)=B_num(j)-dt*(k/2)*sqrt(D_T/(pi*t(j+1)));
    B_num_norm(j+1)=B_num(j+1)/c.B0;
j=j+1;
end


hold on
plot(t,B_num_norm)
legend('Analytic','Numeric')



<<<<<<< HEAD
=======
%iii)d)----------------------------
%Non isothermal case
v.T2=430+273; %[K]
v.T1=400+273; %[K]

D_T1=D_eq(v.T1);
D_T2=D_eq(v.T2);

Ci_T1=Ci_eq(v.T1);
Ci_T2=Ci_eq(v.T2);

k_T1=k_eq(Ci_T1);
k_T2=k_eq(Ci_T2);

clear t
dt=0.001;
B_num_noniso_1(1)=c.B0;
B_num_noniso_1_norm(1)=1;
B_num_noniso_2(1)=c.B0;
B_num_noniso_2_norm(1)=1;

t(1)=0;
j=1;
D_1=D_T1;
k_1=k_T1;
D_2=D_T1;
k_2=k_T1;

while B_num_noniso_2_norm(j)>0
t(j+1)=t(j)+dt;
if B_num_noniso_1_norm(j)<0.7
    D_1=D_T2;
    k_1=k_T2;
end
if B_num_noniso_2_norm(j)<0.3
    D_2=D_T2;
    k_2=k_T2;
end
%Backwards euler
    B_num_noniso_1(j+1)=B_num_noniso_1(j)-dt*(k_1/2)*sqrt(D_1/(pi*t(j+1)));
    B_num_noniso_1_norm(j+1)=B_num_noniso_1(j+1)/c.B0;
    
    B_num_noniso_2(j+1)=B_num_noniso_2(j)-dt*(k_2/2)*sqrt(D_2/(pi*t(j+1)));
    B_num_noniso_2_norm(j+1)=B_num_noniso_2(j+1)/c.B0;
j=j+1;
end

subplot(2,2,2)
plot(t,B_num_noniso_1_norm)
hold on
plot(t,B_num_noniso_2_norm)
axis([0 20 0 1])
title('Plate dissolution, non-isothermal, 1D')
ylabel('scaled volume fraction')
xlabel('time[s]')
grid




%tr1=(pi/D_T(T))*(c.B0/B0r).^2 
%t1star=trl*(k)
%scaledvolf=1-sqrt(t/t1star);
>>>>>>> 5e9f7e34c31fb049a03335680fa527e8af5524e2

%Isokinetic annealing, analytic solution iii)e)----------------------

k_eq = @(C_i) 2*(C_i-v.C_0)/(v.C_p-v.C_0);

k=k_eq(C_i);

B_eq = @(k,t) c.B0 - (k/sqrt(pi))*sqrt(D_T*t);
B_normE(1)=1;
i=1;
t(1)=0;
dt=0.1;
while B_normE>0.7
B_normE(i)= B_eq(k,t(i))/c.B0;
t(i+1)=t(i)+dt;
i=i+1;
end

