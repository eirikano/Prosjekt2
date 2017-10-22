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

%-------------Diffusion------------------
D_eq = @(T) c.D0*exp(-c.Q/(c.R*T));  % [mm^2/s]
D_T=D_eq(v.T_iso);
%----------------------------------------
%-----------------Ci---------------------
Ci_eq = @(T) c.Cstar*exp(-c.dH_0/(c.R*T));
C_i=Ci_eq(v.T_iso);
%----------------------------------------

t=[0.5,1,2,3]*60*60; %hours

xgrid=5000;
xend=5*10^-3;
x=linspace(c.B0,xend,xgrid); 
dt=1; %time step
dx=x(2)-x(1); %distance step


%Analytic solution
C_an_eq = @(x,t) C_i-(C_i - v.C_0).*erf((x-c.B0)./(2*sqrt(D_T*t)));
figure
for i=1:length(t)
    C_an(:,i)= C_an_eq(x,t(i));
    plot(x,C_an(:,i))
    anlegend{i}=[num2str(t(i)/3600) 'hours'];
    hold on
end
grid
title('a)Concentration profile, analytic 1D')
xlabel('Position [�m]')
ylabel('C')
axis([c.B0 3*10^-5 0 C_i])
legend(anlegend)


%iii)a)-------------------------------------------------------------

%Initial conc.
C_i=Ci_eq(v.T_iso);

%Diffusion
D_T=D_eq(v.T_iso);

%Stability
r=(D_T*dt)/(dx^2);
while r>0.5
    dt=dt/4;
    r=(D_T*dt)/(dx^2);
end
t=0;
dt=1;
Cstart(1:length(x),1)=v.C_0;
Cstart(1,1)=C_i;
j=1;
while t<60*60
    Cnext(:,j)=Cnum(v,r,x,Cstart);
    Cstart=Cnext(:,j);
    j=j+1;
    t=t+dt;
end
C_numeric(:,1)=Cnext(:,j-1);
figure
plot(x,C_an(:,2))
hold on
plot(x,C_numeric(:,1))
grid
legend('Analytic','Numeric')
axis([c.B0 3*10^-5 0 C_i])
title('a)Concentration profile, 1D')
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
title('b)Plate dissolution, 1D')
ylabel('scaled volume fraction')
xlabel('time[s]')
grid



%iii)c)Plate and spherical isothermal
clear t
dt=0.01;
B_num(1)=c.B0;
B_num_norm(1)=1;

r_num(1)=c.B0;
r_num_norm(1)=1;

t(1)=0;
j=1;
while (B_num_norm(j)>0)
t(j+1)=t(j)+dt;
%Backwards euler
    %Plate-like precipitates
    B_num(j+1)=B_num(j)-dt*(k/2)*sqrt(D_T/(pi*t(j+1)));
    B_num_norm(j+1)=B_num(j+1)/c.B0;
    
j=j+1;
end


hold on
plot(t,B_num_norm)
legend('Analytic','Numeric')




%iii)d)----------------------------
%Non isothermal case Two-Step
T2=430+273; %[K]
T1=400+273; %[K]

clear t
xgrid=7000;
xend=5*10^-5;
x=linspace(c.B0,xend,xgrid); 
dt=0.001; %time step
dx=x(2)-x(1); %distance step

t(1)=0;
j=1;
D_1=D_eq(T1);
Ci_1=Ci_eq(T1);
D_2=D_eq(T1);
Ci_2=Ci_eq(T1);

B_num_noniso_2(1)=c.B0;
B_num_noniso_2_norm(1)=1;
B_num_noniso_1(1)=c.B0;
B_num_noniso_1_norm(1)=1;

r= @(D) (D*dt)/(dx^2);
r1=r(D_1);
r2=r(D_1);
Cstart1(1:length(x),1)=v.C_0;
Cstart1(1,1)=C_i;
Cstart2=Cstart1;
alreadychanged1=0;
alreadychanged2=0;
while B_num_noniso_2_norm(j)>0
    t(j+1)=t(j)+dt;
if (B_num_noniso_1_norm(j)<0.7)&&(alreadychanged1==0)
    D_1=D_eq(T2);
    Ci_1=Ci_eq(T2);
    r1=r(D_2);
    if r1>0.5
        display('Numerical stability error: decrease dt or increase dx')
    end
    alreadychanged1=1;
end
if (B_num_noniso_2_norm(j)<0.3)&&(alreadychanged2==0)
    D_2=D_eq(T2);
    Ci_2=Ci_eq(T2);
    r2=r(D_2);
    if r2>0.5
        display('Numerical stability error: decrease dt or increase dx')
    end
    alreadychanged2=1;
end
Cstart1(1,1)=Ci_1;
Cstart2(1,1)=Ci_2;
Cnext1(:,j)=Cnum(v,r1,x,Cstart1(:,1));
Cnext2(:,j)=Cnum(v,r2,x,Cstart2(:,1));

Cstart1(:,1)=Cnext1(:,j);
Cstart2(:,1)=Cnext2(:,j);

B_num_noniso_1(j+1)=B_num_noniso_1(j)+((dt*D_1)/(dx*(v.C_p-Ci_1)))*(Cnext1(2,j)-Cnext1(1,j));
B_num_noniso_1_norm(j+1)=(B_num_noniso_1(j+1))/c.B0;

B_num_noniso_2(j+1)=B_num_noniso_2(j)+((dt*D_2)/(dx*(v.C_p-Ci_2)))*(Cnext2(2,j)-Cnext2(1,j));
B_num_noniso_2_norm(j+1)=(B_num_noniso_2(j+1))/c.B0;

j=j+1;
end

subplot(2,2,2)
plot(t,B_num_noniso_1_norm)
hold on
plot(t,B_num_noniso_2_norm)
axis([0 inf 0 1])
title('d)Plate dissolution, two-step, 1D')
ylabel('scaled volume fraction')
xlabel('time[s]')
grid


%iii)e)----------------------------
%Non isothermal case Two-Step
T1=430+273; %[K]
T2=400+273; %[K]

clear t
clear B_num_noniso_2
clear B_num_noniso_1
clear B_num_noniso_2_norm
clear B_num_noniso_1_norm

t(1)=0;
j=1;
D_1=D_eq(T1);
Ci_1=Ci_eq(T1);
D_2=D_eq(T1);
Ci_2=Ci_eq(T1);

B_num_noniso_2(1)=c.B0;
B_num_noniso_2_norm(1)=1;
B_num_noniso_1(1)=c.B0;
B_num_noniso_1_norm(1)=1;

r1=r(D_1);
r2=r(D_1);
Cstart1(1:length(x),1)=v.C_0;
Cstart1(1,1)=C_i;
Cstart2=Cstart1;
alreadychanged1=0;
alreadychanged2=0;
while B_num_noniso_2_norm(j)>0
    t(j+1)=t(j)+dt;
if (B_num_noniso_1_norm(j)<0.7)&&(alreadychanged1==0)
    D_1=D_eq(T2);
    Ci_1=Ci_eq(T2); 
    r1=r(D_2);
    if r1>0.5
        display('Numerical stability error: decrease dt or increase dx')
    end
end
if (B_num_noniso_2_norm(j)<0.3)&&(alreadychanged2==0)
    D_2=D_eq(T2);
    Ci_2=Ci_eq(T2);
    r2=r(D_2);
    if r2>0.5
        display('Numerical stability error: decrease dt or increase dx')
    end
end
Cstart1(1,1)=Ci_1;
Cstart2(1,1)=Ci_2;
Cnext1(:,j)=Cnum(v,r1,x,Cstart1(:,1));
Cnext2(:,j)=Cnum(v,r2,x,Cstart2(:,1));

Cstart1(:,1)=Cnext1(:,j);
Cstart2(:,1)=Cnext2(:,j);

B_num_noniso_1(j+1)=B_num_noniso_1(j)+((dt*D_1)/(dx*(v.C_p-Ci_1)))*(Cnext1(2,j)-Cnext1(1,j));
B_num_noniso_1_norm(j+1)=(B_num_noniso_1(j+1))/c.B0;

B_num_noniso_2(j+1)=B_num_noniso_2(j)+((dt*D_2)/(dx*(v.C_p-Ci_2)))*(Cnext2(2,j)-Cnext2(1,j));
B_num_noniso_2_norm(j+1)=(B_num_noniso_2(j+1))/c.B0;

j=j+1;
end

figure
plot(t,B_num_noniso_1_norm)
hold on
plot(t,B_num_noniso_2_norm)
axis([0 inf 0 1])
title('Plate dissolution, two-step, 1D')
ylabel('scaled volume fraction')
xlabel('time[s]')
grid

%Isokinetic annealing, analytic solution iii)e)----------------------
%% sgsge
D_1=D_eq(T2);
k=k_eq(C_i);

t1star=t1star_eq(D_T,k);
i=1;
dt=0.0005;
t=[0:dt:20];

f_low11=1;
f_low12=1;
bytte07=0;
bytte03=0;
sum=0;

while f_low11(i) > 0 
    
    sum(i+1)=sum(i)+(dt/t1star);
    f_low11(i+1)=1-sqrt(sum(i));
    if (f_low11(i) <= 0.7) && (bytte07==0)
        C_i=Ci_eq(T1);
        D_1=D_eq(T1);
        k=k_eq(C_i);
        t1star=t1star_eq(D_1,k);
        tid_low07 = i;
        bytte07=1;
    end
    i=i+1;
end
tid1=i;


D_T=D_eq(T2);
C_i=Ci_eq(T2);
k=k_eq(C_i);

t1star=t1star_eq(D_T,k);
i=1;

sum=0;
while f_low12 > 0

    sum(i+1)=sum(i)+(dt/t1star);
    f_low12(i+1)=1-sqrt(sum(i));
    if (f_low12(i) <= 0.3) && (bytte03==0)
        C_i=Ci_eq(T1);
        D_1=D_eq(T1);
        k_e=k_eq(C_i);
        t1star=t1star_eq(D_1,k_e);
        tid_low03 = i;
        bytte03=1;
    end
    i=i+1;
    
end
tid2=i;


figure 
plot (t(1:tid1), f_low11, t(1:tid2), f_low12)
hold on
plot
axis([0 inf 0 1])
title('isokinetic solution')
ylabel('scaled volume fraction')
xlabel('time[s]')
