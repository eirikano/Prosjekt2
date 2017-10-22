clear all
clear clc
close all

%-------------Constants------------------
c.D0 = 3.46*10^-5; % [m^2/s]
c.Q = 123800; %[J/mol]
c.R=8.3145; %[J/K*mol]
c.Cstar=2.17*10^3; %[wt%]
c.dH_0=50800; % [J/mol]
c.B0=0.025*10^-6; %[m]
%----------------------------------------

%-------------Variables------------------
v.T_1=400+273; %[K]
v.T_2=430+273;
v.C_p=100;
v.C_0=0;
%----------------------------------------

%-------------Diffusion------------------
D_eq = @(T) c.D0*exp(-c.Q/(c.R*T));
D_T=D_eq(v.T_1);
%----------------------------------------
%-----------------Ci---------------------
Ci_eq = @(T) c.Cstar*exp(-c.dH_0/(c.R*T));
C_i=Ci_eq(v.T_1);
%----------------------------------------

%------------------k---------------------
k_eq = @(C_i) 2*(C_i-v.C_0)/(v.C_p-v.C_0);
k_T=k_eq(C_i);
%----------------------------------------

t=[10^-5,10^-4,10^-3,3]; %seconds

xgrid=300;
xend=7*10^-7;
x=linspace(c.B0,xend,xgrid); 
dx=x(2)-x(1); %distance step


%%
%iii)a)----------------------------
%Analytic solution of diffusion profile for spherical particle
C_an_eq = @(x,t) v.C_0+(C_i - v.C_0).*(c.B0./x).*erfc((x-c.B0)./(2*sqrt(D_T*t)));
figure
subplot(2,1,1)
for i=1:length(t)
    C_an(:,i)= C_an_eq(x,t(i));
    plot(x,C_an(:,i))
    hold on
    anlegend{i}=[num2str(t(i)) 'seconds'];
end
grid
title('a)Concentration profile, analytic 3D')
xlabel('Position [m]')
ylabel('C')
axis([c.B0 10^-7 0 C_i])
legend(anlegend)
%Numeric solution of diffusion profile for spherical particle
dt=10^-2;
Cstart(1:length(x),1)=v.C_0;
Cstart(1,1)=C_i;
while dt>(dx^2)/(2*D_T) %adjusting dt to stability requirement
    dt=dt/2;
end

t=[0:dt:3];
for j =1:length(t) 
    Cnext(:,j)=Cnumspher(v,dt,D_T,x,Cstart);
    Cstart=Cnext(:,j);
end
C_numeric(:,1)=Cnext(:,j-1);

subplot(2,1,2)
plot(x,C_an(:,4))
hold on
plot(x,C_numeric(:,1),'*')
grid
legend(['Analytic ' num2str(t(length(t))) ' seconds'],['Numeric ' num2str(t(length(t))) ' seconds'])
axis([c.B0 10^-5 0 C_i])
title('a)Concentration profile, 1D')
xlabel('Position [m]')
ylabel('C')
leg = strtrim(cellstr(num2str((t./(60^2))'))');

%%
%iii)b)----------------------------
%By numerical solving of equation 16
clear t
t=[0.001:dt:13.5]; %can not divide by zero!
r(1)=c.B0;
for j =1:length(t)
    r(j+1)=r(j)-dt*(((k_T*D_T)/(2*r(j)))+(k_T/2)*sqrt(D_T/(pi*t(j))));
end
figure
plot(t,r(1:length(r)-1)./c.B0)
hold on
grid
title('b and c : Numeric vs short time solution eq. 16')
xlabel('Time [s]')
ylabel('normalized radius')

    
%%
%iii)c)----------------------------
%Approximate solution for short times
clear t
clear r
t(1)=0;
r(1)=c.B0;
j=1;
while t<13.5
    r(j+1)=c.B0-(k_T*D_T*t(j)/(2*c.B0))-(k_T/sqrt(pi))*sqrt(D_T*t(j));
    t(j+1)=t(j)+dt;
    j=j+1;
end

plot(t,r./c.B0)
hold on


%%
%iii)d)----------------------------
%Spherical precipitate numeric
%isothermal case
clear t
clear Cnext
clear Cstart
t(1)=0;
r_num_iso_norm(1)=1;
r_num_iso(1)=c.B0;
Cstart(1:length(x),1)=v.C_0;
Cstart(1,1)=C_i;
j=1;
dr=0;

indexlog=0;
timeindex=0;
while r_num_iso_norm(j)>0
    t(j+1)=t(j)+dt;
    x=linspace(r_num_iso(j),x(length(x))+dr,xgrid); %move the x-vector with the decreasing particle size
    r_num_iso(j+1)=r_num_iso(j)+((dt*D_T)/(dx*(v.C_p-C_i)))*(Cstart(2,1)-Cstart(1,1));
    r_num_iso_norm(j+1)=((r_num_iso(j+1))/c.B0)^3;
    if (r_num_iso_norm(j+1)<0.7)&&(indexlog==0) %for use later
        n=j+1;
        indexlog=1;
    end
    
    if (t(j)>13.5)&&(timeindex==0) %just for comparison with eq 16
        tlog=j;
        timeindex=1;
    end
    
    Cnext(:,j)=Cnumspher(v,dt,D_T,x,Cstart(:,1));
    Cstart(:,1)=Cnext(:,j);
    
    dr=r_num_iso(j+1)-r_num_iso(j);
    j=j+1;
end
%Note! One can optimize the code by using the above concentration profile
%and r_num_iso_norm (fraction) in the non-isothermal case (only valid
%before the temperature increase)
plot(t(1:tlog),r_num_iso(1:tlog)./c.B0)
legend('Numeric eq.16','Short time apprx eq.16','Numeric solution');

figure
subplot(3,1,1)
plot(t,r_num_iso_norm)
axis([0 inf 0 1])
title('d)Spherical dissolution, isotherm, 3D')
ylabel('scaled volume fraction')
xlabel('time[s]')
grid

%Spherical precipitate numeric
%Non isothermal case Two-Step (low to high temp)
x=linspace(c.B0,xend,xgrid);
dt=10^-2;

j=n; %starting from saved index value from isothermal task
D_1=D_eq(v.T_2);
Ci_1=Ci_eq(v.T_2);
D_2=D_eq(v.T_1);
Ci_2=Ci_eq(v.T_1);

while dt>(dx^2)/(2*D_1) %adjusting dt to stability requirement
    dt=dt/2;
end

Cstart1=Cnext(:,n-1);
Cstart2=Cnext(:,n-1);
r_num_noniso_1_norm=r_num_iso_norm(:,1:n);
r_num_noniso_1=r_num_iso(:,1:n);
r_num_noniso_2_norm=r_num_iso_norm(:,1:n);
r_num_noniso_2=r_num_iso(:,1:n);
t=t(1:n);

alreadychanged2=0;
x1=x;
x2=x;
dr1=0;
dr2=0;
while (r_num_noniso_1_norm(j)>0)||(r_num_noniso_2_norm(j)>0)
    t(j+1)=t(j)+dt;    
    x1=linspace(r_num_noniso_1(j),x1(length(x1))+dr1,xgrid); %move the x-vector with the decreasing particle size
    x2=linspace(r_num_noniso_2(j),x1(length(x2))+dr2,xgrid); 
    if (r_num_noniso_2_norm(j)<0.3)&&(alreadychanged2==0)
        D_2=D_eq(v.T_2);
        Ci_2=Ci_eq(v.T_2);
        alreadychanged2=1;
    end
    Cstart1(1,1)=Ci_1;
    Cnext1(:,j)=Cnumspher(v,dt,D_1,x1,Cstart1(:,1));
    Cstart1(:,1)=Cnext1(:,j);
    
    Cstart2(1,1)=Ci_2;
    Cnext2(:,j)=Cnumspher(v,dt,D_2,x,Cstart2(:,1));
    Cstart2(:,1)=Cnext2(:,j);

    r_num_noniso_1(j+1)=r_num_noniso_1(j)+((dt*D_1)/(dx*(v.C_p-Ci_1)))*(Cnext1(2,j)-Cnext1(1,j));
    r_num_noniso_1_norm(j+1)=((r_num_noniso_1(j+1))/c.B0)^3;
    
    r_num_noniso_2(j+1)=r_num_noniso_2(j)+((dt*D_2)/(dx*(v.C_p-Ci_2)))*(Cnext2(2,j)-Cnext2(1,j));
    r_num_noniso_2_norm(j+1)=((r_num_noniso_2(j+1))/c.B0)^3;
    
    dr1=r_num_noniso_1(j+1)-r_num_noniso_1(j);
    dr2=r_num_noniso_2(j+1)-r_num_noniso_2(j);
    j=j+1;
end

subplot(3,1,2)
plot(t,r_num_noniso_1_norm)
hold on
plot(t,r_num_noniso_2_norm)
axis([0 inf 0 1])
title('d)Spherical dissolution, two-step, 3D')
ylabel('scaled volume fraction')
xlabel('time[s]')
grid

%Spherical precipitate numeric
%Non isothermal case Two-Step
%Only change is switching of T1 and T2

clear t
clear r_num_noniso_1
clear r_num_noniso_1_norm
clear r_num_noniso_2
clear r_num_noniso_2_norm
clear Cnext1
clear Cnext2
dt=10^-2;
t(1)=0;
j=1;
D_1=D_eq(v.T_2);
Ci_1=Ci_eq(v.T_2);
D_2=D_eq(v.T_2);
Ci_2=Ci_eq(v.T_2);
while dt>(dx^2)/(2*D_1) %adjusting dt to stability requirement
    dt=dt/2;
end

r_num_noniso_1(1)=c.B0;
r_num_noniso_1_norm(1)=1;
r_num_noniso_2(1)=c.B0;
r_num_noniso_2_norm(1)=1;

Cstart1(1:length(x),1)=v.C_0;
Cstart1(1,1)=Ci_1;
Cstart2=Cstart1;
x1=x;
x2=x;
dr1=0;
dr2=0;
alreadychanged1=0;
alreadychanged2=0;
while (r_num_noniso_1_norm(j)>0)||(r_num_noniso_2_norm(j)>0)
    t(j+1)=t(j)+dt;    
    x1=linspace(r_num_noniso_1(j),x1(length(x1))+dr1,xgrid); %move the x-vector with the decreasing particle size
    x2=linspace(r_num_noniso_2(j),x1(length(x2))+dr2,xgrid); 
    if (r_num_noniso_1_norm(j)<0.7)&&(alreadychanged1==0)
        D_1=D_eq(v.T_1);
        Ci_1=Ci_eq(v.T_1);
        alreadychanged1=1;
        while dt>(dx^2)/(2*D_1) %adjusting dt to stability requirement
            dt=dt/2;
        end
    end
    if (r_num_noniso_2_norm(j)<0.3)&&(alreadychanged2==0)
        D_2=D_eq(v.T_1);
        Ci_2=Ci_eq(v.T_1);
        alreadychanged2=1;
    end
    Cstart1(1,1)=Ci_1;
    Cnext1(:,j)=Cnumspher(v,dt,D_1,x,Cstart1(:,1));
    Cstart1(:,1)=Cnext1(:,j);
    
    Cstart2(1,1)=Ci_2;
    Cnext2(:,j)=Cnumspher(v,dt,D_2,x,Cstart2(:,1));
    Cstart2(:,1)=Cnext2(:,j);

    r_num_noniso_1(j+1)=r_num_noniso_1(j)+((dt*D_1)/(dx*(v.C_p-Ci_1)))*(Cnext1(2,j)-Cnext1(1,j));
    r_num_noniso_1_norm(j+1)=((r_num_noniso_1(j+1))/c.B0)^3;
    
    r_num_noniso_2(j+1)=r_num_noniso_2(j)+((dt*D_2)/(dx*(v.C_p-Ci_2)))*(Cnext2(2,j)-Cnext2(1,j));
    r_num_noniso_2_norm(j+1)=((r_num_noniso_2(j+1))/c.B0)^3;
    
    dr1=r_num_noniso_1(j+1)-r_num_noniso_1(j);
    dr2=r_num_noniso_2(j+1)-r_num_noniso_2(j);
    
    j=j+1;
end

subplot(3,1,3)
plot(t,r_num_noniso_1_norm)
hold on
plot(t,r_num_noniso_2_norm)
axis([0 inf 0 1])
title('d)Spherical dissolution, two-step, 3D')
ylabel('scaled volume fraction')
xlabel('time[s]')
grid


