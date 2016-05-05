close all;
clear all;
clc;
format long;

%1 LoadData%%%%%%%%%%%%%%%%%%%%%%%%
%mydata
nvd1=dlmread('ncJC200.txt');
%data fitting from device JC: pH7, Welec, 200 mM PB pH7 Hall Bar
nvd1n=nvd1(:,1);
nvd1Vg=nvd1(:,2);


%2DefiningConstants%%%%%%%%%%%%%%%%
hbar=6.58211814e-16; h=hbar*2*pi; vF=1000000;  %eVs, Fermi velocity in m/s
kT=8.617e-5*298; eps=8.854e-12; %boltzmann constant Ev, temperature, epsilon_0 
e=1.602e-19; G0=4*e/h; %unit charge, Qcond


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%      parameters     %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%alphaG
%LxL
offset=0;%offest of Vg to align to data
alpha= 1; %coupling coefficient

%Fluid based capacitance
Cdl=.2; %Farads/m^2
M=.1;%concentration in M
d=0.3*(M)^-.5*1/1000000000%in nm Debye screening length Sheehan; analytical chemistry 2012 
ep=5;%permittivity of fluid
Cdlv=eps*ep/d%variable double layer capacitance


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%     calculations    %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nEf=100; nE=500;
Efmin=-1; Efmax=1; 
%creating energy and fermi energy vectors
Efrange=Efmax-Efmin;Ef=Efmin:Efrange/nEf:Efmax;

%Graphene
%Density of states and dispersion relation 
kc=abs(Ef./(hbar*vF)); %lin disp relationship for graphene
DEc=2/pi*Ef./(hbar*vF)^2;%DEc==density of states
nn=(Ef.^2)/(pi*hbar^2*vF^2)/10000;%E in eV so carrier density 1/cm^2
Vg=(Ef/alpha)+offset;

%Determin Capacitance
%Cq=abs(e*2*Ef./(pi*hbar^2*vF^2));
%Ctemp=(Cq.*Cdl);
%Ctot=Ctemp./(Cq+Cdl);
%Vgq=e*10000*nn./Ctot;

nm=(Ef.^2)/(pi*hbar^2*vF^2);
Vge=hbar*vF*pi^.5*nm.^.5+e*nm./Cdl;

nmd=[(Ef.^2)/(pi*hbar^2*vF^2),(Ef.^2)/(pi*hbar^2*vF^2)];
Vgd=[-(hbar*vF*pi^.5*nm.^.5),(hbar*vF*pi^.5*nm.^.5)];
Vgdd=[-(hbar*vF*pi^.5*nm.^.5+e*nm./Cdl),(hbar*vF*pi^.5*nm.^.5+e*nm./Cdl)];
Vgdv=[-(hbar*vF*pi^.5*nm.^.5+e*nm./Cdlv),(hbar*vF*pi^.5*nm.^.5+e*nm./Cdlv)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%       plotting      %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
%plot(nvd1Vg,nvd1n*10000,Vg,nn*10000,Vge,nm);
%plot(nvd1Vg,nvd1n*10000,Vge,nm,Vgd,nmd);
plot(nvd1Vg-.03,nvd1n,Vgd,nmd/10000,'.',Vgdd,nmd/10000,'.',Vgdv,nmd/10000,'.');
xlabel('$V_g$ $(V)$')
ylabel('$n$ x $10^{13}$ ($\frac{1}{cm^2}$)')
ylim([0 60e12])
xlim([-1 1])
set(0,'defaulttextinterpreter','latex')
set(0,'defaultaxesfontname','helvetica')
set(gca,'fontsize',16)
%set(gca,'YTickLabelMode','manual')


