clc
clear all
global C20 A20 la om uu20
load('C20')
load('A20')
load('uu20')
uu=uu20;
u=uu(end);
A=A20;

%% User Entry
tS=400*86400;                      % Total time of simulation (s)
alpha=0.0008;
c1=3/10^6*14.7/101325;              %1/Pa
c2=3/10^6*14.7/101325;              %1/Pa
co=5/10^6*14.7/101325; 
vis=0.002;                          %kg/(m.s)                         
bo=1;
Pi=4000*101325/14.7;                %Pa
rw=1*0.3048;                         %m
re=5000*0.3048;                      %m
Q=100*(0.3048^3)/86400;              %m^3/s
h=30*0.3048;                        %m
k1=1*10^(-15);                      %m^2
phi1=0.2;
k2=21*10^(-15);                   %m^2
phi2=0.1;
om=phi2*c2/(phi1*c1+phi2*c2);
la=alpha*k1*re^2/k2; 
%% ODEs Solution
tDS=(k2*tS)/((phi1*c1+phi2*c2)*vis*re^2);
pD0=zeros(1,2*(length(uu)-1));
[T,Y]= ode45(@ODE20,[0 tDS],pD0);

%% Dimensionless Graphs
figure(1)
for i=1:length(uu)-1
    plot(T,Y(:,i),'linewidth',2);
    hold on
end
hold off
    xlim([0 tDS])
xlabel('Dimensionless time')
ylabel('Dimensionless pressure')
% legend('1st collocation point','2ed collocation point','3th collocation point','4th collocation point','5th collocation point'...
%    ,'6th collocation point','7th collocation point','8th collocation point','9th collocation point','10th collocation point','Location','SouthEast' )

%% Wellbore dimensionless pressure (r=rw)

Pw=0.5/(u-(u^0.5));
for i=1:(length(uu)-1)
Pw=Pw+A(end,i)*Y(:,i);
end
Pw=-Pw/A(end,end);

figure(2)
plot(T,Pw,'linewidth',2);
xlabel('Dimensionless time')
ylabel('Dimensionless wellbore pressure')

%% Actual Graphs
figure(3)
t=((T*(phi1*c1+phi2*c2)*vis*re^2)/k2)/3600; % in hours
p=(Pi-(Q*bo*vis*Y/(2*pi*k2*h)))/101325; %in atomspher
for i=1:(length(uu)-1)
    plot(t,p(:,i),'linewidth',2);
    hold on
end
hold off
xlim([0 max(t)])
xlabel('Time (hr)')
ylabel('Pressure (atm)')
% legend('1st collocation point','2ed collocation point','3th collocation point','4th collocation point','5th collocation point'...
%    ,'6th collocation point','7th collocation point','8th collocation point','9th collocation point','10th collocation point')
%% Wellbore pressure
figure(4)
pw=(Pi-(Q*bo*vis*Pw/(2*pi*k2*h)))*14.7/101325; %in atomspher
plot(t,pw,'linewidth',2);
xlim([0 max(t)])
xlabel('Time (hr)')
ylabel('Wellbore pressure (psi)')