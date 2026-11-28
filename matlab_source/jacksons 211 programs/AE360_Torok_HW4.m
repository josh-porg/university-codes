%Homework  4 Jackson Torok

clc
clear
close all

mu=3.986E5;
e=.1;
a=7000;%km

tstep=.01;%second(s)
tspan=86399;% a day in seconds

data=zeros(tspan/tstep+1,4);

rp=a*(1-e^2)/(1+e);%km
vp=sqrt(mu*(2/rp-1/a));%km/s

x=[rp;0;0;vp];
E=norm(x(3:4))^2/2-mu/norm(x(1:2));
data(1,:)=[norm(x(1:2)),norm(x(3:4)),E,0];

for i=1:tspan/tstep
    
    xn=orb_euler_int(x,tstep);
    E=norm(x(3:4))^2/2-mu/norm(x(1:2));
    data(i+1,:)=[norm(xn(1:2)),norm(xn(3:4)),E,i*tstep];
    
    x=xn;
    
end

figure('Name','Position')
plot(data(:,4)./86400,data(:,1))
title('Magnatude of Position (km) Vs Time (days)')
xlabel('Time (days)')
ylabel('Magnatude of Position (km)')

figure('Name','Velosity')
plot(data(:,4)./86400,data(:,2))
title('Magnatude of Velosity (km/s) Vs Time (days)')
xlabel('Time (days)')
ylabel('Magnatude of Velosity (km/s)')

figure('Name','Specific Mechanical Energy')
plot(data(:,4)./86400,data(:,3))
title('Specific Mechanical Energy (km^2/s^2) Vs Time (days)')
xlabel('Time (days)')
ylabel('Specific Mechanical Energy (km^2/s^2)')



function xnext=orb_euler_int(x,t)
mu=3.986E5;

xdot=[x(3:4);-mu/(norm(x(1:2))^3)*x(1:2)];

xnext=x+xdot*t;
end

