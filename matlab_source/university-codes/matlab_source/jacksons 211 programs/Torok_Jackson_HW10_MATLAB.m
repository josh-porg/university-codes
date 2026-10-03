%Homework 10 Jackson Torok

clc
clear
close all

%% 13.4

data=[1,2494,4157;2,1247,2078;3,831,1386;4,623,1039;5,499,831;6,416,693];

vol=data(:,1);

pres=data(:,2:3);

linguessfor4a = interp1(vol,pres,5.2)

cubguessfor4b = interp1(vol,pres,5.2,'cubic')

%% 13.7

figure('Name','13.7');

plot(vol,pres(:,1),'o')

hold on 

x=1:.2:6;

plot(x,polyval(polyfit(vol,pres(:,1),1),x),'-')

plot(x,polyval(polyfit(vol,pres(:,1),2),x),'-')

plot(x,polyval(polyfit(vol,pres(:,1),3),x),'-')%this one

plot(x,polyval(polyfit(vol,pres(:,1),4),x),'-')

%% 13.16

figure('Name','13.16');

time=0:24;

alt=[0,107.37,210,307.63,400,484.6,550,583.97,580,549.53,570,699.18,850,927.51,950,954.51,940,910.68,930,1041.52,1150,1158.24,1100,1041.76,1050];

subplot(2,2,1:2)
plot(time,alt);
title('Altitude')

subplot(2,2,3)
plot(time(2:25),diff(alt));
title('Velosity')

subplot(2,2,4)
plot(time(3:25),diff(alt,2)); % t=8,13,21
title('Acceleration')

%% 13.18

a=25.48;
b=1.522E-2;
c=-.7155E-5;
d=1.312E-9;

fun=@(T) a+b.*T+c.*T.^2+d.*T.^3;

deltah = quad(fun,300,1000)

%% 13.20

syms y(t)

cond=y(0)==1;

eq=diff(y,t)==1-sin(t);

y(t)=dsolve(eq,cond);

eval(y(4)-y(0))

