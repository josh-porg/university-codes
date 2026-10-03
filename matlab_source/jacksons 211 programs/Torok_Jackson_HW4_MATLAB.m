%Homework 4 Jackson Torok

clc
clear
close all

%% 5.14

theta=[pi/2;pi/4;pi/6];
v=100;
g=9.8;
t=[0:.1:20];

x=t.*v.*cos(theta);

y=t.*v.*sin(theta)-.5.*g.*t.^2;

figure('Name','Question 5.14');

subplot(2,2,1);

plot(t,x(2,:))
title('Horizontal Distance (m) vs Time (s)')
xlabel('Time (s)');
ylabel('Horizontal Distance (m)');
grid on

subplot(2,2,2);

plot(t,y(2,:))
title('Vertical Distance (m) vs Time (s)')
xlabel('Time (s)');
ylabel('Vertical Distance (m)');
grid on

subplot(2,2,3);

plot(x(2,:),y(2,:))
title('Vertical Distance (m) vs Horizontal Distance (m)')
xlabel('Horizontal Distance (m)');
ylabel('Vertical Distance (m)');
grid on

subplot(2,2,4);

plot(x(1,:),y(1,:),'--b',x(2,:),y(2,:),'-g',x(3,:),y(3,:),':k')
title('Vertical Distance (m) vs Horizontal Distance (m)')
xlabel('Horizontal Distance (m)');
ylabel('Vertical Distance (m)');
grid on
legend('\theta=\pi/2','\theta=\pi/2','\theta=\pi/2')

%% 5.15

thet=linspace(0,2*pi,101);

figure('Name','Question 5.15');

subplot(2,2,1);
r=sin(thet).^2 + cos(thet).^2;
polarplot(thet,r)
grid on

subplot(2,2,2);
r=sin(thet);
polarplot(thet,r)
grid on

subplot(2,2,3);
r=exp(thet./5);
polarplot(thet,r)
grid on

subplot(2,2,4);
r=sinh(thet);
polarplot(thet,r)
grid on

%% 5.20

Q=1000;
k0=10;
R=8.314;

T=[300:1:1000];

k=k0.*exp(-Q./(R.*T));

figure('Name','Question 5.20');

subplot(2,1,1);
plot(T,k)
title('Reaction Rate k vs Temperature (K)')
xlabel('Temperature (K)');
ylabel('k');
grid on

subplot(2,1,2);
plot(1./T,log(k))
title('Log(k) vs 1/T')
xlabel('1/T');
ylabel('Log(k)');
grid on

%% 5.21

G=[68,83,61,70,75,82,57,5,76,85,62,71,96,78,76,68,72,75,83,93];

figure('Name','Question 5.21');

subplot(2,2,1);
bar(sort(G))
title('Part a')
subplot(2,2,2);
histogram(G)
title('Part b')
subplot(2,2,3);
histogram(G,'BinEdges',[0,60,60,70,70,80,80,90,90,100])
title('Part c')
subplot(2,2,4);
histogram(G,'BinEdges',[0,60,60,70,70,80,80,90,90,100], 'Normalization','countdensity')
title('Part d')

%% 5.22

[n,~]=histcounts(G,[0,60,70,80,90,100]);

figure('Name','Question 5.22');

subplot(2,1,1);
pie(n)
legend('E','D','C','B','A','Location','eastoutside')
subplot(2,1,2);
pie3(n)
legend('E','D','C','B','A')

%% 5.26

horizontal=x(2,:);
vertical=y(2,:);

figure('Name','Question 5.26');

yyaxis left
plot(t,horizontal)
ylabel('Distance (m)')
yyaxis right
plot(t,vertical)
ylabel('Height (m)')
xlabel('Time (s)')
title('Distance (m) and Height (m) vs Time (s)')

%% 5.29

x1=[0:pi/100:20*pi];

y1=x1.*sin(x1);

z=x1.*cos(x1);

figure('Name','Question 5.29');

subplot(2,2,1);
plot(x1,y1)
xlabel('X')
ylabel('Y')
title('X vs Y')

subplot(2,2,3:4);
polarplot(x1,y1)
title('X vs Y Polar')

subplot(2,2,2);
plot3(x1,y1,z)
xlabel('X')
ylabel('Y')
zlabel('Z')
title('X vs Y vs Z')

%% 5.31

x2=[-5:.5:5];
y2=[-5:.5:5];

[X,Y]=meshgrid(x2,y2);

Z=sin(sqrt(X.^2+Y.^2));

figure('Name','Question 5.31');

subplot(2,3,1);
mesh(Z)
title('Mesh')

subplot(2,3,2);
surf(Z)
title('Surf One Input')

subplot(2,3,3);
surf(X,Y,Z)
title('Surf Three Inputs')

subplot(2,3,4);
surf(Z,'FaceLighting','gouraud')
colormap('jet')
title('Surf New Shading')

subplot(2,3,5);
[C,H]=contour(Z);
clabel(C,H);
title('Contour')

subplot(2,3,6);
surf(Z,'FaceAlpha',.5)
hold on
[C,H]=contour(Z);
clabel(C,H);
title('Surface and Contour Combination')






