clear all

x2=[-1:.025:1];
y2=[-1:.025:1];

[X,Y]=meshgrid(x2,y2);

z=X.^coth(Y);

plot(z)
axis([-10,10,-10,10])