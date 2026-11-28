clear
clc
close all

%% 1a

t=[0:.5:5];

g=32.2;

y0=100;

y=y0-.5*g*(t.^2)

%% 1b

a=.6;
b=.5;
c=.7;
x0=-4;
y0=-2;
z0=-3;
xA=2;
yA=-3;
zA=1;

dA0=((xA-x0)^2+(yA-y0)^2+(zA+z0)^2)^(.5);

d=dA0*sin(acos(((xA-x0)*a+(yA-y0)*b+(zA-z0)*c)/(dA0*(a^2+b^2+c^2)^.5)))

%% 1c

orig=316501.673;

hund=round(orig,2);

thou=round(orig,-3)

%% 4.1

amat=[15, 3, 22; 3, 8, 5; 14, 3, 82];
bmat=[1; 5; 6];
cmat=[12, 18, 5, 2];

dmat=amat(:,3);

emat=horzcat(bmat,dmat);

fmat=vertcat(bmat,dmat);

gmat=vertcat(amat,cmat(1,1:3))

hmat=[amat(1,3), cmat(1,2), bmat(2,1)]

%% 4.3

times=[0:2:24];

Thermocouple=[84.3, 90, 86.7; 86.4, 89.5, 87.6; 85.2, 88.6, 88.3; 87.1, 88.9, 85.3; 83.5, 88.9, 80.3; 84.8, 90.4, 82.4; 85.0, 89.3, 83.4; 85.3, 89.5, 85.4; 85.3, 88.9, 86.3; 85.2, 89.1, 85.3; 82.3, 89.5, 89.0; 84.7, 89.4, 87.3; 83.6, 89.8, 87.2];   

tabmat=horzcat(times',Thermocouple);

[~,max1loc]=max(Thermocouple(:,1));

[~,max2loc]=max(Thermocouple(:,2));

[~,max3loc]=max(Thermocouple(:,3));

[~,min1loc]=min(Thermocouple(:,1));

[~,min2loc]=min(Thermocouple(:,2));

[~,min3loc]=min(Thermocouple(:,3));

maxtimes=[times(max1loc),times(max2loc),times(max3loc)]

mintimes=[times(min1loc),times(min2loc),times(min3loc)]

%% 4.5

ace_data = readtable('ace_data.dat');
ace_data = table2array(ace_data);

years=ace_data(:,1);

ace=ace_data(:,2);

tropical_storms=ace_data(:,3);

hurricanes=ace_data(:,4);

major_hurricanes=ace_data(:,5);

[~,highestace]=max(ace);

[~,mosttrop]=max(tropical_storms);

[~,mosthurr]=max(hurricanes);

[~,mostbighurr]=max(major_hurricanes);

maxyears=[years(highestace), years(mosttrop), years(mosthurr), years(mostbighurr)]

means=[mean(ace), mean(tropical_storms), mean(hurricanes), mean(major_hurricanes)]

medians=[median(ace), median(tropical_storms), median(hurricanes), median(major_hurricanes)]

sorted=sortrows(ace_data,2,'descend')

%% 4.8

rho=[13560, 1000];

g=9.81;

pressure=[0:10000:100000]';

height=pressure./(rho.*g) %mercury is the first column

%% 4.10

amatnull=zeros(size(amat));
bmatnull=zeros(size(bmat));
cmatnull=zeros(size(cmat));

%% 4.11

mag=magic(6);

rowsum=sum(mag,2)

colsum=sum(mag)

diasum=[sum(diag(mag)),sum(diag(flip(mag)))]


