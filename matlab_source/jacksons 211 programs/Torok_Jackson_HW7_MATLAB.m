%Homework 7 Jackson Torok

clc
clear
close all

%% 8.3

data=[1530,116,45,110;1240,114,42,115;2380,118,41,120;1470,124,38,95;3590,126,61,118];

i=find(data(:,2)>115&data(:,2)<125);

fprintf('Engines %u Passed Temp \n',data(i,1))

i=find(data(:,3)>40&data(:,3)<60);

fprintf('Engines %u Passed Humidity \n',data(i,1))

i=find(data(:,4)>100&data(:,4)<200);

fprintf('Engines %u Passed Pressure \n',data(i,1))

i=find(data(:,2)>115&data(:,2)<125&data(:,3)>40&data(:,3)<60&data(:,4)>100&data(:,4)<200);

fprintf('Engines %u Passed Every Test \n',data(i,1))

totpass=100*length(i)/length(data);

fprintf('The Total Pass Rate Was: %4.2f%%\n',totpass)

%% 8.6

figure('Name','Question 8.6')

fplot(@g,[-2*pi 2*pi])

%% 8.9

figure('Name','Question 8.9')

data1=[1979,7.19,16.48;1980,7.83,16.15;1981,7.24,15.65;1982,7.44,16.17;1983,7.51,16.13;1984,7.10,15.65;1985,6.91,16.09;1986,7.53,16.10;1987,7.47,15.99;1988,7.48,16.16;1989,7.03,15.54;1990,6.23,15.90;1991,6.54,15.52;1992,7.54,15.50;1993,6.50,15.90;1994,7.18,15.62;1995,6.12,15.35;1996,7.87,15.16;1997,6.73,15.61;1998,6.55,15.69;1999,6.23,15.45;2000,6.31,15.30;2001,6.74,15.64;2002,5.95,15.46;2003,6.3,15.52;2004,6.04,15.08;2005,5.56,14.77;2006,5.91,14.45;2007,4.29,14.66;2008,4.72,15.27;2009,5.38,15.16;2010,4.92,15.14];

plot(data1(:,1),data1(:,2))
hold on
plot(data1(:,1),data1(:,3))
plot([1979,1979;2010,2010],[15.52,6.51;15.52,6.51])
legend('September','March','March Average','September Average','Location','Southeast')
axis([1979,2010,0,20])

data1(find(data1(:,2)>6.51),1)

data1(find(data1(:,3)>15.52),1)

%% 8.12

x=input('Enter x value: ');
y=input('Enter y value: ');

if (x>y)
    
    disp('x>y')
    
else
    
    disp('x<=y')
    
end

%% 8.13

a=input('Enter a number for Arcsin:');

arcsin(a);

%% 8.18

cr = [130 130 122 126.5 129];
m = menu('What is your major?', 'Civil Engineering', 'Chemical Engineering', 'Computer Engineering', 'Electrical Engineering', 'Mechanical Engineering');

switch m
    
case 1
output=130;

case 2
output = 130;

case 3
output = 122;

case 4
output = 126.5;

case 5
output = 129;

end

fprintf('You''ll need %5.1f credits to graduate \n',output)

%% 8.19

p = menu('What type star would you like?','Five Points','Six Points')
switch p
    
case 1
theta=pi/2:4*pi/5:4.5*pi;
r = ones(1,length(theta));
figure;
polarplot(theta,r)

case 2
theta=pi/2:2*pi/3:2*pi+pi/2;
r=ones(1,length(theta));
figure;
polarplot(theta,r)
hold on
theta = pi/6:2*pi/3:2*pi+pi/6;
polarplot(theta,r)
hold off

end

%% 8.20

l = menu('Did you stay in the long term or short term lot?','Long','Short');

if l ==1
    
    w = input('How many weeks?');
    d = input('How many days?');
    h = input('How many hours? Round up');
    
    hc = 2 + (h-1)*1;
    
    if hc >9
        
        hc = 9;
        
    end
    
        dc = d *9;
        
    if dc>60
        
        dc = 60;
    end
    
    total = w*60 + dc + hc
    
else
        
    d = input('How many days? ');
    h = input('How many hours? ');
    m = input('How many minutes? ');
    
    dc = d*32;
    
    m = h*60+m;
    
    mc = 2 + (ceil((m-30)/20)*1);

    if mc<2

        mc=2;

    elseif mc>32
        
        mc = 32;
        
    end
    
    total = dc+mc
    
end


%% Functions


function out=g(x)

    if x<-pi
       out=-1; 
    elseif x>=-pi && x<=pi
       out=cos(x);
    elseif x>pi
       out=-1;
    end
    
end

function out=arcsin(x)

    if (x<1 && x>-1)
        out = asin(x)
    else
        out = 'Arcsin is non-real at the current input'
    end

end

