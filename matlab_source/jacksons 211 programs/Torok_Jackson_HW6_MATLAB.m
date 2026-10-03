%Homework 4 Jackson Torok

clc
clear
close all

%% 7.3

a=input('Area of the Base \n');

h=input('Height of Cone \n');

v=1/3*a*h;

disp(['The volume of the cone is: ', num2str(v)])

%% 7.6

age=input('\nWhat is your age?\n');

disp(['Your age is ', num2str(age)])

%% 7.12

yen=5:5:125;

y2d=yen*0.0092;

table(yen',y2d','VariableNames',{'Yen','Dollar'})

euro=1:2:60;

e2d=euro*1.19;

table(euro',e2d','VariableNames',{'Euro','Dollar'})

dol=1:10;

d2e=dol/1.19;

d2p=dol/1.39;

d2y=dol/0.0092;

table(dol',d2e',d2p',d2y','VariableNames',{'Dollar','Euro','Pound','Yen'})

%% 7.7

arr=input('Please input an array of numbers\n');

l=length(arr);

disp(['You entered a ', num2str(l), ' column long array'])

%% 7.13

last={'Smith','Jones','Webb','Anderson'};

first={'Fred','Kathy','Milton','John'};

Age=[6,22,92,45];

height=[47,66,62,72];

weight=[82,140,110,190];

table(last',first',Age',height',weight','VariableNames',{'Last Name','First Name','Age','Height','Weight'})

%% 7.17

t=0:pi/100:2*pi;

r=t./t;

[x,y]=pol2cart(t,r);

plot(x,y)

[px,py]=ginput(2);

hold on

plot(px,py)

dist=sqrt((px(2)-px(1))^2+(py(2)-py(1))^2)
