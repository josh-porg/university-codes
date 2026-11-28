%Homework 3 Jackson Torok

clc
clear
close all

% Question 1

F1=[10,5,6];
F2=[-8,3,9];

CosTheta = max(min(dot(F1,F2)/(norm(F1)*norm(F2)),1),-1);
ThetaInDegrees = real(acosd(CosTheta))

%% Question 2

syms x y z
a=3*x+4*y+8*z==15;
b=-x+y+4*z==23;
c=-5*x+4*y+7*z==5;

sol=solve(a,b,c);

d=7*sol.x(1)-2*sol.y(1)+sol.z(1) %The system does not work as this output differs from the expected


%% Question 3

m=[15,25,10,50,30];
loc=[10,15,5;15,5,-10;20,25,15;5,5,5;25,35,20];
C=[dot(m',loc(:,1)),dot(m',loc(:,2)),dot(m',loc(:,3))]./sum(m)

%% Question 4

L=10/12;
F1=750;
F2=150;
F3=500;
theta1=30;
theta2=120;
theta3=190;

Ma=L*(F1*sind(theta1-theta1)+F2*sind(theta2-theta1)+F3*sind(theta3-theta1))

%% 10.14

A=[2,-1;2,5];
B=[4,2;2,1];
C=[2,0,0;1,2,2;5,-4,0];

deta=det(A)
detb=det(B)
detc=det(C)

if deta~=0

    inva=inv(A)
    
end

if detb~=0

    invb=inv(A)
    
end

if detc~=0

    invc=inv(A)
    
end

%% 10.16

l=10/12;
h=5/12;
F=35;
thetaFor=55;

len=sqrt(l^2+h^2);

thetaPos=180-atand(h/l);

Mbracket=len*F*sind(thetaPos-thetaFor)

%% 10.18

amat=[-2,1;1,1];
aansmat=[3;10];
bmat=[5,3,-1;3,2,1;4,-1,3];
bansmat=[10;4;12];
cmat=[3,1,1,1;1,-3,7,1;2,2,-3,4;1,1,1,1];
cansmat=[24;12;17;0];

aleftans=amat\aansmat
ainvans=inv(amat)*aansmat

bleftans=bmat\bansmat
binvans=inv(bmat)*bansmat

cleftans=cmat\cansmat
cinvans=inv(cmat)*cansmat

%% 5.1

eq1=@(x) exp(x);
eq2=@(x) sin(x);
eq3=@(x) 5*x^2+2*x+4;
eq4=@(x) sqrt(x);

subplot(2,2,1);

fplot(eq1)
title('exp(x)')
xlabel('X');
ylabel('Y');
axis([0,10,0,10000])
grid on

subplot(2,2,2);

fplot(eq2)
title('sin(x)')
xlabel('X');
ylabel('Y');
axis([0,10,-1,1])
grid on

subplot(2,2,3);

fplot(eq3)
title('5*x^2+2*x+4')
xlabel('X');
ylabel('Y');
axis([0,10,0,500])
grid on

subplot(2,2,4);

fplot(eq4)
title('sqrt(x)')
xlabel('X');
ylabel('Y');
axis([0,10,0,4])
grid on

%% 5.4

equation1=@(x) sin(x);
equation2=@(x) sin(2*x);
equation3=@(x) sin(3*x);

figure(2)

fplot(equation1,'--r')

hold on

fplot(equation2,'-b')

fplot(equation3,':g')

axis([-pi,pi,-1,1])

