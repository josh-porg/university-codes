% Jackson Torok Exam 1
clear
close all
clc

%% Question 1:

%Code in the command window will be deleted if the window is closed,
%whereas an m-file will be saved and can be run anytime.

%% Question 2:

%1. Ending each line in a semicolon to suppress unwanted output.
%2. Commenting alongside your code to allow others to read it without
%coding knowledge and a large struggle
%3. Save often
%4. Keep the code simple and readable
%4.5. Indent and space correctly to aid in said readability

%% Question 3:

%1. A variable may not start with a number
%2. It cannot contain an operator
%3. It may not be the same as the program name
%4. It cannot exceed 64 characters

%% Question 4:

%Matlab reads code top to bottom in each section. For example: Matlab will
%read a program top-down, and if there is a function call, it will start
%reading that function top to bottom where ever it is, whether it is lower
%in the program or in a separate tab.
%Matlab follows order of operations in each line of code. A funtion call
%or array subsection call is treated as if it is in parentheses

%%Question 5:

%For matrix multiplication, the number of columns of the first matrix must
%match the number of rows in the second matrix. For element-wise
%multiplication, the size of both matrices must be the same. For matrix
%multiplication, Matlab uses * and /, for element-wise, Matlab uses .* and
%./  .

A=[2,6;5,4];
B=[7,9;1,2];

matmult=A*B

elemult=A.*B

%% Question 6:

aa=[1,2,3,9,8,7;4,5,6,6,5,4;7,8,9,3,2,1;4,5,6,6,5,4;7,8,9,3,2,1;1,2,3,9,8,7];

b1=aa(1:5,1)

b2=aa(1:5,3:6)

%% Question 7:

% P(4 of a kind)=number of ways to draw a 4 of a kind / total 4 card combinations

P = 13/nchoosek(52,4);

disp(['Chance that 4 cards drawn from a 52 card deck make a 4 of a kind:  ' num2str(P*100) ' %']);

%%Questtion 8:

v1=[3,2,4];

v2=[1,3,4];

theta=acosd(dot(v1,v2)/(norm(v1)*norm(v2)))

%% Question 9:

a=[2,8,15;25,9,32;14,28,51];

b=[16;17;38];

c=[7,19,34,12];

d=a(:,3)

e=dot(b,d)

f=cross(b,d)

g=sum(b+d)

h=horzcat(b,d)

i=e+c(1,3)+a(2,2)

j=det(a)

%% Question 10:

syms x y z

eq1=3*x+4*y-5*z==-4;

eq2=-x+2*y+z==6;

eq3=2*x-3*y+4*z==8;

sol=solve(eq1,eq2,eq3)

sol.x+sol.y+sol.z %I am assuming there was a typo and this should equal 6

%% Question 11:

g=9.8;
v=100;
thet=40;

tland=2*(v*sind(thet))/g

thet=[0:1:90];

[range,I]=max(((v*sind(thet))/g).*v.*cosd(thet))

thetamax=thet(I)

v=100;

t25=linspace(0,2*((v*sind(25))/g),200);
t45=linspace(0,2*((v*sind(45))/g),200);
t75=linspace(0,2*((v*sind(75))/g),200);

x25=t25.*v.*cosd(25);
y25=t25.*v.*sind(25)-.5.*g.*t25.^2;

x45=t45.*v.*cosd(45);
y45=t45.*v.*sind(45)-.5.*g.*t45.^2;

x75=t75.*v.*cosd(75);
y75=t75.*v.*sind(75)-.5.*g.*t75.^2;

plot(x25,y25);
hold on
plot(x45,y45);
plot(x75,y75);

xlabel('Distance (m)');
ylabel('Height (m)');

legend('25 Degrees','45 Degrees','75 Degrees');




