%Exam 2 Jackson Torok

clc
clear
close all

%% 1

% The function declaration, name, output, and possible inputs

% Global variables are easy to create and allow comunication between
% programs, but they can easily be rewritten which  makes debugging
% difficult

% The name of a function reflects that of a variable. They must start with 
% a letter, and can contain letters, digits, or underscores. No operators
% or paranthesis are allowed.

%% 2

% ==  equal to if the two sides are equal, then the output is true 
% ~=  not equual to  if the two sides are not eequal, then the output is true

% true & false ---> False

%% 3

figure('Name','Question 3')

x=linspace(-2,2,100);

y=linspace(-2,2,100);

[X,Y]=meshgrid(x,y);

Z=sin(3.*X+2.*Y).*exp(-X.^2-Y.^2);

subplot(2,2,1)
mesh(X,Y,Z)
xlabel('x')
ylabel('y')
zlabel('z')
title('Mesh')

subplot(2,2,2)
surf(X,Y,Z)
xlabel('x')
ylabel('y')
zlabel('z')
title('Surface')

subplot(2,2,3)
hold on
surfc(X,Y,Z)
view([30,35])
xlabel('x')
ylabel('y')
zlabel('z')
title('Surface and Contour')

subplot(2,2,4)
contour(X,Y,Z)
xlabel('x')
ylabel('y')
title('Contour')

%% 4

data=[2001,20,19;2002,23,24;2003,19,25;2004,30,29;2005,25,30;2006,15,40;2007,32,29];

years=data(:,1);

in_students=data(:,2);

out_students=data(:,3);

temp=in_students-out_students;

i=find(temp>0);

for k=1:length(i)
    fprintf('In year %4u, there are %2u in-state students and %2u out-of-state students \n',years(i(k)),in_students(i(k)),out_students(i(k)))
end

%% 5

figure('Name','Question  5')

x=linspace(0,10,40);
f=zeros(length(x),1);

for i=1:length(x)

    f(i)=q5(x(i));

end

plot(x,f)
xlabel('x')
ylabel('F(x)')
title('F(x) vs x')

%% 6

x1=pi/6;
real=cos(pi/6);
guess=0;
error=real-guess;
iter=0;

while abs(error)>1E-5

    guess=guess+(((-1)^iter)*((x1)^(2*iter))/(factorial(2*iter)));

    error=real-guess;

    iter=iter+1;

end

fprintf('\n \nThe result is %8.6f with an error of %3.3e taking %u iterations',guess,error,iter)

%% function for 5

function f=q5(x)

    if x<0
    f=0;
    elseif x>=0&&x<=6
    f=x.^2;
    elseif x>6&&x<=8
    f=36;
    elseif x>8   
    f=4.5.*x; 
    end
    
end




