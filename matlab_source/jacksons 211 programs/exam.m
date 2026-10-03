
% Final Exam Jackson Torok
clc
clear
close all
%% 1
% 1. Save regularly
% 2. Comment and explain code alongside
% 3. Name variables and functions acording to their functionality
% 4. Avoid global varibles when possible due to the ease of overwriting
% accidents and difficulty of debugging
%% 2
%The higher order approximations tend to become exponentially worse the
%further from the intital point that you get.
%% 3
r=input('Please input a radius for a circle: ');
Area = pi*r^2
%% 4
zeros(2)
ones(2)
eye(2)
%% 5
% interpolation takes the data points and approximates a line that
% intersects all of them, allowing you to approximate an intermediate
% value. Regression tries to find a best fit function that shows the
% approximate relationship between the x and y data sets.
% Interpolation because it hits every data point exactly
%% 6
%(a < = x & b > = 2 & i > 1 | (b * x ~ = 0))
%((-3<=0 & 2>=2) & (1>1 | (2*0 ~= 0 )
%((true & true) & (false | false)
%(true & false)
%
% it is false
a=-3;
b=2;
x=0;
i=b/2;
logans = (a<=x & b>= 2 & i > 1 | (b * x ~= 0)); %verification
%% 7
A=[2 4 6; 1 3 5; 7 8 9];
B = [
'teacher'; 'student'
];
C = 1:3;
D = {A B C};
size = size(D); % ans a it is a 1x3 cell
paracall = D(1,2);
curlycall = D{1,2};
% calling with () will get you a 1x1 cell with the corresponding
% information. Calling with {} will get you the array within that cell.
dfirst = D{1,1};
dthird = D{1,3};
c7 = dfirst(2,3)+dthird(2)
%% 8
if (x<0)
fx=-x^2;
elseif (x>=0 && x<=2)
fx=x^2;
elseif (x>=2 && x<=6)
fx=x^3;
else
fx=x;
end
%% 9
pimat=pi*ones(4,6)
%% 10
[t,y] = ode45(@odefun,[0 1],[0 0])
%% 11
% The first john is a string which is handled as singular value, whereas
% the second one is a char array which is handled as a vector of individual
% characters. You can put them in an array together, but the char version
% will be converted into a string. If you want to keep each variable's
% format, you will need to put them in a cell array.
%% 12
a12=fibonacci(1,1,10)
b12=fibonacci(5,8,10)
c12=a12.*b12
%% 13
apoly=0;
q=2;
piguess=0;
tol=1E-6;
while (pi-piguess>tol)
x=linspace(-1,1,q);
y=sqrt(1-x.^2);
apoly=trapz(x,y);
piguess = 2*apoly;
q=q+1;
end
piguess = vpa(piguess)
q %this is the aproximate number of sides it must have
%% 14
% this would have been 10 times easier just using the erf() function,
% oh well
x=linspace(0,3,100);
u=1;
for i=x
errf(u)=(2/sqrt(pi))*integral(@(s) exp(-(s.^2)),0,i);
u=u+1;
end
plot(x,errf(:))
hold on
plot(x,gradient(errf(:),3/100))
ttl=title('Error Function');
xl=xlabel('Values of x');
yl=ylabel('erf(x)');
grid on
erftl=text(2.5,1.05,'erf(x)');
ddxtl=text(1.5,.15,'d/dx');
ttl.FontName='Times New Roman';
erftl.FontName='Times New Roman';
ddxtl.FontName='Times New Roman';
xl.FontName='Times New Roman';
yl.FontName='Times New Roman';
ttl.FontSize=16;
figure(2)
p=polyfit(x,errf(:),3);
plot(x,errf(:))
hold on
plot(x,polyval(p,x))
ttl=title('Error Function');
xl=xlabel('Values of x');
yl=ylabel('erf(x)');
legend('erf(x)','Polynomial Approx','Location','southeast')
grid on
ttl.FontName='Times New Roman';
xl.FontName='Times New Roman';
yl.FontName='Times New Roman';
ttl.FontSize=16;
% Answer 15 is in the other 2 files
%% functions
function dydt = odefun(t,y)
dydt=t^2+y;
end
function arr=fibonacci(fir,sec,num)
tem=[fir,sec];
for (i=3:(num))
tem(i)= tem(i-1)+tem(i-2);
end
arr=tem;
end