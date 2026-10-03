%Homework 8 Jackson Torok

clc
clear
close all

fileID = fopen('report.txt','w');

%% 9.5

p=primes(100);
mult=zeros(1,length(p)-1);

for i=1:length(p)-1
    
    mult(i)=p(i)*p(i+1);
    
end

%% 9.6

q=figure('Name','Question 9.6','Visible','off');

f=input('Input the first 2 numbers in the Fibbonachi sequence as a horizontal vector\n');
count=input('How many terms do you want\n');

for i=3:count
    
   f(i)=f(i-2)+f(i-1);
    
end

i=1:i;

polarplot(i,f(:))

set(q,'Visible','on')

%% 9.7

w=figure('Name','Question 9.7','Visible','off');

f=input('Input the first 2 numbers in the Fibbonachi sequence as a horizontal vector\n');
count=input('How many terms do you want\n');

i=3;

while i<=count
    
   f(i)=f(i-2)+f(i-1);
   
   i=i+1;
    
end

i=1:i-1;

figure(w);
polarplot(i,f(:))

set(w,'Visible','on')

%% 9.8

x=input('Input the first 2 numbers in the Fibbonachi sequence as a horizontal vector\n');

k=3;

x(k)=x(k-2)+x(k-1);

while abs((x(k)/x(k-1))-(x(k-1)/x(k-2)))>.001
    
   k=k+1;
   
   x(k)=x(k-2)+x(k-1);

end

    phi=x(k)/x(k-1)
    
    error=(x(k)/x(k-1)-x(k-1)/x(k-2))

%% 9.13

in=5;

mycalc=my_sqrt(in,2,.0001)

myerror=abs(my_sqrt(in,2,.0001)-sqrt(in))

%% 9.15

mycalc=my_sin(2)

myerror2=abs(my_sin(2)-sin(2))

%% 9.16

l=1;
picalc=1;

while l<=3000
    
    l=l+1;

    picalc(l)=picalc(l-1)+((-1)^(l-1))/(2*l-1);

    if (abs(picalc(l)-picalc(l-1))<.001)
        
        break
        
    end

end

my_pi=4*picalc(l)

myerror3=my_pi-pi

%% 9.18

data=[3590.66,3614.17,3622.14,3620.55,3636.91,3604.42,3578.69,3593.57;
      3590.66,3612.05,3620.16,3614.95,3635.28,3601.47,3575.55,3592.23;
      3589.77,3610.43,3619.41,3610.73,3635.33,3598.96,3574.76,3591.02;
      3594.09,3611.26,3620.50,3611.93,3635.76,3596.53,3577.56,3590.18;
      3610.81,3629.09,3625.96,3623.13,3636.83,3599.44,3589.38,3597.27;
      3631.05,3640.49,3638.82,3648.98,3633.90,3600.07,3609.19,3613.54;
      3633.00,3641.14,3636.52,3660.86,3628.45,3594.17,3608.05,3612.62;
      3629.55,3637.50,3634.55,3655.34,3623.62,3589.64,3605.82,3609.07;
      3626.90,3635.37,3633.66,3653.01,3621.56,3591.25,3605.53,3606.01;
      3623.82,3633.52,3634.08,3650.27,3619.46,3590.88,3605.57,3606.44;
      3621.90,3631.10,3630.31,3645.67,3615.10,3587.90,3601.87,3605.47;
      3617.89,3626.22,3626.54,3639.75,3609.82,3584.43,3597.75,3600.80];

Year=[2008,2009,2010,2011,2012,2013,2014,2015];

avgperyear=[mean(data(:,1)),mean(data(:,2)),mean(data(:,3)),mean(data(:,4)),mean(data(:,5)),mean(data(:,6)),mean(data(:,7)),mean(data(:,8))];

avg=mean(avgperyear);

exceedperyear=zeros(1,length(Year));

for o=1:size(data,2)

    for p=1:size(data,1)
        
        if (data(p,o)>avg)
            
            fprintf(fileID,'Month %2u in %4u went above the 8 year average\n',p,Year(o));
            
            exceedperyear(o)=exceedperyear(o)+1;
            
        end
    end
end

exceedperyear

fclose(fileID);
%% functions

function out=my_sqrt(num,guess,error)

    k=1;
    g(k)=guess;
    
    y(k)=1/num*g(k)^2;
    g(k+1)=g(k)/8*(15-y(k)*(10-3*y(k)));

    while abs(g(k+1)-g(k))>=error
        
        k=k+1;
        
        y(k)=1/num*g(k)^2;
        
        g(k+1)=g(k)/8*(15-y(k)*(10-3*y(k)));
        
    end
    
    out=g(k+1);

end

function out=my_sin(x)

    k=2;

    s=[x,x-(x^3)/factorial(3)];

    while true
        
        if (abs(s(k)-s(k-1))<=.001)
            break
        end
        
        k=k+1;
        
        s(k)=s(k-1)+((-1)^(k-1))*(x^(2*k-1))/factorial(2*k-1);
        
    end

    out=s(k);
    
end
