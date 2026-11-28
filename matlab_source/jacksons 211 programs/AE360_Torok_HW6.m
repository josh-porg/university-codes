%Homework 6 Jackson Torok

clc
clear
close all

deg=0;
m=[0;30;90;135;180;200;350].*pi/180;
e=[.01;.1;.5;.9];
tolrad=1E-6;
q=1;

for i=1:length(e)
    
    for j=1:length(m)
        
        if (m(j)<0 && m(j)>-pi) || m(j)>pi
        
            E0=m(j)-e(i);
            
        else
            
            E0=m(j)+e(i);
        
        end
        
        data(q,:)=[m(j)*180/pi,e(i),kep_new(E0,e(i),m(j),tolrad)*180/pi];
        
        q=q+1;
        
    end
    
end

table(data(:,1),data(:,2),data(:,3),'VariableNames',{'M (deg)','e','E (deg)'})

function Efinal=kep_new(E,e,m,tolrad)

l=2;

E(l)=E(l-1)+(m-E(l-1)+e*sin(E(l-1)))/(1-e*cos(E(l-1)));

l=l+1;

E(l)=E(l-1)+(m-E(l-1)+e*sin(E(l-1)))/(1-e*cos(E(l-1)));

    while abs((E(l)-E(l-1))-(E(l-1)-E(l-2)))>tolrad

       l=l+1;

       E(l)=E(l-1)+(m-E(l-1)+e*sin(E(l-1)))/(1-e*cos(E(l-1))); 

    end

    Efinal=E(l);
    
end


