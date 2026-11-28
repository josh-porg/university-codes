%Homework 9 Jackson Torok

clc
clear
close all

mu=3.98600441E5;

ri=[1131.34,-2282.343,6672.423];
vi=[-5.64305,4.30333,2.4289];
[e,a,i,om,w,ta]=convert(mu,ri,vi);
x=[ri';vi'];
tspan=86399;% a day in seconds
tolrad=1E-10;
P=2*pi*sqrt(a^3/mu);
k = floor(tspan/P);

if (ta>=0 && ta<pi)
    
    Ei=acos((e+cos(ta))/(1+e*cos(ta))); 
    q=0;
    
else
    
    Ei=acos((e+cos(ta))/(1+e*cos(ta)))+pi;
    q=1;
    
end

mi=Ei-e*sin(Ei);
m=mi+sqrt(mu/a^3)*tspan-2*pi*k;

e1=orb_euler_int(mu,x,tspan,.1)';
e2=orb_euler_int(mu,x,tspan,1)';
e3=orb_euler_int(mu,x,tspan,10)';
e4=orb_euler_int(mu,x,tspan,60)';
e5=orb_euler_int(mu,x,tspan,300)';

re1 = e1(1:3);
ve1 = e1(4:6);
re2 = e2(1:3);
ve2 = e2(4:6);
re3 = e3(1:3);
ve3 = e3(4:6);
re4 = e4(1:3);
ve4 = e4(4:6);
re5 = e5(1:3);
ve5 = e5(4:6);

Ef=kep_new(Ei,e,m,tolrad);

if (q==0)
    
    taf=acosd((cos(Ef)-e)/(1-e*cos(Ef)));
    
else
    
    taf=acosd((cos(Ef)-e)/(1-e*cos(Ef)))+pi;
    
end


[rfk,vfk]=convertback(mu,e,a,i,om,w,taf);

table(re1',rfk,ve1',vfk,'VariableNames',{'Euler R','Kepler R','Euler V','Kepler V'})

MEk=norm(vfk)^2/2-mu/norm(rfk);
AMk=norm(cross(rfk,vfk)); % assuming m=1

MEe=norm(re1)^2/2-mu/norm(ve1);
AMe=norm(cross(re1,ve1)); % assuming m=1

table(MEe,MEk,AMe,AMk,'VariableNames',{'Euler ME','Kepler ME','Euler |AM|','Kepler |AM|'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xf=orb_euler_int(mu,x,tspan,tstep)

    for i=1:tspan/tstep

        xdot=[x(4:6);-mu/(norm(x(1:3))^3)*x(1:3)];
        xn=x+xdot.*tstep;

        x=xn;

    end
    xf=x;
end

function [mage,a,i,om,w,ta]=convert(mu,r,v)

    e=1/mu*((norm(v)^2-mu/norm(r))*r-dot(r,v)*v);

    mage=norm(e);

    a=-mu/(norm(v)^2-(2*mu)/norm(r));

    h=cross(r,v);

    i=acosd(h(3)/(norm(h)));

    n=cross([0,0,1],h);

    if(n(2)>=0)
        om=acosd(n(1)/norm(n));
    else
        om=360-acosd(n(1)/norm(n));
    end

    if(e(3)>=0)
        w=acosd(dot(n,e)/(norm(n)*mage));
    else
        w=360-acosd(dot(n,e)/(norm(n)*mage));
    end

    if(dot(r,v)>=0)
        ta=acosd(dot(e,r)/(mage*norm(r)));
    else
        ta=360-acosd(dot(e,r)/(mage*norm(r)));
    end

end

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

function [R,V]=convertback(mu,mage,a,i,om,w,ta)

p=a*(1-mage^2);

Rpqw=(p/(1+mage*cosd(ta))).*[cosd(ta),sind(ta),0];

Vpqw=sqrt(mu/p).*[-sind(ta),mage+cosd(ta),0];

TransformMat=[cosd(om)*cosd(w)-sind(om)*sind(w)*cosd(i), -cosd(om)*sind(w)-sind(om)*cosd(w)*cosd(i), sind(om)*sind(i);...
              sind(om)*cosd(w)+cosd(om)*sind(w)*cosd(i), -sind(om)*sind(w)+cosd(om)*cosd(w)*cosd(i), -cosd(om)*sind(i);...
              sind(w)*sind(i),                           cosd(w)*sind(i),                            cosd(i)];

R=TransformMat*Rpqw';

V=TransformMat*Vpqw';

end
