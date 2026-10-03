clear
clc
close all

disp(['Question 3:']);

r=[-7953.8073,-4174.5370,-1008.9496];
v=[3.6460035,-4.9118820,-4.9193608];

mu=3.98600441E5;

e=1/mu*((norm(v)^2-mu/norm(r))*r-dot(r,v)*v)

mage=norm(e)

a=-mu/(norm(v)^2-(2*mu)/norm(r))

h=cross(r,v)

i=acosd(h(3)/(norm(h)))

n=cross([0,0,1],h)

if(n(2)>=0)
    om=acosd(n(1)/norm(n))
else
    om=360-acosd(n(1)/norm(n))
end

if(e(3)>=0)
    w=acosd(dot(n,e)/(norm(n)*mage))
else
    w=360-acosd(dot(n,e)/(norm(n)*mage))
end

if(dot(r,v)>=0)
    ta=acosd(dot(e,r)/(mage*norm(r)))
else
    ta=360-acosd(dot(e,r)/(mage*norm(r)))
end


%% Question 4

a=a;  %This is so if other values want to be tested it is easier to set
i=i;
om=om;
w=w;
ta=round(ta);
mu=mu;
mage=mage;

p=a*(1-mage^2)

Rpqw=(p/(1+mage*cosd(ta))).*[cosd(ta),sind(ta),0]

Vpqw=sqrt(mu/p).*[-sind(ta),mage+cosd(ta),0]

TransformMat=[cosd(om)*cosd(w)-sind(om)*sind(w)*cosd(i), -cosd(om)*sind(w)-sind(om)*cosd(w)*cosd(i), sind(om)*sind(i);...
              sind(om)*cosd(w)+cosd(om)*sind(w)*cosd(i), -sind(om)*sind(w)+cosd(om)*cosd(w)*cosd(i), -cosd(om)*sind(i);...
              sind(w)*sind(i),                           cosd(w)*sind(i),                            cosd(i)]

R=TransformMat*Rpqw'

V=TransformMat*Vpqw'




