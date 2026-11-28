clear all

%2.3%

a=5^2
b=(5+3)/(5*6);
c=sqrt(4+6^3)
d=9+6/12+7*5^(2+3)
e=1+5*3/6^2+2^(2-4)/1/5.5

%%
%2.4%

r=5;
Area=pi*r^2

r=10;
SurfArea=4*pi*r^2

r=2;
Volume=4/3*pi*r^3

%%
%2.8%

r=3;
h=[1,5,12];

VolCyl=pi*r^2*h

h=12;
b=[2,4,6];
AreaTri=1/2*b*h
VolTriPrism=AreaTri*6

%%
%2.12%

v=linspace(1,20,20);
v=linspace(0,2*pi,21);
v=linspace(4,20,15)
v=logspace(1,3,10)

%%
%2.14%

g=9.8;
t=linspace(0,100,101);
d=1/2*g*t.^2;
tab=table(t',d','VariableNames',["Time","Distance"])

%%
%2.17%

G=6.673E-11;
Me=6E24;
Mm=7.4E22;
d=3.9E8;

F=G*Me*Mm/(d^2)

%%
%3.4%

k0=1200/60;
Q=8000;
R=1.987;
T=linspace(100,500,9);

k=k0*exp(-Q./(R.*T));

tab1=table(T',k','VariableNames',["Temperature (K)","Rate Constant k (s^-1)"])

%%
%3.10%

theta=round(linspace(0,2*pi,64),1);

tab2=table(theta',sin(theta'),cos(theta'),tan(theta'),'VariableNames'...
     ,["Theta (rad)","sin(theta)","cos(theta)","tan(theta)"])

%%
%3.17%

G=[68,83,61,70,75,82,57,5,76,85,62,71,96,78,76,68,72,75,83,93];
size=length(g)
G=sort(G)
avg=mean(G)
med=median(G)
mod=mode(G)
stddev=std(G)
%the median is more representitive of the data because there are outliers
%that distort the mean while leaving the median unaffected

%%
%3.21%

norm=normrnd(70,2,241,1);
t=linspace(0,120,241);
plot(t,norm);

maxtemp=max(norm);
mintemp=min(norm);

maxtemptime=t(norm==maxtemp);
mintemptime=t(norm==mintemp);
disp(['Max Temp: ',num2str(maxtemp),' at: ',num2str(maxtemptime),' minute(s)']);
disp(['Min Temp: ',num2str(mintemp),' at: ',num2str(mintemptime),' minute(s)']);

%%
%3.23%

R=5;
f=15E3;
w=2*pi*f;
v=10;
c=1E-9;
l=200E-3;

Zc=1/(w*c*1j);
Zl=w*l*1j;
Z=Zc+Zl+R;

I=v/Z;
mag=abs(I);
ang=angle(I);
