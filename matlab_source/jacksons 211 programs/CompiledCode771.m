% Written by: Jackson Torok %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close, clear, clc;

%% Declare variables

% "velocities"
syms v1 vt v2 c cstar M2 M Mt
% masses
syms mf m0 mpl mdot mprop mengine MRtot MRengine zeta
% engine dimensions
syms rt r2 dt d2 At A2
% pressures
syms p1 pt p2 p3 p
% Forces
syms F Fpress Fmoment 
% Time
syms t Isp Itot
% Misc
syms g0 wdot maxaccel accellim sf k


%% Define equations
% ideal rocket equations      note:c=v2 when fpress=0 -> Isp=v2/g0
eqs=[ c==v2+(p2-p3)*A2/mdot , mdot==wdot/g0 , F==c*mdot ...
    , Fpress==(p2-p3)*A2 , Fmoment==v2*wdot/g0 , mdot==(m0-mf)/t ...
    , At==pi/4*dt^2 , A2==pi/4*d2^2 , d2==2*r2 , dt==2*rt ...
    , MRengine==(mf-mpl)/(m0-mpl) , MRtot==mf/m0 , zeta==(m0-mf)/(m0-mpl) ...
    , maxaccel==F/mf , cstar==p1*At/mdot , Itot==F*t , accellim==sf*g0 , MRengine==mengine/mprop ...
    , mdot==(mprop-mengine)/t , Isp==c/g0 , p1==p*(1+1/2*(k-1)*M^2)^(k/(k-1)) , A2/At==Mt/M2*((1+1/2*(k-1)*M2^2)/(1+1/2*(k-1)*Mt^2))^((k+1)/(2*k-2))];

%% Set known values



%% Solve

% eval(eqs)

sol=solve(subs(eqs));

% for i=1:length(1)
%     sol(i,:)=struct2table(vpasolve(subs(eqs)));
% end

%% Output
disp(sol)

%% Plot
% 2 Y-axis figure
f=figure;
left_color = [.75 .25 0];
right_color = [0 .45 .45];
set(f,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(alt,sol.F/1000)
ylabel('Thrust (lb x10^3)')
xlabel('Altitude (ft)')
ylim([190 220])
yyaxis right
plot(alt,sol.Isp)
ylabel('Specific Impulse (sec)')
ylim([250 290])
title('Minuteman Rocket Performance')

%% Past Interesting Homeworks

% Problem 3 Knowns 
% alt=3.28084*[0,1000,3000,5000,10000,25000,50000,75000,100000,130000,160000,200000,300000,400000,600000,1000000];
% pr=[1,.887,.66919,.53313,.26151,.025158,.00078735,2.0408E-5,3.15983E-7,1.2341E-8,2.9997E-9,8.3628E-10,8.6557E-11,1.4328E-11,8.1056E-13,7.4155E-14];
% p2=8.66*144; %psf
% F=194600; %lbs
% A2=1642/144; %ft^2
% Isp=254; %sec
% g0=32.174; %ft/s^2   note: slug=lbm/g0
% c=Isp*g0;
% mdot=F/c;
% v2=(F-(p2-101325*pr(1)*0.020885)*A2)/mdot;
% syms Isp F c
% Interesting work-arounds:



