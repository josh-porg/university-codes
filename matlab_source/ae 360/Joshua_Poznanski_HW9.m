%% Driver Final

close all; clear; clc;

% constants 
mu = 3.986004418*10^14; % in meters cubed per seconds squared
% setup
Ri = [1131.34; -2282.343; 6672.423]; % units: km
Vi = [-5.64305; 4.30333; 2.42879]; % units km/s

% convert km to m
Ri = Ri * 1000;
Vi = Vi * 1000;
%% a) euler integration of position and velocity
% TODO: use the given values for time
t = 24*60*60; % convert 24 hours to seconds (hr*60 = min, min*60 = s) % TODO: use the given values
dt = .1; % TODO: use the given values
[R,V,E,h] = integrate2Body3D(Ri,Vi,100000,dt);
%subplot(1,3,1); plot(R); title("R"); subplot(1,3,2); plot(V); title("V"); subplot(1,3,3); plot(E); title("E");%for debug
%figure; subplot(1,4,1); plot(R); title("R"); subplot(1,4,2); plot(V); title("V"); subplot(1,4,3); plot(E); title("E"); subplot(1,4,4); plot(h); title("h");
%figure; plot(h); title("h");

%% b) RV to COE, kepler's euqation, and COE to RV
[a,e,i,RA,w,nu] = rv2coe(Ri,Vi); %TODO: get create  version of this equation

e = norm(e); % we will use e as a scalar (it's magnitude) not as a vetor
n = sqrt(mu/a^3);
Ei = atan2(sqrt(1-e^2)*sin(nu), e+cos(nu)); % auto quadrent check
%Ei = acos(e+cos(nu) / (1+e*cos(nu))); % need manual quadrent check TODO: implement quad check if using this
Mi = Ei - e*sin(Ei);
Mf = Mi + n*t; %Mi+n*deltaT

%Eold = Mf; %by hand method
%Ef = Mf + e*sin(Eold); %by hand method
if (-pi<Mf && Mf<0 || Mf>pi)
    Eold = Mf - e;
else
    Eold = Mf + e;
end
Ef = Eold - (Mf-Eold+e*sin(Eold)) / (-1+e*cos(Eold)); %newton method
while(abs(Ef-Eold) > 1e-6)
    Eold = Ef;
    %Ef = Mf + e*sin(Ef); %by hand method 
    Ef = Eold - (Mf-Eold+e*sin(Eold)) / (-1+e*cos(Eold)); %newton method
end
Mf;
nuf = atan2(sqrt(1-e^2)*sin(Ef), cos(Ef)-e); % auto quadrent check
%by hand method nuf = 1.3285 ef = 1.3207
%newton raphson nuf = 1.3285 ef = 89.2853 mod 2pi = 1.3207

[Rn, Vn] = coe2rv(a,e,i,RA,w,nuf);

En = (1/2)*norm(Vn)^2-mu/norm(Rn);
hn = norm(cross(Rn,Vn));

%% c) comparison of a) and b)
Ei = (1/2)*norm(Vi)^2-mu/norm(Ri);



disp("comparison of final vlaues")
disp("euler")
disp("R = " + R(end))
disp("V = " + V(end))
disp("E = " + E(end))
disp("h = " + h(end))
disp("Newton")
disp("R = " + norm(Rn))
disp("V = " + norm(Vn))
disp("E = " + En)
disp("h = " + hn)
disp("initial E")
disp("E = " + Ei)
disp("change in E")
fprintf("Euler: %d\n",abs(Ei-E(end)));
fprintf("Newton: %d\n",abs(Ei-En));
fprintf("Newton is %2.1f percent more acurate than euler\n",abs(Ei-E(end))/abs(Ei-En))

fprintf("Differance in Position: %d\n",abs(R(end)-norm(Rn)));
fprintf("Differance in Velocity: %d\n",abs(V(end)-norm(Vn)));
fprintf("Differance in Specific Energy: %d\n",abs(E(end)-En));
fprintf("Differance in Specific Angular Momentum: %d\n",abs(V(end)-norm(Vn)));

function [R,V,E,h] = integrate2Body3D(Ri,Vi,t,dt)
    mu = 3.986004418*10^14;
    
    timesteps = t/dt;
    
    R = zeros(1,timesteps);
    V = zeros(1,timesteps);
    h = zeros(1,timesteps);
    X = [Ri(1); Ri(2); Ri(3); Vi(1); Vi(2); Vi(3)]; %X = [x; y; z; vx; vy; vz];
    %Xd = [V(1); V(2); V(3); ax, ay, az]; %Xd = [vx; vy; vz; ax, ay, az];
    
    for i = 1:timesteps
    %using the full calculation formula instead of norm (it is alegadly faster)
    R(i) = sqrt(X(1)^2 + X(2)^2 + X(3)^2); % magnitude of position
    V(i) = sqrt(X(4)^2 + X(5)^2 + X(6)^2); % magnitude of velocity
    h(i) = norm(cross([X(1),X(2),X(3)],[X(3),X(4),X(5)]));
    
    %caclulate change in state at current time step
    fxt = [X(4); X(5); X(6); -(mu/R(i)^3)*X(1); -(mu/R(i)^3)*X(2); -(mu/R(i)^3)*X(3)]; % fxt = [vx; vy; vz; -(mu/r^3)*x; -(mu/r^3)*y; -(mu/r^3)*z];
    X = X + + fxt*dt; %Xnew = Xold + fxt*dt; % advance time 
    end
    E = (1/2).*V.^2-mu./R;
end
