clc, clear
%% Known/Chosen Values
ToD=input("What is the time of day in military time? (i.e. 6:30pm is 18:30, formatted as 18.50) ")
if ToD<6 || ToD>18
    error("Time of day is outside of time range")
else
    theta = ToD2theta(ToD);
end
Resolution = input("What is desired resolution? (m) ")
lambda = 0.0000005; %meters
wd = 3000 * 0.00006; %meters
Respon=10^(10); %V/W
MinDetec = 300; %mV/mu-meters
r = 0.3; %Reflectance of the ground
E_sun = 1370; %W/m^2
tau_optics = 0.96;
gamma = 1.4;
R = 287.05; %Gas Constant of Air in J/kg-K
%% Sensor Calculation
if Resolution <= 3.5
    da=input("What is diameter of aperture? (m) ")
    if da > 0.1
        error("Aperture diameter is too wide")
    elseif da<0
        error("Aperture diameter is too narrow")
    else
    end
    A_aperture = 0.25 * pi * (da^2); % compute appature area 
    H_da = res(Resolution, lambda); % Calculates the ratio of Altitude and Aperture Diameter given a desired resolution and chosen wavelength
    H = H_da*da %meters
    
    % display to the user what type of vehicle is based on altitude
    if H >= 160000
        disp("Vehicle is a satellite")
    elseif H<=20000
        disp("Vehicle is an aircraft")
    else 
        error("Vehicle not at functional height.")
    end

    f = input("What is the focal length? (m) ")
    if f > 10
        error("Focal length too large")
    elseif f < 1
        error("Focal length too small")
    end


    SW = Swath(H, wd, f) %meters
    if SW < 9656.064 %6 mi in meters
        error("Swath Width too narrow")
    end


    A_detected = 0.0006 * SW;
    Lr = ReflecRad(E_sun, r, theta)
    E_sensor = falloff(Lr, A_aperture,f, theta)
    S = signal(Respon,E_sensor, A_detected, theta)                            %NEED TO FIND A WAY TO INCORPORATE TAU_OPTICS INTO CODE

else
    error("Resolution is too high.")
end



%% Function List
function H_da = res(Resolution,lambda)  % Calculates the ratio of Altitude and Aperture Diameter given a desired resolution and chosen wavelength
H_da=Resolution/(1.22*lambda);
end

function SW = Swath(H,wd,f) % Calculates the Swath Width from the Altitude, detector width, and chosen focal length
SW=H*(wd/f);
end

function Lr = ReflecRad(E_sun, r, theta)    % Calculates the radience reflected from a subject using the Sun's irradience/exitance, the reflectance value of the subject, and the time of day in degrees theta
Lr=(E_sun*r*cosd(theta))/pi;
end

function E_sensor = falloff(Lr, A_aperture, f, theta)    % Calculates the Irradience entering the sensor from the Refected radiance, the area of the aperture, and the radius of the detector
E_sensor=(Lr*A_aperture*((cosd(theta))^4))/(f^2);   %For right now it is being assumed that r_0 from the lense falloff eq is = to focal length f
end

function S = signal(Respon, E_sensor, A_detected, theta)    % Calculates signal from the response time, sensor irradience, detected area, and time of day in degrees theta
S = Respon* E_sensor*A_detected*cosd(theta);
end

function V= v_ac(gamma, R, T)   % Converts atmospheric data to velocity of aircraft
V=sqrt(gamma*R*T);
end

function theta = ToD2theta(ToD) % Converts Time of Day to degrees theta
DegperHour=180/12; %conversion between hour and degrees
theta=ToD*DegperHour-180; %Converts the Time of Day (ToD) to degrees with 0 degrees being Noon (12:00)
end

function [T,p,rho]=Atmo(H)      % Atmospheric Model     
R=287.05; %Gas Constant of Air in J/kg-K
R_univ=8.3144598; %Universal gas constant of air in J/mol-K
    if H<=11000    %Troposphere
            T_SL= 288.15; % Temp @ Sea Level in Kelvin
            p_SL= 101325; % Pressure @ Sea Level in Pascals
        T=T_SL+((-0.0065)*H);  %Calculates the Temperature using the product of the height multiplied by the lapse rate and summing it with the base temp
        p= p_SL * (1-(0.0065/(T_SL))*(H-0))^((9.81*0.0289644)/(R_univ*0.0065));   %Calculates the pressure of the atmosphere using the Barometric Formula (assuming temp is variable)
        rho=p/(R*T);        %Density found using the Ideal Gas Law

    elseif 11000<H && H<=20000   %Tropopause
            p_TropPause=22632;   %base pressure of atmosphere in tropopause
            T_TropPause=216.65;    %base temp of atmosphere in tropopause
        T=T_TropPause-(0*H);
        p=p_TropPause* exp((-9.81*0.0289644*(H-11000))/(R_univ*T_TropPause));
        rho=p/(R*T);
    end
end