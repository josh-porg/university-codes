clc, clear
%% Known/Chosen Values
ToD = input("What is the time of day in military time? (i.e. 6:30pm is 18:30, formatted as 18.50) ")
if ToD<6 || ToD>18
    error("Time of day is outside of time range")
else
    theta = ToD2theta(ToD);
end
Resolution = input("What is desired resolution? (m) ")
lambda = 0.0000005; % meters
wd = 3000 * 0.00006; % meters
Respon=10^(10); % V/W
MinDetec = 300; % mV/mu-meters
r = 0.3; % Reflectance of the ground
E_sun =1370; % W/m^2
tau_optics = 0.96;
gamma = 1.4;
R = 287.05; % Gas Constant of Air in J/kg-K
%% Sensor Calculation
if Resolution <= 3.5
    da = input("What is diameter of aperture? (m) ")
    if da > 0.1
        error("Aperture diameter is too wide")
    elseif da<0
        error("Aperture diameter is too narrow")
    else
    end
    A_aperture = 0.25 * pi * (da^2); % compute appature area
    H_da = res(Resolution, lambda); % Calculates the ratio of Altitude and Aperture Diameter given a desired resolution and chosen wavelength
    H = H_da*da % meters
    
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
    SW = Swath(H, wd, f) % meters
    if SW < 9656.064 % 6 mi in meters
        error("Swath Width too narrow")
    end
    theta_Rayleigh=1.22*(lambda/da)
    A_detected = 0.00006^2; % m^2
    Lr = ReflecRad(E_sun, r, theta)
    E_sensor = falloff(Lr, A_aperture,f, theta_Rayleigh, tau_optics)
    S = signal(Respon,E_sensor, A_detected)    % V/micrometer
    if S < 0.3
        error("Signal too low for detection")
    end
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

function E_sensor = falloff(Lr, A_aperture, f, theta, tau)    % Calculates the Irradience entering the sensor from the Refected radiance, the area of the aperture, and the radius of the detector
E_sensor=((Lr*A_aperture*((cosd(theta))^4))/(f^2))*tau;   % For right now it is being assumed that r_0 from the lense falloff eq is = to focal length f
end

function S = signal(Respon, E_sensor, A_detected)    % Calculates signal from the response time, sensor irradience, detected area, and time of day in degrees theta
S = Respon* E_sensor*A_detected;
end

function theta = ToD2theta(ToD) % Converts Time of Day to degrees theta
DegperHour=180/12; % Conversion between hour and degrees
theta=ToD*DegperHour-180; % Converts the Time of Day (ToD) to degrees with 0 degrees being Noon (12:00)
end