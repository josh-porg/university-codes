classdef RocketBullet
    %RocketBullet an object to contian rocket bullet data
    %   Detailed explanation goes here

    properties
        Name
        Hydrolic_Diameter
        Nozzle_Exit_Diameter
        Projectile_Length
        Nose_Length
        Nose_Tip_Diameter
        Boattail_Diameter
        Nozzle_Exit_Area
        Projectile_Mass % wet mass
        Impulse_Specific
        Propellant_Mass % mass of propellant
        Motor_On_Time
        Motor_Off_Time
    end

    methods
        function obj = RocketBullet(Name, Hydrolic_Diameter,Nozzle_Exit_Diameter,Projectile_Length,Nose_Length,Nose_Tip_Diameter,Boattail_Diameter,Projectile_Mass,Impulse_Specific,Propellant_Mass,Motor_On_Time,Motor_Off_Time)
            %RocketBullet Construct an instance of this class
            %   Detailed explanation goes here
            obj.Name = Name;
            obj.Nozzle_Exit_Area = pi * (Nozzle_Exit_Diameter/2)^2;
            obj.Hydrolic_Diameter = Hydrolic_Diameter;
            obj.Nozzle_Exit_Diameter = Nozzle_Exit_Diameter;
            obj.Projectile_Length = Projectile_Length;
            obj.Nose_Length = Nose_Length;
            obj.Nose_Tip_Diameter = Nose_Tip_Diameter;
            obj.Boattail_Diameter = Boattail_Diameter;
            obj.Projectile_Mass = Projectile_Mass;
            obj.Impulse_Specific = Impulse_Specific;
            obj.Propellant_Mass = Propellant_Mass;
            obj.Motor_On_Time = Motor_On_Time;
            obj.Motor_Off_Time = Motor_Off_Time;
        end
    end
end