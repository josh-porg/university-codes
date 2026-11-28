classdef Material
    %Material basic characteristics of a materiale
    %   Detailed explanation goes here

    properties
        Name
        E
        F_TU
        F_TY
        rho
        cost_density
    end

    methods
        function obj = Material(Name, E,F_TU,F_TY,rho,cost_density)
            %Material Construct an instance of this class
            %   Detailed explanation goes here
            obj.Name = Name;
            obj.E = E;
            obj.F_TU = F_TU;
            obj.F_TY = F_TY;
            obj.rho = rho;
            obj.cost_density = cost_density;
        end
    end
end