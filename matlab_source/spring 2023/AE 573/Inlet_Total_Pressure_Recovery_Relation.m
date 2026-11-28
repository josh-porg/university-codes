classdef Inlet_Total_Pressure_Recovery_Relation < Relation
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        relationship = str2sym('x2/x3 - x1');
        variables = ["pi_d", "pt_0", "pt_2"];
        variableNames = ["Inlet Total Rressure Recovery", "Free-Stream Total Presure", "Engine Face Total Pressure"];
    end
end

