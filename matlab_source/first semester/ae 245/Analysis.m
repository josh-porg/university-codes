classdef Analysis
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
%     methods
%         function obj = untitled()
%             
%         end
%         
%         function outputArg = method1(obj,inputArg)
%            
%         end
%     end
    
    methods (Static)
        function error = calculateError(demanded, actual)
            error = (demanded - actual) / demanded;
        end
    end
end

