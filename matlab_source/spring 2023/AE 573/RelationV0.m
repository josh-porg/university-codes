classdef (Abstract) RelationV0
    %Relation Defines a relationship between variable and provides a
    % method for solving for them
    %   Detailed explanation goes here
    
    properties (Abstract, Constant)
        relationship
        variables
        variableNames
    end
    
    methods (Static)
%         function obj = untitled(inputArg1,inputArg2)
%             %UNTITLED Construct an instance of this class
%             %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
%         end
        
        function out = solve(var_Matrix)
            %Relation_Solver Solves the relation for the unknown variable
            %   Detailed explanation goes here
            %implementingClass = meta.class.fromName(mfilename('Inlet_Total_Pressure_Recovery_Relation'));
            %foo = implementingClass.PropertyList(1).ConstantValue
            relation = relationship;
            sym_Matrix = symvar(relation); % create a matrix of symbolic variable used in the relation
            
            for i = [1:length(var_Matrix)] % itterate through the given variables
                if(isnan(var_Matrix(i)) == false) % if it is known
                    % substitute the known value in place of its respecitve symbolic variable
                    relation = subs(relation,sym_Matrix(i),var_Matrix(i)); 
                else %it is out unknown variable
                    %out_index = i;
                    solve_var = sym_Matrix(i); % mark this variable to be returned
                end
            end
            
            out = solve(relation,solve_var); % return the newly solved for variable
        end
    end
end

