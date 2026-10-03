classdef Relation
    %RELATION class to store an equation and what vairalbes it constains. will symbolicly
    %solve for 1 unknown variable if all others are given.
    %   Version 2
    
    properties
        equationName % name of the equation that is being represented (string)
        relationship % equation (symbolic)
        variables % variables in the equation in the order defined by matlab (string array)
        varAssumptions % dictionary of assumptions about variables to be used in solving key (stirngs) -> value (string array)
        equationAssumptions % assumtion under which equation is valid (string array)
        equationTags % tags by which the equation can be found while looking for it in a list of equations (string array)
    end
    
    
    methods
        function self = Relation(relationshipString,options)
            arguments
                relationshipString string
                options.variableAssumptions dictionary = dictionary(); % default value is an empty dictionary
                options.equationName string = relationshipString % default value for equation name is just the equation is represents
                options.equationAssumptions (1,:) string = [""] % defaults to the set contianing only the empty string
                options.equationTags string = [""] % defaults to the set contianing only the empty string
            end

            self.relationship = str2sym(relationshipString);
            self.variables = string(symvar(self.relationship));

            self.varAssumptions = options.variableAssumptions;

            self.equationName = options.equationName;
            self.equationAssumptions = options.equationAssumptions;
            self.equationTags = options.equationTags;

            % check that only valid keys are in dictionary
            if self.varAssumptions.isConfigured
                keys = self.varAssumptions.keys;
                for i = 1:length(keys)
                    key = keys(i);
                    if (ismember(key,self.variables) || strcmp("equation")) == false %this does not corespond to a variable and is not an assumption equation
                        ME = MException("key %s not a variable or 'equation'",key); %raise an expection
                        throw(ME)
                    end
                end
            end

        end
        
        function out = solve(self,var_Matrix) % TODO: use optional arguments to add assumptions
            %Relation_Solver Solves the relation for the unknown variable
            %   Detailed explanation goes here

            relation = self.relationship;
            sym_Matrix = symvar(relation); % create a matrix of symbolic variable used in the relation
            
            if self.varAssumptions.isConfigured % do we have any mathematical assumptions about variables? (could use num entries =0 but that seems a bit worse)
                %make assumptions
                keys = self.varAssumptions.keys;
                for i = 1:length(keys)
                    key = keys(i);
                    assumps = self.varAssumptions(key); % get the assumptions associated with the key
                    assumps = assumps{1}; % unwrap it from cell array to string array

                    if strcmp(key,"equation") % if the assumption is in the form of an equation
                        for assumption = assumps()
                            assumeAlso(assumption) % then directly assume the equations
                        end
                    else % thus the assumption applies directly to a variable
                        assumptionVar = str2sym(key); % get the variable it is refering to
                        for assumption = assumps()
                            assumeAlso(assumptionVar,assumption) % and apply the relevant assumptions
                        end
                    end
                end
            end


            %assume(sym_Matrix,'real')

            
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

