classdef GasLaw
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        equationOfState
        cpEquation
        formulas
    end

    methods
        function self = GasLaw(equationOfState,cpEquation)
            %Creates equations from a gas law and maxwell relations
            %   Detailed explanation goes here

            % store inputs
            self.equationOfState = equationOfState;
            self.cpEquation = cpEquation;

            % create symbolic versions of inputs
            EoS = str2sym(equationOfState);

            %% begin to generate partial derivatives
            syms p v T p_f(v,T) v_f(p,T) T_f(p,v)

            % partials of pressure
            %press with respect to temp holding vol const
            dglpt = diff(subs(EoS,p,p_f),T);
            dglptv = rhs(isolate(dglpt,diff(p_f,T)));

            % press with respect to vol holding temp constant
            dglpt = diff(subs(EoS,p,p_f),v);
            dglpvt = rhs(isolate(dglpt,diff(p_f,v)));

            %remove function notation
            formulas.Partial.PrTcV = subs(dglptv,p_f,p);
            formulas.Partial.PrVcT = subs(dglpvt,p_f,p);

            %now derivatives of tempreature
            % temp with respect to press holding vol constant
            dgltp = diff(subs(EoS,T,T_f),p);
            dgltpv = rhs(isolate(dgltp,diff(T_f,p)));

            % temp with respect to vol holding press constant
            dgltv = diff(subs(EoS,T,T_f),v);
            dgltvp = rhs(isolate(dgltv,diff(T_f,v)));

            %remove function notation
            formulas.Partial.TrPcV = subs(dgltpv,T_f,T);
            formulas.Partial.TrVcP = subs(dgltvp,T_f,T);

            %partial derivatives of volume
            % vol with respect to temp holding press constant
            dglvt = diff(subs(EoS,v,v_f),T);
            dglvtp = rhs(isolate(dglvt,diff(v_f,T)));

            % vol with respect to press holding temp constant
            dglvp = diff(subs(EoS,v,v_f),p);
            dglvpt = rhs(isolate(dglvp,diff(v_f,p)));

            %remove function notation
            formulas.Partial.VrTcP = subs(dglvtp,v_f,v);
            formulas.Partial.VrPcT = subs(dglvpt,v_f,v);

            %check for satisfaction of reciorocity relation (TODO: check cyclic relation as well)
            if isequal(1/formulas.Partial.TrVcP,formulas.Partial.VrTcP) == false % if does not satisfy reciprocity relation
                disp("does not satisfy resiprocity relation for volume and temperature at constant prssure")
            end
            if isequal(1/formulas.Partial.PrVcT,formulas.Partial.VrPcT) == false % if does not satisfy reciprocity relation
                disp("does not satisfy resiprocity relation for pressure and volume at constant temperature")
            end
            if isequal(1/formulas.Partial.TrPcV,formulas.Partial.PrTcV) == false % if does not satisfy reciprocity relation
                disp("does not satisfy resiprocity relation for pressure and temperature at constant volume")
            end

            %% specific heat
            formulas.cp = str2sym(cpEquation);

            %% maxwell relations and other derived relations
            %maxwell relations
            formulas.Partial.SrVcT = formulas.Partial.PrTcV;
            formulas.Partial.SrPcT = -formulas.Partial.VrTcP;

            % temp and specific volume are independent
            formulas.Partial.UrVcT = T * formulas.Partial.PrTcV - p; % this is internal pressure
            formulas.Partial.UrPcT = formulas.Partial.UrVcT * formulas.Partial.VrPcT; %from the above and basic math

            % temp and pressure are independent
            formulas.Partial.HrPcT = v - T * formulas.Partial.VrTcP;
            formulas.Partial.HrVcT = T * formulas.Partial.PrTcV + v * formulas.Partial.PrVcT; %from the above and basic math

            %from specific heats relationship
            formulas.cv = formulas.cp + T*(formulas.Partial.VrTcP^2 * formulas.Partial.PrVcT);

            formulas.Partial.HrTcP = formulas.cp; % by definition
            formulas.Partial.UrTcV = formulas.cv; % by definition
            formulas.Partial.SrTcP = formulas.cp / T; % temp and press are independent
            formulas.Partial.SrTcV = formulas.cv / T; % temp and vol are independent
            formulas.Partial.TrVcu = 1/formulas.cv * (p - T*formulas.Partial.PrTcV); % joule coefficient
            formulas.Partial.TrPch = 1/formulas.cp * (T*formulas.Partial.VrTcP - v); % joule-thomson coefficient
            formulas.Partial.PrVcS = formulas.cp/formulas.cv * formulas.Partial.PrVcT; % from maxwell funny buisness on ratio of specific heats

            %% Deltas
            %entropy
            formulas.Delta.SrTV = int(formulas.Partial.SrTcV,T) + int(formulas.Partial.SrVcT,v);
            formulas.Delta.SrTP = int(formulas.Partial.SrTcP,T) + int(formulas.Partial.SrPcT,p);
            %enthalpy
            formulas.Delta.HrTP = int(formulas.Partial.HrTcP,T) + int(formulas.Partial.HrPcT,p);
            %internal energy
            formulas.Delta.UrTV = int(formulas.Partial.UrTcV,T) + int(formulas.Partial.UrVcT,v);

            %TODO: experiment with using definite integrals

            self.formulas = formulas; % assign to object
        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end