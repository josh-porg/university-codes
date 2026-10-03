classdef Ply
    %Ply Summary of this class goes here
    %   Othotropic ply: note function desciption refrences tensile
    %   properties. Ply does not distinguish between tinsile and
    %   compressive properties

    properties
        % Material Properies
        % tensile stiffness properties
        E_1 % Tensile Youngs modulus in principal direction (Pa)
        E_2 % Tensile Youngs modulus in transverse direction (Pa)

        % shear stiffness properties
        G_12 % Shear modulus in plane (Pa)
        nu_12 % possion's ratio in plane (Pa)

        % Stress allowables
        stressAllowables % [F_1_t, F_2_t, F_1_c, F_2_c, tau_12]
        % where
        % F_1_t is Tensile failure stress in principal direction (Pa)
        % F_2_t is Tensile failure stress in transverse direction (Pa)
        % F_1_c is Compressive failure stress in principal direction (Pa)
        % F_2_c is Compressive failure stress in transverse direction (Pa

        % Stress allowables
        strainAllowables % [epsilon_1_t, epsilon_2_t, epsilon_1_c, epsilon_2_c, gamma_12]
        % where
        % epsilon_1_t is Tensile failure strain in principal direction (Pa)
        % epsilon_2_t is Tensile failure strain in transverse direction (Pa)
        % epsilon_1_c is Compressive failure strain in principal direction (Pa)
        % epsilon_2_c is Compressive failure strain in transverse direction (Pa)
        

        Q % Reduced Stiffness Matrix

        % invariants ( used for rotations )
        U_1
        U_2
        U_3
        U_4
        U_5

        failureStressesDefined
        failureStrainsDefined
        combinedFailureDefined

    end

    properties(Dependent)
        % Stress allowable properties
        F_1_t % Tensile failure stress in principal direction (Pa)
        F_2_t % Tensile failure stress in transverse direction (Pa)

        F_1_c % Compressive failure stress in principal direction (Pa)
        F_2_c % Compressive failure stress in transverse direction (Pa)


        % Strain allowable properties
        epsilon_1_t % Tensile failure strain in principal direction (Pa)
        epsilon_2_t % Tensile failure strain in transverse direction (Pa)

        epsilon_1_c % Compressive failure strain in principal direction (Pa)
        epsilon_2_c % Compressive failure strain in transverse direction (Pa)
    end

    methods
        function this = Ply(E_1, E_2, G_12, nu_12, varargin)
            %Ply Construct an instance of this class
            %   Creates a Ply for use in CLT. stress alowwables not
            %   currently implemented if stress allowables unknown put in
            %   an arbitrary number

            % Initialize failure properties
            this.stressAllowables = zeros([1,5]);
            this.strainAllowables = zeros([1,5]);
            this.failureStressesDefined = false;
            this.failureStrainsDefined = false;
            this. combinedFailureDefined = false; % currently this will always be false because I is not yet implemented

            % Parse optional arguments
            if nargin > 4
                for i = 1:2:length(varargin)
                    switch varargin{i}
                        case 'StressAllowables'
                            this.stressAllowables = varargin{i+1};
                            this.failureStressesDefined = true;
                        case 'StrainAllowables'
                            this.strainAllowables = varargin{i+1};
                            this.failureStrainsDefined = true;
                        otherwise
                            error('Ply:InvalidArgument', 'Unknown parameter name: %s', varargin{i});
                    end
                end
            end

            % mechanical properties
            this.E_1 = E_1;
            this.E_2 = E_2;
            this.G_12 = G_12;
            this.nu_12 = nu_12;

            % compute poisson's ratio in 21 direction 
            nu_21 = nu_12 .* E_2 ./ E_1;
            % compute reduced stiffnesses of ply
            Q_11 = E_1 ./ (1 - nu_12 .* nu_21);
            Q_22 = E_2 ./ (1 - nu_12 .* nu_21);
            Q_12 = nu_12 .* E_2 ./ (1 - nu_12 .* nu_21);
            Q_66 = G_12;

            % compute reduced stifness matrix
            this.Q = [Q_11, Q_12, 0;...
                 Q_12, Q_22, 0;...
                 0   , 0   , Q_66];

            % compute invariants ( for rotation into global space )
            this.U_1 = (3*Q_11 + 3*Q_22 + 2*Q_12 + 4*Q_66) / 8;
            this.U_2 = (Q_11 - Q_22) / 2;
            this.U_3 = (Q_11 + Q_22 - 2*Q_12 - 4*Q_66) / 8;
            this.U_4 = (Q_11 + Q_22 + 6*Q_12 - 4*Q_66) / 8;
            this.U_5 = (Q_11 + Q_22 - 2*Q_12 + 4*Q_66) / 8;

        end

        function Q_bar = Qbar(this,theta)
            %Qbar Returns the stiffness matrix in gloabl space of a ply
            %at angle theta
            %   Takes: angle of ply, theta (radians)
            %   Returns: stiffness matrix in gloabl space, Q_bar ()
            
            % compute stiffnesses
            Q_bar_11 = this.U_1 + this.U_2 .* cos(2*theta) + this.U_3 .* cos(4*theta);
            Q_bar_12 = this.U_4 - this.U_3 .* cos(4*theta);
            Q_bar_22 = this.U_1 - this.U_2 .* cos(2*theta) + this.U_3 .* cos(4*theta);
            Q_bar_16 = 1/2 * this.U_2 .* sin(2*theta) + this.U_3 .* sin(4*theta);
            Q_bar_26 = 1/2 * this.U_2 .* sin(2*theta) - this.U_3 .* sin(4*theta);
            Q_bar_66 = this.U_5 - this.U_3 .* cos(4*theta);
            
            % assemble stiffness matrix
            Q_bar = [Q_bar_11, Q_bar_12, Q_bar_16;...
                     Q_bar_12, Q_bar_22, Q_bar_26;...
                     Q_bar_16, Q_bar_26, Q_bar_66];
        end

        function [hasFailed, Failures, FailureRatios] = TestFailure(this, theta, epsilon, RunningLoads, failureCriteria, plyNumber)
            % TODO: allow it to passed agnositc to the failure criteria and
            % let the mehtod work it out itself

            % rotate the strain into plyCoordinates
            % trasnformation matrix from laminate strains to lamina strains
            T_2 = [
                cos(theta).^2, sin(theta).^2, -sin(theta).*cos(theta);
                sin(theta).^2, cos(theta).^2, sin(theta).*cos(theta);
                2*sin(theta).*cos(theta), -2*sin(theta).*cos(theta), cos(theta).^2 - sin(theta).^2
                ];

            % rotate the laminate strains into principal coordinates
            epsilon_pricipal = T_2 \ epsilon; % ply strians in principal coordinates

            % Input Validation
            if (failureCriteria == "Max_Strain" && this.failureStrainsDefined == false)... % strains used but allowable not defined
                    || (failureCriteria == "Max_Stress" && this.failureStressesDefined == false) % stresses used but allowable not defined
                error("Ply:UndefinedFailureParameters", "Failure Cireteria for %s failure theory not defined", failureCriteria); % failure criteria have not been defined
            elseif(failureCriteria ~= "Max_Strain" && failureCriteria ~= "Max_Stress")
                error("Ply:NotImplemented", "%s Failure theory not implemented", failureCriteria); % Explain the other failure theory are not implemented
            end

            switch(failureCriteria)
                case "Max_Strain"
                    if(this.failureStrainsDefined == false)
                        error("Ply:UndefinedFailureParameters", "Failure Cireteria for %s failure theory not defined", failureCriteria);
                    end
                    actual = epsilon_pricipal; % use the strain as the actual
                    % use strain allowables
                    allowable_tensile = this.strainAllowables(1:2);
                    allowable_compressive = this.strainAllowables(3:4);
                    allowable_shear = this.strainAllowables(5);
                case "Max_Stress"
                    if(this.failureStressesDefined == false)
                        error("Ply:UndefinedFailureParameters", "Failure Cireteria for %s failure theory not defined", failureCriteria);
                    end

                    % compute stress from strains
                    sigma_pricipal = this.Q * epsilon_pricipal; % compute principal stresses
                    actual = sigma_pricipal; % use the stress as the actual
                    % use stress alowables
                    allowable_tensile = this.stressAllowables(1:2);
                    allowable_compressive = this.stressAllowables(3:4);
                    allowable_shear = this.stressAllowables(5);
                case "Combined"
                    error("Ply:NotImplemented", "%s Failure theory not implemented", failureCriteria);
                otherwise
                    error("Ply:NotImplemented", "%s Failure theory not implemented", failureCriteria);
            end
            

            % acept omisssion of ply number
            %             switch nargin
            %                 case 4 % no failure criteria
            %                     F_1_t = 0;
            %                     F_2_t = 0;
            %                     F_1_c = 0;
            %                     F_2_c = 0;
            %                     warning("Ply:Defaults","Ply failure stress values not provided");
            %                     this.failureStressesDefined = false;
            %                 case 8 % failure stresses fully defined
            %                     this.failureStressesDefined = true;
            %             end

            %Initialize Failures arrays
            Failures = []; % array to store all the failures

            %strain ratios (allowable/actual)
            R_t(1:2) = allowable_tensile ./ actual(1:2);
            R_c(1:2) = allowable_compressive ./ actual(1:2);

            % for shear strain (which doesn't depend on sign) take abs
            %failureRatios(epsilon_t, epsilon_c, gamma_12, epsilon_pricipal)
            R_tot = [R_t'; R_c'; abs(allowable_shear / actual(3))]; % full list of failure modes
            hasFailed = any(R_tot < 1 & R_tot > 0); % strain ratio < 1 is failure strain ratio < 0 is wrong direction


            % what to do if we have failed
            if hasFailed
                % describe failure
                FailureDir = R_tot < 1 & R_tot > 0;

                %failureModes
                failureCode = find(FailureDir);

                % produce failure objects
                for failuretype = failureCode'
                    % compute failure loads
                    failureLoad = R_tot(failuretype) * RunningLoads;
                    % store failure in a failure object
                    Failures = [Failures, PlyFailure(plyNumber, failuretype, failureCriteria, failureLoad)];
                end
            end

            FailureRatios = R_tot;
        end

    end

    methods (Static)
        function [E_1, E_2, nu_12, G_12] = mechanicalProperties(E_f, G_f, nu_f, E_m, G_m, nu_m, V_f, xi_transverse, xi_shear)
            %mechanicalProperties Determines the mechanical properties of a
            %ply composed of matrix and fiber with known characteristics
            %and a known volume fraction.
            %
            %   TAKES: Young's modulus of the fiber, E_f (Pa); Shear
            %   modulus of the fiber, G_f (msi); Possion's ratio of the
            %   fiber, nu_f (~); Young's modulus of the matrix, E_m (Pa);
            %   Shear modulus of the matrix, G_m (msi); Possion's ratio of
            %   the matrix, nu_m (~); Fiber volume fraction, V_f (~).
            %
            %   OPTIONAL ARGUMENTS: Halpin-Tsai coefficient for transverse
            %   modulus, xi_transverse (~); Halpin-Tsai coefficient for
            %   shear modulus, xi_shear (~).
            %
            %   RETURNS: ply axial Young's modulus, E_1 (Pa); ply
            %   transverse Young's modulus, E_2 (Pa); ply Poissions ratio,
            %   nu_12 (~); ply shear modulus, G_12 (Pa).
            %
            %   ASSUMES: matrix volume fraction V_m + fiberber volume frication V_f = 1.
            %   true for plys consisting of 1 fiber and 1 matrix material with negligible
            %   interphase.
            %
            %   NOTE: For circular or square cross section xi_transverse =
            %   2. For rectangualr cross section xi_transverse = 2 * a/b.
            %
            %   NOTE: if only one value of xi is provided it is assumed to
            %   be xi_transverse
            switch(nargin)
                case 7 % did the user ommit xi?
                    xi_shear = 1; % value is approximately 1 in most cases
                    xi_transverse = 2; % assume cicular or square cross section
                    warning("Ply:Defaults","xi values not provided: using defaults")
                case 8 % did the user ommit xi?
                    xi_shear = 1; % value is approximately 1 in most cases
                    warning("Ply:Defaults","xi_shear value not provided: using default")
                case 9 % all values provided we are happy
                    % do notheing
                otherwise
                    error("invalid number of input args: %g", nargin);
            end

            % now calculate the mechanical properties

            % mechanics of matrials approach for axial modulus and poisson's ratio
            E_1 = E_f.*V_f + E_m.*(1-V_f); % axial modulus of ply (Pa)
            nu_12 = nu_f.*V_f + nu_m.*(1-V_f); % poisson's ratio of ply (~)

            % Halpin-Tsai equations for transverse and shear moduli
            E_2 = Ply.halpinTsai(E_f, E_m, V_f, xi_transverse); % Transverse Young's modulus of ply (Pa)
            G_12 = Ply.halpinTsai(G_f, G_m, V_f, xi_shear); % Shear modulus of ply (Pa)
        end

    end

    methods (Static, Access = private)
        function modulus = halpinTsai(modulus_f, modulus_m, V_f, xi)
            %halpinTsai Used to determine the transverse or shear modulus
            %of a ply
            % TAKES: fiber modulus, modulus_f (Pa); matrix modulus,
            % modulus_m (Pa); Fiber volume fraction V_f, (~); Halpin-Tsai
            % coeffcient, xi (~).
            % RETURNS: Ply transverse or shear modulus, modulus (Pa).
            % NOTE: Only use for transevers or shear modulus!
            eta = (modulus_f./modulus_m - 1) ./ (modulus_f./modulus_m + xi);
            modulus = modulus_m .* (1 + xi.*eta.*V_f) ./ (1 - eta.*V_f);
        end
    end


end