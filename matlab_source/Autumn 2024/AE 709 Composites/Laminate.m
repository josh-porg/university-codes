classdef Laminate
    %Ply Summary of this class goes here
    %   Laminate properties determined from ply properties and stacking
    %   sequence.
    %   ASSUMES: classical lamination theory for thin laminates with plane
    %   stress and plane strain, where a line originally straight and
    %   perpendeicular to the neutral axis remains perpendicular to the
    %   neutral axis when the laminate is extended and bent. the line is
    %   also assumed to have a constant length.

    properties
        % Material Properies
        % tensile stiffness properties
        E_x % Tensile Youngs modulus in principal direction (Pa)
        E_y % Tensile Youngs modulus in transverse direction (Pa)

        % shear stiffness properties
        G_xy % Shear modulus in plane (Pa)
        nu_xy % possion's ratio in plane (Pa)

        % CLT Matricies - consider making these dependent properties and
        % using getters and setters
        %A % A Stiffness Matrix
        %B % B Stiffness Matrix
        %D % D Stiffness Matrix
        % this contains the ABD matricies
        laminateStiffness
        laminateCompliance % inverse of the stiffness matrix


        % Layup Properties

        Plies % an array of all the Ply objects that will be used in this laminate
        
        % Ply Material Matrix assigns the Ply material to each layer. array of n-elements where the
        % value of the ith element corresponds to the index of the material in the
        % Plies matrix that should be applied to the ith ply of the laminate
        Ply_Material_Matrix % what Ply to assign to each layer (Ply)

        n_Plys % number of plys in the laminate (~)

        symmetry % Laminate symmetery: 0 -> asymmetric; 1 -> odd symmetric; 2 -> even symmetric; -1 odd anti-symmetric -2; even anti-symmetric (~)

        t_ply % array of the thicknesses of each ply in order of topmost ply to bottom most ply
        t_laminate; % total thickness of the laminate

        plyAngles % array of ply angles (Radians)

        z % an array of the locations of each ply surface where neutral axis is at zero, the upper surface of the top-most ply is at largest negative number and bottom surface of the bottommost ply is at the most positive z value
    end
    properties(Dependent)
        A
        B
        D
    end

    methods
        function this = Laminate(Plies, Ply_Material_Matrix, t_ply, plyAngles, symmetry)
            % Ply Construct an instance of this class
            %    Creates A Laminate Using CLT.
            % TAKES: Plies, row vector of Ply objects containin the set of
            % plys that will be used (Ply); Ply_Material_Matrix, a vector
            % in which an element at index i with value j assigns the ith
            % layer from the top of laminate to use the material in index j
            % of the Plies array, if only material is used an empty array
            % can be passed (~); t_ply thickness of each ply either an
            % array containing the thicknes of each layer or a scalar for 
            % uniform layer thicknesses (m); plyAngles, a vector containing
            % the angle of each ply, ordered from top-most ply to
            % bottom-most ply (radians); of laminate, symmetery: 0 ->
            % asymmetric; 1 -> odd symmetric; 2 -> even symmetric; -1 odd
            % anti-symmetric -2; even anti-symmetric (odd symmetry must not
            % have an anlge ply as the mirror ply).
            %
            % ASSUMES: classical lamination theory for thin laminates with plane
            %   stress and plane strain, where a line originally straight and
            %   perpendeicular to the neutral axis remains perpendicular to the
            %   neutral axis when the laminate is extended and bent. the line is
            %   also assumed to have a constant length.

            % warn the user about the assumptions in Clasical Lamination
            % Theory.
            warning("Laminate:CLT", "ASSUMES: classical lamination theory for thin laminates with plane stress and plane strain," + ...
                " where a line originally straight and perpendeicular to the neutral axis" + ...
                " remains straight and perpendicular to the neutral axis when the laminate is extended and bent." + ...
                " The line is also assumed to have a constant length.");

            % input validation

            % it has been observed that a common user error is to enter an
            % array in degrees rather than radians we will check if the
            % user has entered any common angle mistakes and warn them
            if(sum(abs(plyAngles) > pi/2) > 0) % these angles any angle greater than pi radians is most likey an (improperly converted) angle in degrees.
                 error("Laminate:degrees", "Input Laminate angles %s exeeds than +-pi radians. keep Laminate anlges between zero and pi radians (inclusive)", mat2str(plyAngles(plyAngles > pi/2)));
            end

            % check for invalid odd symetry
            if(abs(symmetry) == 1 && (plyAngles(end) ~= 0 && plyAngles(end) ~= deg2rad(90)) ) % does the an odd symmetric laminate have an angle ply on the symetry plane?
                % remind the user that it is impossilbe to have odd
                % symmetry with an angle ply as the mirror ply
                % NOTE: these are all in radians (this may supprise user)
                if(symmetry == 1) % if it's odd symetric its just weird but the answer will be the same. (failure criteria may be different requiring a split of the ply)
                    warning("Possibly hazardous outcome: Ply Anlges <%s> are risky because of odd symmetry about an angle ply (%g)", mat2str(plyAngles), plyAngles(end))
                else % if it's odd antisymetric it is impossilbe
                    error("Ply Anlges <%s> are invalid because you cannot have odd antisymmetry about an angle ply (%g)", mat2str(plyAngles), plyAngles(end));
                end
            end

            % ensure that the arays passed are row vectors
            reshape(Ply_Material_Matrix, 1, []); % reshape into a row vector
            reshape(plyAngles, 1, []); % reshape into a row vector
            reshape(t_ply, 1, []); % reshape into a row vector


            
            % stroe component Ply objects
            this.Plies = Plies;

            % assign matrials to each ply if not already configured
            if(isempty(Ply_Material_Matrix)) % is an empty matrix
                Ply_Material_Matrix = ones(1,(length(plyAngles))); % assume only one ply material used
            elseif length(Ply_Material_Matrix) ~= length(plyAngles) % check for invalid matrix length
                error("length of ply materialx matix %g, unequal to length of angle matrix %g", length(Ply_Material_Matrix), length(plyAngles)); % warn of the issue
            end % if we have made it thus far then the angle matrix must be valid

            % assign thicknesses if not already computed
            if isscalar(t_ply) % is only one thickness is specified
                t_ply = ones(1,length(plyAngles)) * t_ply; % assume all plies are the same thickness
            elseif length(t_ply) ~= length(plyAngles) % each ply has its own thicknesses
                error("length of thickness matix %g, unequal to length of angle matrix %g", length(t_ply), length(plyAngles)); % warn of the issue
            end

            % compute symetries
            this.symmetry = symmetry;

            % apply the correct adjustments to relevant matricies for symmetry
            switch symmetry
                case 0 % laminate is asymmetric
                    % assign the given input values
                    this.plyAngles = plyAngles;
                    this.t_ply = t_ply;
                    this.Ply_Material_Matrix = Ply_Material_Matrix;
                case 1 % laminate has odd symmetry
                    % add a mirror to the the laminate about the center of
                    % the final ply
                    this.plyAngles = horzcat(plyAngles, fliplr(plyAngles(1:end-1)));
                    this.t_ply = horzcat(t_ply, fliplr(t_ply(1:end-1)));
                    this.Ply_Material_Matrix = horzcat(Ply_Material_Matrix, fliplr(Ply_Material_Matrix(1:end-1)));
                case 2 % laminate has even symmetry
                    % mirror the plies duplicating the central element
                    this.plyAngles = horzcat(plyAngles, fliplr(plyAngles(1:end)));
                    this.t_ply = horzcat(t_ply, fliplr(t_ply(1:end)));
                    this.Ply_Material_Matrix = horzcat(Ply_Material_Matrix, fliplr(Ply_Material_Matrix(1:end)));
                case -1 % laminate has odd anti-Symmetry
                    % compute the anti-symmetric angle compliment matrix
                    antiSymmetricAngles = -fliplr(plyAngles(1:end-1)); % flip the order and switch the sign of all plies
                    antiSymmetricAngles(antiSymmetricAngles == -deg2rad(90)) = deg2rad(90); % -90 degree plies are equivalent to 90 degree plies
                    % mirror the plies maintianing a single central element
                    this.plyAngles = horzcat(plyAngles, antiSymmetricAngles);
                    this.t_ply = horzcat(t_ply, fliplr(t_ply(1:end-1)));
                    this.Ply_Material_Matrix = horzcat(Ply_Material_Matrix, fliplr(Ply_Material_Matrix(1:end-1)));
                case -2 % laminate has even anti-Symmetry
                    % compute the anti-symmetric angle compliment matrix
                    antiSymmetricAngles = -fliplr(plyAngles(1:end)); % flip the order and switch the sign of all plies
                    antiSymmetricAngles(antiSymmetricAngles == -deg2rad(90)) = deg2rad(90); % -90 degree plies are equivalent to 90 degree plies
                    % mirror the plies maintianing a single central element
                    this.plyAngles = horzcat(plyAngles, antiSymmetricAngles);
                    this.t_ply = horzcat(t_ply, fliplr(t_ply(1:end)));
                    this.Ply_Material_Matrix = horzcat(Ply_Material_Matrix, fliplr(Ply_Material_Matrix(1:end)));
                otherwise
                    error("symmetry %g, not valid. symmetry must be 0 -> asymmetric; 1 -> odd symmetric; 2 -> even symmetric; -1 -> odd anti-symmetric; -2 -> even anti-symmetric");
            end



            this.n_Plys = length(this.plyAngles);

            % for readablitlity
            t_ply = this.t_ply;
            Ply_Material_Matrix = this.Ply_Material_Matrix;
            plyAngles = this.plyAngles;

            % compute the CLT matricies (A,B,D)

            % compute matrix of ply boundary locations
            if isscalar(t_ply) % is only one thickness is specified
                z = [0:this.n_Plys] * t_ply; % assume all plies are the same thickness
            elseif length(t_ply) == this.n_Plys % each ply has its own thicknesses
                z = zeros(1, this.n_Plys+1); % initialize ply boundary locations
                for i = 2:this.n_Plys+1
                    z(i) = z(i-1) + t_ply(i-1); % the location of each subsequent boundary is the location of the previous boundary + the ply thickness
                end
            else % check for invalid matrix size
                error("length of thickness matix %g, unequal to length of angle matrix %g", length(t_ply), length(plyAngles)); % warn of the issue
            end

            % compoute the laminate thickness
            this.t_laminate = max(z);

            % account for the nuertal axis
            z_NA = max(z)/2; % assume that the neural axis is at half the thickness of the laminate
            z = z - z_NA; % adjust the coordiantes of the ply boundaries to be relative to the neutral axis
            this.z = z; % store the value in the object

            % compute CLT matricies

            % compute the A maxtrix
            A = zeros(3,3); % preallocate A matrix
            for k = 1:this.n_Plys
                A = A + Plies(Ply_Material_Matrix(k)).Qbar(plyAngles(k)) * (z(k+1) - z(k));
            end

            % compute the B matrix
            B = zeros(3,3); % preallocate B matrix
            for k = 1:this.n_Plys
                B = B + (1/2) * Plies(Ply_Material_Matrix(k)).Qbar(plyAngles(k)) * (z(k+1)^2 - z(k)^2);
            end

            % compute the D matrix
            D = zeros(3,3); % preallocate D matrix
            for k = 1:this.n_Plys
                D = D + (1/3) * Plies(Ply_Material_Matrix(k)).Qbar(plyAngles(k)) * (z(k+1)^3 - z(k)^3);
            end

            % use known simplifications to remove numerical error in special cases
            % somewhat risky but it makes the B matrix term not include 1*10^-13 terms
            % which may be important for inversion of the matrix

%             % remove numerical error for symmetric laminates
%             if(symmetry == 1 || symmetry == 2) % is the laminate symmetric
%                 B = zeros(3,3); % the B matrix is zero
%             end
%             
%             % remove numerical error for anti-symmetric laminates
%             if(symmetry == -1 || symmetry == -2) % is the laminate symmetric
%                 A(1,3) = 0; A(3,1) = 0; % simplify the A matrix
%                 A(2,3) = 0; A(3,2) = 0; % simplify the A matrix
%                 D(1,3) = 0; D(3,1) = 0; % simplify the D matrix
%                 D(2,3) = 0; D(3,2) = 0; % simplify the D matrix
%             end
%             
%             % remove numerical error for no angle ply laminates
%             hasAnglePlies = sum(abs((plyAngles ~= pi/2) .* plyAngles),"all") > 0; % determine if there are any angle plies
%             if(hasAnglePlies == false) % there are no angle plies?
%                 B(1,2) = 0; B(2,1) = 0; % simplify the B matrix
%                 B(1,3) = 0; B(3,1) = 0; % simplify the B matrix
%                 B(2,3) = 0; B(3,2) = 0; % simplify the B matrix
%                 B(3,3) = 0; % simplify the B matrix
%                 D(1,3) = 0; D(3,1) = 0; % simplify the D matrix
%                 D(2,3) = 0; D(3,2) = 0; % simplify the D matrix
%             end


            % compute laminate Bulk Response
            % assemble the total relation between running loads and strains
            this.laminateStiffness = [A,B;...
                B,D];
            this.laminateCompliance = inv(this.laminateStiffness);

            this.E_x = 1 / (this.t_laminate * this.laminateCompliance(1,1));
            this.E_y = 1 / (this.t_laminate * this.laminateCompliance(2,2));
            this.G_xy = 1 / (this.t_laminate * this.laminateCompliance(3,3));
            this.nu_xy = -this.laminateCompliance(1,2)/this.laminateCompliance(1,1);

        end

        function [failedPlies, fails, failRats] = Load(this, RunningLoads, failureCriteria)
            % create vecotors of failed plies
            failedPlies = zeros(size(this.plyAngles)); % list of whether a ply as failed
            fails = []; % list of all the fialures that have occured
            failRats = zeros([length(this.plyAngles),2,5]); % list of all the failure ratios

            % get Laminate neutral axis strains and curvatures due to applied load
            laminateStrainsAndCurvatures = this.laminateCompliance() * RunningLoads;

            epsilon_0 = laminateStrainsAndCurvatures(1:3); % laminate neutral axis strains
            kappa = laminateStrainsAndCurvatures(4:end); % laminate neitral axis curvatures

            % compute ply Strains
            for i = 1:2 % perfom once for upper surface and once for lower surface
                z_surf = this.z(i:end+i-2); % the distances of the upper (and lower) surfaces (depending on itteration

                % convert into ply strains at surface in laminate coordinates
                epsilon_LamCoords = cell2mat(arrayfun(@(z) epsilon_0 + z * kappa, z_surf, 'UniformOutput', false));

                % convert into ply strains in material coordinates
                for k = 1:length(this.plyAngles)
                    epsilon_Lam = epsilon_LamCoords(:,k); % laminate coordinate strains
                    theta = this.plyAngles(k);

                    % determine failure
                    [hasFailed, Failures, FailureRatios] = this.Plies(this.Ply_Material_Matrix(k)).TestFailure(theta, epsilon_Lam, RunningLoads, failureCriteria, k);

                    % save failure data
                    failedPlies(k) = failedPlies(k) || hasFailed; % if it has fialed and it had not already
                    fails = [fails, Failures];
                    failRats(k,i,:) = FailureRatios;
                end
            end

        end


        % Getters and setters

        % CLT Matricies
        function A = get.A(this)
            A = this.laminateStiffness(1:3,1:3);
        end

        function B = get.B(this)
            B = this.laminateStiffness(1:3,4:6);
        end

        function D = get.D(this)
            D = this.laminateStiffness(4:6,4:6);
        end

        function this = set.A(this, A)
            error("setting A not implemented");
        end

        function this = set.B(this, B)
            error("setting B not implemented");
        end

        function this = set.D(this, D)
            error("setting D not implemented");
        end

    end
end