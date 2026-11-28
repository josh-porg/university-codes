classdef DynamicModeDecomposer
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    %   maybe offer exact vs projected option
    %   consider taking the pod modes by with method of snapshots as seen
    %   on page 29 of dmd book (p.42 of pdf)
    
    properties
        truncationPercentage = .5 / 100
        
        rank
        
        singularValues
        
        truncatedSingularValues
        
        dataMatrix %may be better just to take a snapshot matrix input
        
        meanMatrix
        
        countOfTime %depricted and should be avoided
        
        snapshotMatrix
        
        dataDimentions % original dimentions of raw datamatrix
        
        modes %DMD modes of the data (often referd to as PHI)
        
        descreteTimeEigenValues % descrete time eigenValues (often refered to a lambda)
        
        continuousTimeEigenValues % continuous time eigenValues OMEGA = log(LAMBDA)/delta_t
        
        b % initial magnitude of the modes given by initial condition
        
        delta_t % dt = time step advancing X1 to X2 (X to X’
        
        timeDynamics
        
        reconstructedSnapshots
        
        reconstructedDataMatrix
        
        backgroundModes
        
        backgroundEigenValues %continuous time eig vals
        
        forgroundModes
        
        forgroundEigenValues %continuous time eig vals
        
        reconstructedBackground
        
        reconstructedBackgroundData
        
        reconstructedForground
        
        reconstructedForgroundData
        
        bBackground
        
        bForground
        
        backgroundTimeDynamics
        
        forgroundTimeDynamics
        
        architecture
        
    end
    
    methods
        function this = DynamicModeDecomposer(rawDataMatrix,options)
            %DynamicModeDecomposer Construct an instance of this class
            %   Detailed explanation goes here
            % dt = time step advancing X1 to X2 (X to X’
            arguments
                rawDataMatrix;
                options.rank {mustBeScalarOrEmpty} = -1;
                options.delta_t {mustBeScalarOrEmpty} = 1;
                options.backgroundSeperation {mustBeNumericOrLogical} = false;
                options.architecture = 'standard';
                
            end
            
            if ndims(rawDataMatrix) == 2 %it is or can be treated as a snapshot matrix
                %this.snapshotMatrix = tall(rawDataMatrix);
                this.snapshotMatrix = rawDataMatrix;
                this.countOfTime = size(rawDataMatrix,2);
            else %it's acutally data
                
                this.dataMatrix = rawDataMatrix; %save the original data matrix (see if this can be avoided) %range 0-255
                this.dataDimentions = size(rawDataMatrix); %save size of raw data matrix

                %calculate the count of time
                %could also use the number of colomns in snapshotmatrix
                if length(this.dataDimentions) == 3 %3D case
                    this.countOfTime = length(rawDataMatrix(1,1,:)); %calculate the count of time
                elseif length(this.dataDimentions) == 4 %4D case
                    this.countOfTime = length(rawDataMatrix(1,1,1,:)); %calculate the count of time
                else
                    throw(MException("datamatrix not regonised as having 3 or 4 dimentions"));
                end
                
                this.snapshotMatrix = this.data2snapshots(this.dataMatrix); %convert the data into snapshots %range 0-255
                %this.snapshotMatrix; this.meanMatrix = this.removeMean(this.snapshotMatrix); %remove the mean %snapshotmatrix 0-255 % mean matrix 0-255
                
            end
            
            %make sure rank is reasonable
            if options.rank <= 0
                this.rank = intmax('uint64');
                %display('negative or 0 rank was provided. using full rank.');
            else
                this.rank = options.rank;
            end
            
            
            this.delta_t = options.delta_t;
            
            
            %do the dmd
            this = this.calculateModes(this.rank);
            
            this.modes; %work out why this is here and/or remove
            
            firstSnapshot = this.snapshotMatrix(:,1); %first snapshot range = 252
            %possibly try restroing the mean before setting intial
            %condition (which will probably go horribly wrong) so maybe
            %never remove mean in the first place
            
            %DIFFERENCE FROM THERE MOTHOD THEY USE B = PHI\X1 
            %however in there maths they say that morre penrose spudoinverse and multiplication should work
            %this.b = this.modes\firstSnapshot; %this should also work
            this.b = pinv(this.modes)*firstSnapshot; % range real (this.b) = -1.2425*10^-4 -> 2.0370*10^3
            
            %OMEGA  =log(LAMBDA)/delta_t
            this.continuousTimeEigenValues = log(this.descreteTimeEigenValues)/this.delta_t;
            %real omega min = -.7552 max = -.0027
            
            if options.backgroundSeperation == true
                this = this.seperateBackground();
            end
            
            this.architecture = options.architecture;
            
            this = this.reconstruct();
            
            
        end
        
        function this = reconstruct(this)
            %RECONSTRUCT reconstructs on modes
            
            seperateBackground = ~isempty(this.backgroundModes) && ~isempty(this.forgroundModes);
            
            this.timeDynamics = zeros(this.rank, this.countOfTime-1);
            t = (0:this.countOfTime-2) * this.delta_t; % -2 because we start at t = 0
            for i = (1:this.countOfTime-1)
                this.timeDynamics(:,i) = (this.b .* exp(this.continuousTimeEigenValues * t(i))); % = e^(OMEGA*t) * b
            end
            %reconstruct on modes phi
            this.reconstructedSnapshots = this.modes * this.timeDynamics; % = PHI * e^(OMEGA*t) * b
            %add the mean back in
%             this.reconstructedSnapshots = this.restoreMean(this.reconstructedSnapshots,this.meanMatrix(:,1:end-1));
            
            
            if seperateBackground
                this.backgroundTimeDynamics = zeros(numel(this.backgroundEigenValues), this.countOfTime-1);
                this.forgroundTimeDynamics = zeros(numel(this.forgroundEigenValues), this.countOfTime-1);
                
                firstSnapshot = this.snapshotMatrix(:,1);
                
                this.bBackground = pinv(this.backgroundModes) * firstSnapshot;
                this.bForground = pinv(this.forgroundModes) * firstSnapshot;
                
                t = (0:this.countOfTime-2) * this.delta_t;
                
                for i = 1:this.countOfTime-1
                    %these may need to be .* elementwise multiplication
                    this.backgroundTimeDynamics(:,i) = this.bBackground .* exp(this.backgroundEigenValues .* t(i));
                    this.forgroundTimeDynamics(:,i) = this.bForground .* exp(this.forgroundEigenValues .* t(i));
                end
                
                this.reconstructedBackground = this.backgroundModes * this.backgroundTimeDynamics;
                this.reconstructedForground = this.forgroundModes * this.forgroundTimeDynamics;
            end
            
            %convert back into datamatrix
            if isempty(this.dataMatrix)
                %do nothing
            elseif length(this.dataDimentions) == 3 %3D case
                this.reconstructedDataMatrix = this.snapshots2data(this.reconstructedSnapshots, [this.dataDimentions(1),this.dataDimentions(2),this.dataDimentions(3)-1]);
                if seperateBackground
                    this.reconstructedBackgroundData = this.snapshots2data(this.reconstructedBackground, [this.dataDimentions(1),this.dataDimentions(2),this.dataDimentions(3)-1]);
                     this.reconstructedForgroundData = this.snapshots2data(this.reconstructedForground, [this.dataDimentions(1),this.dataDimentions(2),this.dataDimentions(3)-1]);
                end
            elseif length(this.dataDimentions) == 4 %4D case
                this.reconstructedDataMatrix = real(this.snapshots2data(this.reconstructedSnapshots, [this.dataDimentions(1),this.dataDimentions(2),this.dataDimentions(3),this.dataDimentions(4)-1]));
                if seperateBackground
                    this.reconstructedBackgroundData = real(this.snapshots2data(this.reconstructedBackground, [this.dataDimentions(1),this.dataDimentions(2),this.dataDimentions(3),this.dataDimentions(4)-1]));
                    this.reconstructedForgroundData = real(this.snapshots2data(this.reconstructedForground, [this.dataDimentions(1),this.dataDimentions(2),this.dataDimentions(3),this.dataDimentions(4)-1]));
                end
            else
                throw(MException("datamatrix not regonised as having 3 or 4 dimentions"));
            end
        end
        
        function modes = getModes(this, rank, format) %TODO GET RID OF RANK ARGUMENT
            %GETMODES rank number of modes, in format 'data' or 'snapshot',
            %default is data
            
            if nargin < 3 && ~isempty(this.dataDimentions)
                format = "data";
            elseif nargin < 3
                format = "snapshot";
            end
            
            %{
            %make sure rank is reasonable
            if rank <= 0
                rank = intmax('uint64');
                display('negative or 0 rank was provided. using full rank.');
            end
            %}
            
            if isempty(this.modes)
                disp("modes not previously calculated. they have been calculated but not stored")
                this = this.calculateModes(rank);
                
            end
            
            if format == "snapshot"
                modes = this.modes;
            else %currently only 3d case
                
                %disp(size(this.modes)) % for debugging
                
                if length(this.dataDimentions) == 3 %currently only 3d case
                    modes = zeros(this.dataDimentions(1),this.dataDimentions(2),size(this.modes,2)); %preallocate mdoes to prevent resizing in loop
                    modes(:,:,1) = this.snapshots2data(this.modes(:,1),[this.dataDimentions(1),this.dataDimentions(2),1]);
                    
                    for i = 1:size(this.modes,2)
                        %display('itteration %i', i);
                        modes(:,:,i) = this.snapshots2data(this.modes(:,i),[this.dataDimentions(1),this.dataDimentions(2),1]);
                    end
                    
                elseif length(this.dataDimentions) == 4 %4D case
                    modes = zeros(this.dataDimentions(1),this.dataDimentions(2),this.dataDimentions(3),size(this.modes,2)); %preallocate mdoes to prevent resizing in loop
                    modes(:,:,:,1) = this.snapshots2data(this.modes(:,1),[this.dataDimentions(1),this.dataDimentions(2),this.dataDimentions(3),1]);
                    
                    for i = 1:size(this.modes,2)
                        %display('itteration %i', i);
                        modes(:,:,:,i) = this.snapshots2data(this.modes(:,i),[this.dataDimentions(1),this.dataDimentions(2),this.dataDimentions(3),1]);
                    end
                else
                    throw(MException("datamatrix not regonised as having 3 or 4 dimentions"));
                end
                
            end
        end
        
        function this = calculateModes(this, rank)
            %CALCULATE MODES
            [this.modes, this.descreteTimeEigenValues, this.singularValues, this.rank] = this.decompose(rank); % do decomposition to get
            %if rank was set to infinity update it now to be the full rank
            %of the system
            this.truncatedSingularValues = this.singularValues(1:this.rank); %needs to be assigned
            
        end
            
        function [modes, descreteTimeEigenValues, singularValues, rank] = decompose(this, rank)
            %DECOMPOSE decomposes data, using trunkates number of modes to
            %rank
            %   get two snapshotmatrices x and xprime such that x includes
            %   snapshots 1 -> n-1 and xprime includes snapshots 2 -> n.
            %   compute the SVD of x (because we are about to take the
            %   psudo inverse. create the small matrix s_tilde that
            %   approximates the much larger matrix A. 
            %   S_tilde = U'xprimeVSigma' transpose. solve the eignen problem
            %   AW = WLAMBDA
            
            if rank <= 0
                exception = MException('value of rank (%i) was 0, or negative.', rank);
                throw(exception);
            end
            
            
            x = this.snapshotMatrix(:,1:end-1); % 0-255
            xprime = this.snapshotMatrix(:,2:end); % 0-255
            
            %switch this.architecture;
            %case 'standard'
            [U, SIGMA, V] = svd(x,'econ'); 
            %[U, SIGMA, V] = gather(U, SIGMA, V);
            singularValues = diag(SIGMA);
            
            %try comparing with matlab's rank function
            
            %define cutoff for zero
            effectiveZeroCutoff = sum(SIGMA,'all')*this.truncationPercentage; %might be slightly faster if i use singular values instead of sigma
            %find at what point the values of SIGMA go to 0
            nonzeroSingularValueIndexes = singularValues > effectiveZeroCutoff;
            %SingularValues = singularValues(nonzeroSingularValueIndexes); % we don't want to do this because we wuld rather sae the information
            SIGMA = SIGMA(nonzeroSingularValueIndexes,nonzeroSingularValueIndexes);
            
            %%%%%%%%%%% TODO: USE NNZ TO FINE NON ZERO VALUES
            %find function will find non zero indecies.
            
            U = U(:,nonzeroSingularValueIndexes);
            V = V(:,nonzeroSingularValueIndexes);
            autoTruncateRank = sum(nonzeroSingularValueIndexes,'all');
            
            
            %SIGMA_lowRank = SIGMA(nonzeroSingularValueIndexes); %cpuld be
            %clever and do this. may even be more clever to do this on the
            %diag of Sigma stored as matrix "singular values" then to use
            %that truncate each matrix
            
            
            
            %maybe make this an if statment and provide a debugging report
            rank = min(min(rank, size(U,2)), autoTruncateRank); % if we get fewer modes than the requested rank, use the number of modes
            
            %truncate matixes at the desired rank
            U_lowRank = U(:,1:rank); % max = .1214 min = -.0969 % same max min relative error .3426
            SIGMA_lowRank = SIGMA(1:rank,1:rank); % max = 1.047 *10^5 min = 0 %same range as full
            V_lowRank = V(:,1:rank); % dont't be confused by this being V and not V' % max = .5447 min = -.5759 %some loss of range (<1/6)
            
            %singularValues = diag(SIGMA_lowRank);
            
            %continue with dmd
            % s_tilde = U' * xprime * V / SIGMA; % full thing % s_tilde 3579^2; all real
            s_tilde = U_lowRank' * xprime * V_lowRank / SIGMA_lowRank; %s tilde gets negative value which gives complex eigs
            %s_tilde = U_lowRank' * xprime * V_lowRank * (1./SIGMA_lowRank); %this doesnt idk why not %because SIGMA is singular we can invert it
            %
            %TODO: use a cuttoff for 0
            
            %end
            
            %s_tilde range = 3.2467
            
            % find eigenvalues and eigenvectors
            [W,LAMBDA] = eig(s_tilde); % i might want to specify vector % consider just calculating the larges few with eigs
            % eig of s tilde is returning complex numbers for both w and
            % lambda
            
            %range W ~= 1
            %range LAMBDA ~= 1
            
            %lambda is diagonal matrix
            % use diag(LAMBDA) to get descreteTimeEigenValues
            
            modes = xprime * V_lowRank / SIGMA_lowRank * W; %exact dmd modes
            %apparently faster to use division rather than inverse
            %modes = xprime * V_lowRank * inv(SIGMA_lowRank) * W; %dmd modes
            
            %range modes ~= .05 + .05i
            
            %order from largest to smallest
            
            % sort the eigen vectors by size and create a corresponding 
            % change in idex so we can perform the same switch to the eigen vectors
            %%%[W, WIndex] = sort(W,"descend"); 
            
            descreteTimeEigenValues = diag(LAMBDA);
            %%%LAMBDA = LAMBDA(WIndex);
            
            %{
            modes = U * LAMBDA;
            
            mode = zeros(this.dataDimentions(1),this.dataDimentions(2));
            for i = 1:quantity
                mode(i) = U * LAMBDA(i);
            end
            %}
            
            %dmd modes are 
            %modes = xprime*v*sigma
        end
        
        function this = seperateBackground(this)
            background = find(abs(this.continuousTimeEigenValues) < 1e-2); % find modes which aproximatly no change in time  
            forground = setdiff(1:this.rank, this.continuousTimeEigenValues); % get all other modes
            
            this.backgroundEigenValues = this.continuousTimeEigenValues(background);
            this.backgroundModes = this.modes(:,background);
            
            this.forgroundEigenValues = this.continuousTimeEigenValues(forground);
            this.forgroundModes = this.modes(:,forground);
        end
        
    end
    
    methods (Static)
        
        function snapshotMatrix = data2snapshots(dataMatrix)
            %convertToSnapshotMatrix reshapes a dataMatrix into a snapshot
            %matrix
            %data and returns the snapshot matrix
            %   takes data matrix of from S(i,j,t)
            %   where i is x-pos    , j is y-pos, t is time
            %   outputs the snapshot matrix in which each coloumn
            %   represents a snapshot at time t
            
            %TODO: generalize this to work on n-dimentional data
            
            %note countOfTime = dataShape(ndims(datamatrix)) 
            dimentionality = ndims(dataMatrix); % the number of dimention of the data matrix
            countOfTime = size(dataMatrix,ndims(dataMatrix)); %time is assumed to be the last dimention of data matrix
            
            if (dimentionality == 3) % data is a 3d matrix 
                snapshotMatrix = reshape(permute(dataMatrix,[3,2,1]), countOfTime, [])'; % reshapes the datamatrix to have a row for each time
            elseif (dimentionality == 4) % data is a 4d matrix 
                snapshotMatrix = reshape(permute(dataMatrix,[4,1,2,3]), countOfTime, [])'; % reshapes the datamatrix to have a row for each time
            else
                throw(MException("datamatrix not regonised as having 3 or 4 dimentions"));
            end
            
        end
        
        function dataMatrix = snapshots2data(snapshotMatrix, dataDimentions)
            %TODO: consider using hilbert curve
            snapshotMatrix = snapshotMatrix';
            
            dimentionality = size(dataDimentions,2);
            
            %undo transformations done in data2snapshots
            if (dimentionality == 3) % data is a 3d matrix 
                dataMatrix = permute(reshape(snapshotMatrix,dataDimentions(3),dataDimentions(2),[]),[3,2,1]);
            elseif (dimentionality == 4) % data is a 4d matrix 
                dataMatrix = ipermute(reshape(snapshotMatrix,dataDimentions(4),dataDimentions(1),dataDimentions(2),[]),[4,1,2,3]);
            else
                throw(MException("data matrix not regonised as having 3 or 4 dimentions"));
            end
            
        end
        
        function [meanFreeSnapshotMatrix, meanMatrix] = removeMean(snapshotMatrix)
            % remove mean
            
            %ineffetive method
            %{
            dataMean = mean(snapshotMatrix,1); %mean calculates mean given array and dimention: 1=columns
            meanMatrix = repmat(dataMean, size(snapshotMatrix,1), 1); %repmat repeats the given value into the given arry size thus this create a matrix of the average of a suitable size to subtract form  datamatrix
            meanFreeSnapshotMatrix = snapshotMatrix - meanMatrix; %removes the mean from each row of data to give snapshots
            %}
            
            average = mean(snapshotMatrix,2);
            meanMatrix = average * ones(1, size(snapshotMatrix, 2));
            meanFreeSnapshotMatrix = snapshotMatrix - meanMatrix;
            
        end
        
        function meanRestoredSnapshotMatrix = restoreMean(snapshotMatrix, meanMatrix)
            % size(snapshotMatrix)% for debug
            % snapshotMatrix % for debug
            % size(meanMatrix) % for debug
            % meanMatrix % for debug
            meanRestoredSnapshotMatrix = snapshotMatrix + meanMatrix; %reintroduces the mean from each row of data to give snapshots
        end
        
    end
    
    methods (Access = private)
        function foo(this)
            %Temporary to maintain syntacticly correct class while being
            %written and tested. TODO: remove when othe private methods are
            %introduced
        end
    end
end

