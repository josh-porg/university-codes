classdef ProperOrthogonalDecomposer
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        maxEigs %TODO: implement
        
        rawDataMatrix
        
        dataMatrix
        
        meanMatrix
        
        countOfTime
        
        snapshotMatrix
        covarianceMatrix
        spatialModes
        timeCoefficients
        eigenValues
        
        rawDataDimentions % original dimentions of raw datamatrix
    end
    
    methods
        function this = ProperOrthogonalDecomposer(rawDataMatrix)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            
            %save the data
            this.rawDataMatrix = rawDataMatrix;
            this.rawDataDimentions = rawDataMatrix.size(); %save size of raw data matrix
%             % reshape and store as dataMatrix
%             % currently done in Getsnapshot
%             this.dataMatrix = reshape(permute(this.rawDataMatrix,[3,2,1]),this.countOfTime,[]); % reshapes the datamatrix to have a row for each time
            
            %calculate the count of time
            this.countOfTime = length(rawDataMatrix(1,1,:));
            
            %remove mean and convert to snapshots
            this = this.getSnapshotMatrix();
            %this = this.getSnapshotMatrix("Y");
            
            %calculate covariave matrix
            this = this.getCovarianceMatrix();
            
            %create matrix of spatial modes phi of eigen vectors ordered by size of eigen values
            this = this.getSpatialModes();
            
            %multipy U' by phi to get time coefficients A
            this = this.getTimeCoefficients();
            
        end
        
%         function this = ProperOrthogonalDecomposer(rawDataMatrix, maxEigs)
%             %UNTITLED3 Construct an instance of this class
%             %   Detailed explanation goes here
%             
%             %save the data
%             this.rawDataMatrix = rawDataMatrix;
%             
% %             % reshape and store as dataMatrix
% %             % currently done in Getsnapshot
% %             this.dataMatrix = reshape(permute(this.rawDataMatrix,[3,2,1]),this.countOfTime,[]); % reshapes the datamatrix to have a row for each time
%             
%             %calculate the count of time
%             this.countOfTime = length(rawDataMatrix(1,1,:));
%             
%             %remove mean and convert to snapshots
%             this = this.getSnapshotMatrix();
%             
%             %calculate covariave matrix
%             this = this.getCovarianceMatrix();
%             
%             %create matrix of spatial modes phi of eigen vectors ordered by size of eigen values
%             this = this.getSpatialModes();
%             
%             %multipy U' by phi to get time coefficients A
%             this = this.getTimeCoefficients();
%             
%         end
        
        
        function U_tilde_k = reconstruct(this,mode)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            % reconstruct
            U_tilde_k = (this.timeCoefficients(:,mode) * this.spatialModes(:,mode)');
        end
        
        function restoredMatrix = restoreFormat(this, matrixOfSnapshots)
            restoredMatrix = permute(reshape(matrixOfSnapshots,size(this.rawDataMatrix,3),size(this.rawDataMatrix,2),[]),[3,2,1]);

        
        end
        
        function meanRestoredSnapshotMatrix = getMeanRestoredSnapshotMatrix(this, modeMatrix)
            meanRestoredSnapshotMatrix = modeMatrix + this.meanMatrix; %reintroduces the mean from each row of data to give snapshots
        end
        
        function print(this)
            "rawDataMatrix"
            this.rawDataMatrix
            "dataMatrix"
            this.dataMatrix
            "meanMatrix"
            this.meanMatrix
            "countOfTime"
            this.countOfTime
            "snapshotMatrix"
            this.snapshotMatrix
            "covarianceMatrix"
            this.covarianceMatrix
            "spatialModes"
            this.spatialModes
            "timeCoefficients"
            this.timeCoefficients
            "eigenValues"
            this.eigenValues
        end
        
    end
    
    methods (Static)        
%         function outputArg = composeData(dataMatrices)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             
%         end
    end
    
    methods (Access = private)
        
        function this = getSnapshotMatrix(this)
            %getSnapshotMatrix rshapes dataMatrix and removes the temporal mean of the
            %data and returns the snapshot matrix
            %   takes data matrix of from S(i,j,t)
            %   where i is x-pos    , j is y-pos, t is time
            %   outputs the snapshot matrix of data - mean
            
            
            % reshape and store as dataMatrix
            
            if (ndims(data) == 3) %this is a 3d array 
                this.dataMatrix = reshape(permute(this.rawDataMatrix,[3,2,1]),this.countOfTime,[]); % reshapes the datamatrix to have a row for each time
            elseif (ndims(data) == 4)
                this.dataMatrix = reshape(permute(this.rawDataMatrix,[4,1,2,3]),this.countOfTime,[]); % reshapes the datamatrix to have a row for each time
            else
                throw(MException("datamatrix not regonised as having 3 or 4 dimentions"));
            end
            
            
            % remove mean
            dataMean = mean(this.dataMatrix,1); %mean calculates mean given array and dimention: 1=columns
            this.meanMatrix = repmat(dataMean, this.countOfTime, 1); %repmat repeats the given value into the given arry size thus this create a matrix of the average of a suitable size to subtract form  datamatrix
            this.snapshotMatrix = this.dataMatrix - this.meanMatrix; %removes the mean from each row of data to give snapshots
            
        end
        
        %TODO move to static private
%         function meanRestoredSnapshot = restoreMean(matrixOfSnapshots, mean)
%             meanRestoredSnapshot = matrixOfSnapshots + mean; %reintroduces the mean from each row of data to give snapshots
%         end
        
        function this = getCovarianceMatrix(this)
            %calculate covariave matrix
            this.covarianceMatrix = (1 / (length(this.rawDataMatrix))) * (this.snapshotMatrix'*this.snapshotMatrix); 
        end
        
        function this = getSpatialModes(this)
            %find eigenvalues and eigenvectors
            [this.spatialModes, LAMBDA] = eig(this.covarianceMatrix, 'vector');
            %order from largest to smallest
            [lambda, lambdaIndex] = sort(LAMBDA, "descend");
            this.eigenValues = lambda;
            %create matrix of spatial modes phi of eigen vectors ordered by size of eigen values
            this.spatialModes = this.spatialModes(:,lambdaIndex);
        end
        
        function this = getTimeCoefficients(this)
            %multipy U' by phi to get time coefficients A
            this.timeCoefficients = this.snapshotMatrix*this.spatialModes;
        end
        
    end
end

