function [snapshotMatrix, meanMatrix] = getsnapshotmatrix(dataMatrix)
% Depricated - getsnapshotmatrix rshapes dataMatrix and removes the temporal mean of the
%data and returns the snapshot matrix
%   takes data matrix of from S(i,j,t)
%   where i is x-pos    , j is y-pos, t is time
%   outputs the snapshot matrix of data - mean

countOfTime = length(dataMatrix(1,1,:));

% reshape
dataMatrix = reshape(permute(dataMatrix,[3,2,1]),countOfTime,[]); % reshapes the datamatrix to have a row for each time

% remove mean
dataMean = mean(dataMatrix,1); %mean calculates mean given array and dimention: 1=comumns
meanMatrix = repmat(dataMean, countOfTime, 1); %repmat repeats the given value into the given arry size thus this create a matrix of the average of a suitable size to subtract form  datamatrix 
snapshotMatrix = dataMatrix - meanMatrix; %removes the mean from each row of data to give snapshots

end
