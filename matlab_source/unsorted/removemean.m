function [snapshotMatrix] = removemean(dataMatrix)
%removeMean removes the mean of the data
%   takes data matrix of from S(i,j,t)
%   where i is x-pos, j is y-pos, t is time
%   outputs the snapshot matrix of data - mean

countOfTime = length(dataMatrix(1,1,:));

dataMean = mean(dataMatrix,1); %mean calculates mean given array and dimention: 1=comumns
meanSubtractionMatrix = repmat(); %repmat repeats the given value into the given arry size
snapshotMatrix = -1;
end

