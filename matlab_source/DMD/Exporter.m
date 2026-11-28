function Exporter(caseName,variableNames)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
fileSuffixes = ["results.mat";"reconstruction.mat"];
fileContents = ["modes","reconstruction"];

names = strings(length(variableNames),length(fileSuffixes));

for i = 1:length(variableNames)
    for j = 1:length(fileSuffixes)
        names(i,j) = caseName + variableNames(i) + fileSuffixes(j);
    end
end

for i = 1:length(variableNames)
        saveModesForExport(caseName,variableNames(i));
end

for i = 1:size(names,1)
    for j = 1:size(names,2)
        export_snapshots(names(i,j), caseName+"_"+variableNames(i)+"_"+fileContents(j));
    end
end
end

