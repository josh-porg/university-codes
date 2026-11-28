function saveModesForExport(caseName,variableName)
%saveModesForExport saves a file of the modes in a programaticly exportable
%format
%   Detailed explanation goes here
%saveModes for export
load(caseName+variableName+"results.mat","modes");
%load("dmd700000"+variable+"results.mat",'modes');
p = modes;
save(caseName+variableName+"_modes.mat","p","-v7.3");
end



% %saveModes for export
% variable = 'P';
% load('baz.mat');
% %load("dmd700000"+variable+"results.mat",'modes');
% p = modes;
% save("dmd700000"+variable+"modes.mat","p","-v7.3");