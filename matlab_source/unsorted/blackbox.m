function out = blackbox(index)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%countOfDigits = 10000


%100 = 479
%50 = 248
%47 = 235
%44 = 214
%40 = 195
%39 = 188
%38 = 179
%37 = 178
%35 = 166
%25 = 118


%6 rows per inch
%inches 10 + 1.5 + 1 + .5 + 16 + 1 = 30
%inches * (row/inch) = 180
digits(index);

qnumber = vpa(pi);

%pidigits = zeros(countOfDigits,1);

charp = char(qnumber);
charp(2) = []; % get rid of the decimal point

qdigits = arrayfun(@str2num,charp);

% for i = 1:countOfDigits
%     %fprintf('digit in position %i = %i\n', i, floor(pinumber));
%     %pidigits(i) = max(1,floor(pinumber));
%     pidigits(i) = floor(pinumber);
%     pinumber = pinumber - floor(pinumber);
%     pinumber = pinumber*10;
% end

% pididgitsums(1) = zeros(countOfDigits,0);
% pididgitsums(1) = 3;
% for i=2:countOfDigits
%     pididgitsums(i) = pididgitsums(i-1) + pidigits(i);
% end

qsum = sum(qdigits);
%plot(pidigits,'o');
out = qsum;
clear('qsum','qdigits','charp','qnumber','digits');
end

