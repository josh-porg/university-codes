function HW_1_func()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
clear; close; clc;
dim = 10e8; %10^8 is generally pretty good 
a = zeros(dim,1);
b = ones(dim,1);
c = zeros(dim,1);

tic
c = a + b;
toc

tic
c = a - b;
toc

tic
c = a .* b;
toc

tic
c = b .\ a;
toc

tic
c = sin(a);
toc


tic
c = log(a);
toc

disp(c)
end

