%Homework 9 Jackson Torok

clc
clear
close all

%% 11.1

n = 1:10E6;
dp = vpa(sum(1./n))
sp = vpa(sum(1./single(n))) % single causes round off errors

%% 11.2

nd = 1:10;
expected = sum(1./nd)

ni = int8(1:10);
ians = sum(1./ni)

imat = 1./ni % lower values round to 0

%% 11.4

doublea = 5+3i;
singlea = single(5+3i);

doublea^100
singlea^100

% when the value in single exceeds the data capacity, inf is assigned. This
% happened before the double because the single has less memory 

%% 11.6

charnum = char('85');
elementcount = length(charnum)

numequiv = double('8')

numequiv = double('5')

%% 11.7

name1 = char('Emeliann  ','Isabel  ','Kayleigh  ','Michelle  ','Matt  ');

bd =[3,6,2002; 2,18, 2002; 3,22,2002; 5,20,2002; 12,4,2001];
bds = num2str(bd);
disp([name1,bds])

%% 11.8

name2 = string({'Emeliann';'Isabel';'Kayleigh';'Michelle';'Matt'});
fprintf(' \n%10s %3.0f %3.0f %5.0f',[name2';bd'])

%% 11.15

A = [1,2;3,4];
B = [10,20;30,40];
C = [3,6;9,12];

ABC(:,:,1) = A;
ABC(:,:,2) = B;
ABC(:,:,3) = C;

disp(ABC)

Column_A1B1C1(:,:) = ABC(:,1,:)

r2 = squeeze([ABC(2,1,:),ABC(2,2,:)]);
Row_A2B2C2 = squeeze(ABC(2,:,:))'

onetwothree = ABC(1,2,3)

%% 11.17

load test_results
data(:,:,1,1) = test1year1;
data(:,:,1,2) = test2year1;
data(:,:,1,3) = test3year1;
data(:,:,2,1) = test1year2;
data(:,:,2,2) = test2year2;
data(:,:,2,3) = test3year2;

s1q2y1t3 = data(1,2,1,3)

s1q1t2 = data(1,1,:,2)

s2t1y2 = data(2,:,2,1)

q3t2 = data(:,3,:,2)

%% 11.19

Al = 'aluminum';
Cu = 'copper';
Fe = 'iron';
Mo = 'molybdenum';
Co = 'cobalt';
elements = {Al,Cu,Fe,Mo,Co}'



