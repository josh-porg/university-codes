%Homework 5 Jackson Torok

clc
clear
close all

%% 6.1

n=10:100;

numgrain=num_grain(n);

plot(n,numgrain)
title('Number of grains in a x100 mangified square inch of metal VS ASTM Grain Size')
xlabel('Grain Size')
ylabel('Number of Grains')

%% 6.6

r=[7926,4217];

h=[0:100:10000]./5280;

dist=distance(r,h);

tab=table(dist(:,1),dist(:,2),'VariableNames',{'Horizon Distance on Earth (mi)','Horizon Distance on Mars (mi)'});

%% 6.7
figure(2)

t=0:.5:30;

hi=height(t);

plot(t,hi)

[~,ind]=max(hi);

tfall=t(ind)

%% 6.14
figure(3)

height_handle=@height;

fplot(height_handle,[0,60])

tground=fzero(height_handle,30)

