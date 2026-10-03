load('NN_missile_guidance.mat');

pd = zeros([1,50]);
for i = 1:300
    Pd(i) = simAndPlot_requested(xBest, true, false, false);
end

mean(pd)