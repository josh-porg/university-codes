
% enablePlottting
enablePlottting = false; % should it p.lot as it goes?
add_noise = true; % add noise to sensor inputs for training?

%%%%%%%%%%
% TODO: used seeded random evasion. store the seed and evaluate the entire
% generation on the same seeds. regegerate seeds for next generation
%%%%%%%%%%

%%%%%%%%%%
%TODO: use len_hist = x(end); and make x one larger than needed
%%%%%%%%%%

% -------------------------------------------------------------------------
% 1. Create a feedforward network with two hidden layers of 32 neurons each
% -------------------------------------------------------------------------
inputsize = 8+20;
layerSizes = [2,2,2,2,2,2,2,2,2,2,2,2].^5;
outputsize = 20;

inputsize = 13;
layerSizes = [40,40,40,40,40];
layerSizes = [9,9,9];
%layerSizes = [15,15,15,15];

outputsize = 4;

if outputsize > 1
    inputsize = inputsize + 2*outputsize-1;
end

nVars = layerSizes(1)*inputsize + layerSizes(1);
for i = 2:length(layerSizes)
    layerSize = layerSizes(i);
    nVars = nVars + layerSizes(i)*layerSizes(i-1) + layerSizes(i); % weights + bias
end
nVars = nVars + outputsize*layerSizes(end) + outputsize; % weights + bias
nVars = nVars+1+1; % add one for the smoothing size, one for explosion angle and 7 of the other smoothing sizes

% nVars = layerSizes(1)*2 + layerSizes(1) ...     % weights1 + bias1
%       + layerSizes(2)*layerSizes(1) + layerSizes(2) ... % weights2 + b2
%       + layerSizes(2) + 1;                       % W3 + b3

% % 4. Set GA bounds (e.g. ±1 around zero)
lb = -1*ones(1, nVars);
ub =  1*ones(1, nVars);

% 5. Create a fitness-function handle (passing extra args)
fitnessFcn = @(x) computeFitness(x, enablePlottting, add_noise);

% 6. Configure GA options
popSize = 60; % how many individual per generation
elitePercent = .10; % save this fraction of population (the best)

opts = optimoptions('ga', ...
    'PopulationSize', popSize, ...
    'MaxGenerations', 70, ...
    'UseParallel', true ...
    , ...
    'EliteCount', ceil(elitePercent*popSize), ...
    'Display', 'iter');

% 7. Run GA!
[xBest, fvalBest] = ga(fitnessFcn, nVars, [], [], [], [], lb, ub, [], opts);

% 8. Build the optimized network

bestNN = @(n,Active) mySimpleNN(xBest, [n;Active], layerSizes);
% netOptimized = setwb(netTemplate, xBest);
% disp(['Optimized max Pd = ', num2str(1 - fvalBest)]);

maxPd = simulateAndPlot(xBest, true, false, false);

save("NN_missile_guidance_2", "xBest")

% plot a few to entertain
for i = 1:50
    maxPd = simulateAndPlot(xBest, true, false, true);
end


function Pd_avg = computeFitness(x, enablePlottting, add_noise)
    Pd_avg = 0;
    n_inters = 25;
    % do it 3 times because the running is random
    for itteration = 1:n_inters
        Pd_avg = Pd_avg + simulateAndPlot(x, enablePlottting, add_noise, false);
    end

    Pd_avg = Pd_avg/n_inters;
    close all;
end
