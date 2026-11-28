% Enable plotting during training
enablePlotting = false; % Plot during training
add_noise = true; % Add noise to sensor inputs for training

% Neural Network Architecture
inputSize = 6; % Number of input features
layerSizes = [25, 25, 25]; % Hidden layer sizes
layerSizes = [7, 7, 7]; % Hidden layer sizes
outputSize = 1; % Number of outputs (first is turn rate command)

% Create the RL environment
env = MissileGuidanceEnv(enablePlotting, add_noise, inputSize, layerSizes, outputSize);

% Create the actor network (policy)
actorNetwork = [
    featureInputLayer(inputSize, 'Normalization', 'none', 'Name', 'state')
    % fullyConnectedLayer(layerSizes(1), 'Name', 'fc1')
    % reluLayer('Name', 'relu1')
    % fullyConnectedLayer(layerSizes(2), 'Name', 'fc2')
    % reluLayer('Name', 'relu2')
    % fullyConnectedLayer(layerSizes(3), 'Name', 'fc3')
    % reluLayer('Name', 'relu3')
    fullyConnectedLayer(layerSizes(1), 'Name', 'fc1')
    reluLayer('Name', 'tanh0')
    fullyConnectedLayer(layerSizes(2), 'Name', 'fc2')
    reluLayer('Name', 'tanh2')
    fullyConnectedLayer(layerSizes(3), 'Name', 'fc3')
    reluLayer('Name', 'tanh3')
    fullyConnectedLayer(outputSize, 'Name', 'output')
    tanhLayer('Name', 'tanh1') % Bounded between -1 and 1
];

%actorOptions = rlRepresentationOptions('LearnRate', 1e-4, 'GradientThreshold', 1);

actor = rlContinuousDeterministicActor(actorNetwork, env.getObservationInfo, env.getActionInfo);

% Create the critic network (value function)
criticNetwork = [
    featureInputLayer(inputSize, 'Normalization', 'none', 'Name', 'state')
    fullyConnectedLayer(layerSizes(1), 'Name', 'fc1')
    reluLayer('Name', 'relu1')
    fullyConnectedLayer(layerSizes(2), 'Name', 'fc2')
    reluLayer('Name', 'relu2')
    fullyConnectedLayer(layerSizes(3), 'Name', 'fc3')
    reluLayer('Name', 'relu3')
    fullyConnectedLayer(1, 'Name', 'output')
];

%new
% Create critic network with proper input layers
statePath = [
    featureInputLayer(inputSize, 'Name', 'stateInput', 'Normalization', 'none')
    fullyConnectedLayer(128, 'Name', 'stateFC1')
    reluLayer('Name', 'stateRelu1')
    fullyConnectedLayer(64, 'Name', 'stateFC2')
];

actionPath = [
    featureInputLayer(outputSize, 'Name', 'actionInput', 'Normalization', 'none')
    fullyConnectedLayer(64, 'Name', 'actionFC1')
];

commonPath = [
    concatenationLayer(1, 2, 'Name', 'concat')
    fullyConnectedLayer(64, 'Name', 'commonFC1')
    reluLayer('Name', 'commonRelu1')
    fullyConnectedLayer(32, 'Name', 'commonFC2')
    reluLayer('Name', 'commonRelu2')
    fullyConnectedLayer(1, 'Name', 'output')
];

% Assemble the network
criticNetwork = layerGraph(statePath);
criticNetwork = addLayers(criticNetwork, actionPath);
criticNetwork = addLayers(criticNetwork, commonPath);

% Connect the layers
criticNetwork = connectLayers(criticNetwork, 'stateFC2', 'concat/in1');
criticNetwork = connectLayers(criticNetwork, 'actionFC1', 'concat/in2');

% Plot to verify network structure (optional)
% plot(criticNetwork)

% Create critic representation
criticOpts = rlRepresentationOptions('LearnRate',1e-3,'GradientThreshold',1);
critic = rlQValueFunction(criticNetwork,...
    env.getObservationInfo(),...
    env.getActionInfo(),...
    'Observation',{'stateInput'},...
    'Action',{'actionInput'});
%new

%criticOptions = rlRepresentationOptions('LearnRate', 1e-4, 'GradientThreshold', 1);
%critic = rlValueRepresentation(criticNetwork, env.getObservationInfo, 'Observation', {'state'}, criticOptions);
%critic = rlQValueFunction(criticNetwork, env.getObservationInfo, env.getActionInfo, 'Observation', {'state'}, 'Action', {'Action'});

% Create DDPG agent
agentOptions = rlDDPGAgentOptions(...
    'SampleTime', 0.1, ...
    'TargetSmoothFactor', 1e-3, ...
    'ExperienceBufferLength', 1e6, ...
    'MiniBatchSize', 128, ...
    'DiscountFactor', 0.99, ...
    'NoiseOptions', rl.option.GaussianActionNoise(...
        'StandardDeviation', 0.3, ...
        'StandardDeviationDecayRate', 1e-4, ...
        'StandardDeviationMin', 0.1));

agent = rlDDPGAgent(actor, critic, agentOptions);

% Training options
trainOpts = rlTrainingOptions(...
    'MaxEpisodes', 1000, ...
    'MaxStepsPerEpisode', 6001, ...
    'ScoreAveragingWindowLength', 20, ...
    'StopTrainingCriteria', 'AverageReward', ...
    'StopTrainingValue', -0.1, ...
    'SaveAgentCriteria', 'EpisodeReward', ...
    'SaveAgentValue', -0.2, ...
    'SaveAgentDirectory', 'savedAgents', ...
    'Plots', 'training-progress', ...
    'Verbose', true);

% ===== 1. DDPG Agent Options =====
agentOptions = rlDDPGAgentOptions(...
    'SampleTime', 0.1, ...
    'TargetSmoothFactor', 1e-3, ...       % Tau for target network updates
    'ExperienceBufferLength', 1e6, ...    % Replay buffer size
    'MiniBatchSize', 128, ...             % Batch size for learning
    'DiscountFactor', 0.99, ...           % Gamma
    'NoiseOptions', rl.option.GaussianActionNoise(...
        'StandardDeviation', 0.5, ...      % Higher initial exploration
        'StandardDeviationDecayRate', 1e-5,... % Slower decay
        'StandardDeviationMin', 0.01), ...  % Minimum noise
    'ActorOptimizerOptions', rlOptimizerOptions(...
        'LearnRate', 1e-4, ...            % Actor learning rate
        'GradientThreshold', 1.0), ...     % Clip gradients
    'CriticOptimizerOptions', rlOptimizerOptions(...
        'LearnRate', 1e-3, ...            % Critic learning rate (higher than actor)
        'GradientThreshold', 1.0));

% ===== 2. Training Options =====
trainOpts = rlTrainingOptions(...
    'MaxEpisodes', 1000, ...
    'MaxStepsPerEpisode', 60000, ...
    'ScoreAveragingWindowLength', 20, ...
    'StopTrainingCriteria', 'AverageReward', ...
    'StopTrainingValue', -0.1, ...        % Stop if avg reward > -0.1
    'SaveAgentCriteria', 'EpisodeReward', ...
    'SaveAgentValue', -0.2, ...           % Save agents with reward > -0.2
    'SaveAgentDirectory', 'savedAgents', ...
    'Plots', 'training-progress', ...
    'Verbose', true);

% ===== 3. Train and Save =====
agent = rlDDPGAgent(actor, critic, agentOptions);
trainingStats = train(agent, env, trainOpts);
save('trainedMissileAgent.mat', 'agent');

% ===== 4. Evaluation =====
simOpts = rlSimulationOptions('MaxSteps', 6001);
experiences = sim(env, agent, simOpts);

% ===== 5. Visualization =====
for i = 1:10
    simulateAndPlotRL(env, agent, true, false)
end

% % Train the agent
% trainingStats = train(agent, env, trainOpts);
% 
% % Save the trained agent
% save('trainedMissileAgent.mat', 'agent');
% 
% % Evaluate the trained agent
% simOpts = rlSimulationOptions('MaxSteps', 600);
% experiences = sim(env, agent, simOpts);
% 
% % Visualize some simulations
% for i = 1:10
%     simulateAndPlotRL(env, agent, true, false)
% end

%% Helper function to simulate and plot with RL agent
function maxPd = simulateAndPlotRL(env, agent, enablePlotting, add_noise)
    % Reset the environment
    obs = reset(env);
    done = false;
    maxPd = 0;
    
    % Run simulation
    while ~done
        % Get action from agent
        action = getAction(agent, obs);

        if isa(action(1), "cell")
                action = cell2mat(action);
        end
        
        % Step the environment
        [nextObs, reward, done, info] = step(env, action);
        
        % Update max Pd
        maxPd = max(maxPd, -reward); % Reward is negative Pd
        
        % Update observation
        obs = nextObs;
    end
end

