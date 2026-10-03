%% deal with the training data
% Load the saved simulation data
load('simulation_data.mat'); % Loads all_n_history, all_Active_history, all_Pd_history, all_Phi_history

% Combine features into input matrix (add other features as needed)
X = [all_n_history, all_Active_history, all_Phi_history];%, all_Phi_history];
y = all_Pd_history;

% Remove any rows with NaN values (if they exist)
valid_samples = ~any(isnan(X), 2) & ~isnan(y);
X = X(valid_samples, :);
y = y(valid_samples);

% Split data into training (60%), validation (20%), and test (20%) sets
rng(42); % For reproducibility
n_samples = size(X, 1);
shuffled_indices = randperm(n_samples);

train_end = floor(0.6 * n_samples);
val_end = floor(0.8 * n_samples);

train_indices = shuffled_indices(1:train_end);
val_indices = shuffled_indices(train_end+1:val_end);
test_indices = shuffled_indices(val_end+1:end);

X_train = X(train_indices, :);
y_train = y(train_indices);

X_val = X(val_indices, :);
y_val = y(val_indices);

X_test = X(test_indices, :);
y_test = y(test_indices);

%% configure the network

%shape of network
layerSizes = [8, 8, 8,8,8,8]; % Three hidden layers with 3 neurons each
nOutputs = 1; % Number of output neurons

% training options
options.learningRate = 0.01;
options.maxEpochs = 500;
options.minError = 1e-5;
options.batchSize = 32; % Mini-batch size
options.momentum = 0.9;
options.weightDecay = 0.001; % L2 regularization

% configure training data
load n

%Train the network
[trainedWeights, trainingError] = trainNN(X_train, y_train, layerSizes, nOutputs, options);


% Evaluate on validation set
val_predictions = zeros(size(y_val));
for i = 1:size(X_val, 1)
    val_predictions(i) = mySimpleNN(trainedWeights, X_val(i,:)', layerSizes, nOutputs);
end

% Calculate validation RMSE
val_rmse = sqrt(mean((val_predictions - y_val).^2));
fprintf('Validation RMSE: %.4f\n', val_rmse);

% Final evaluation on test set (only do this once!)
test_predictions = zeros(size(y_test));
for i = 1:size(X_test, 1)
    test_predictions(i) = mySimpleNN(trainedWeights, X_test(i,:)', layerSizes, nOutputs);
end

% Calculate test RMSE
test_rmse = sqrt(mean((test_predictions - y_test).^2));
fprintf('Test RMSE: %.4f\n', test_rmse);

% Optional: Plot predictions vs actual
figure;
subplot(2,1,1);
plot(y_val, val_predictions, 'o');
xlabel('Actual Pd');
ylabel('Predicted Pd');
title('Validation Set Performance');
grid on;

subplot(2,1,2);
plot(y_test, test_predictions, 'o');
xlabel('Actual Pd');
ylabel('Predicted Pd');
title('Test Set Performance');
grid on;

function [trainedWeights, trainingError] = trainNN(trainingData, trainingTargets, layerSizes, nOutputs, options)
    % TRAINNN Trains the neural network using backpropagation
    %
    % Inputs:
    %   trainingData - NxD matrix where N is number of samples, D is input dimension
    %   trainingTargets - NxM matrix where M is number of outputs
    %   layerSizes - Vector specifying number of neurons in each hidden layer
    %   nOutputs - Number of output neurons
    %   options - Struct with training parameters:
    %       .learningRate - Initial learning rate
    %       .maxEpochs - Maximum number of training epochs
    %       .minError - Minimum error threshold for early stopping
    %       .batchSize - Size of mini-batch (0 for full batch)
    %       .momentum - Momentum coefficient (0 for no momentum)
    %       .weightDecay - L2 regularization coefficient
    %
    % Outputs:
    %   trainedWeights - Final weight vector after training
    %   trainingError - Vector of training errors across epochs

    % Initialize weights
    initialWeights = initializeWeights(layerSizes, nOutputs, size(trainingData, 2));
    
    % Set default options if not provided
    if nargin < 5
        options = struct();
    end
    if ~isfield(options, 'learningRate'), options.learningRate = 0.01; end
    if ~isfield(options, 'maxEpochs'), options.maxEpochs = 1000; end
    if ~isfield(options, 'minError'), options.minError = 1e-5; end
    if ~isfield(options, 'batchSize'), options.batchSize = 0; end
    if ~isfield(options, 'momentum'), options.momentum = 0.9; end
    if ~isfield(options, 'weightDecay'), options.weightDecay = 0.001; end
    
    % Initialize variables for training
    currentWeights = initialWeights;
    velocity = zeros(size(currentWeights)); % For momentum
    trainingError = zeros(options.maxEpochs, 1);
    nSamples = size(trainingData, 1);
    
    % Determine batch configuration
    if options.batchSize <= 0 || options.batchSize >= nSamples
        nBatches = 1;
        batchSize = nSamples;
    else
        nBatches = ceil(nSamples / options.batchSize);
        batchSize = options.batchSize;
    end
    
    % Training loop
    for epoch = 1:options.maxEpochs
        % Shuffle data for each epoch
        shuffledIndices = randperm(nSamples);
        
        epochError = 0;
        
        for batch = 1:nBatches
            % Get current batch
            batchStart = (batch-1)*batchSize + 1;
            batchEnd = min(batch*batchSize, nSamples);
            batchIndices = shuffledIndices(batchStart:batchEnd);
            
            batchData = trainingData(batchIndices, :)';
            batchTargets = trainingTargets(batchIndices, :)';
            
            % Forward pass
            [output, layerOutputs] = forwardPass(currentWeights, batchData, layerSizes, nOutputs);
            
            % Compute error
            error = output - batchTargets;
            batchError = 0.5 * sum(error(:).^2) / numel(batchIndices);
            epochError = epochError + batchError;
            
            % Backward pass (compute gradient)
            grad = backpropagation(currentWeights, batchData, batchTargets, layerSizes, nOutputs, layerOutputs);
            
            % Add regularization gradient
            grad = grad + options.weightDecay * currentWeights;
            
            % Update weights with momentum
            velocity = options.momentum * velocity - options.learningRate * grad;
            currentWeights = currentWeights + velocity;
        end
        
        % Store training error
        trainingError(epoch) = epochError / nBatches;
        
        % Check for early stopping
        if trainingError(epoch) < options.minError
            trainingError = trainingError(1:epoch);
            fprintf('Training converged after %d epochs\n', epoch);
            break;
        end
        
        % Optional: learning rate decay
        % options.learningRate = options.learningRate * 0.995;
        
        % Display progress
        if mod(epoch, 10) == 0
            fprintf('Epoch %d, Error: %f\n', epoch, trainingError(epoch));
        end
    end
    
    trainedWeights = currentWeights;
end

function [output, layerOutputs] = forwardPass(w, inp, layerSizes, nOutputs)
    % Extended forward pass that also returns all layer outputs
    inputSize = size(inp, 1);
    nLayers = length(layerSizes);
    idx = 0;
    layerInput = inp;
    layerOutputs = cell(nLayers + 1, 1); % +1 for output layer
    layerOutputs{1} = inp; % Store input as first "layer"
    
    % Hidden layers
    for i = 1:nLayers - 1
        currentLayerSize = layerSizes(i);
        
        % Extract weights and bias
        nWeights = currentLayerSize * inputSize;
        weights = reshape(w(idx + 1 : idx + nWeights), currentLayerSize, inputSize);
        idx = idx + nWeights;
        bias = reshape(w(idx + 1 : idx + currentLayerSize), currentLayerSize, 1);
        idx = idx + currentLayerSize;
        
        % Forward pass with tanh activation
        layerInput = tanh(weights * layerInput + bias);
        layerOutputs{i+1} = layerInput; % Store after activation
        
        inputSize = currentLayerSize;
    end
    
    % Output layer (linear activation)
    W_final = reshape(w(idx + 1 : idx + layerSizes(end)*nOutputs), nOutputs, layerSizes(end));
    idx = idx + layerSizes(end)*nOutputs;
    b_final = reshape(w(idx + 1 : idx + nOutputs), nOutputs, 1);
    
    output = W_final * layerInput + b_final;
    layerOutputs{end} = output;
end

function grad = backpropagation(w, inp, target, layerSizes, nOutputs, layerOutputs)
    % BACKPROPAGATION Computes the gradient using backpropagation
    nLayers = length(layerSizes);
    
    % Initialize gradient
    grad = zeros(size(w));
    
    % Extract all weights and biases from the weight vector
    [weights, biases] = unpackWeights(w, layerSizes, nOutputs, size(inp, 1));
    
    % Output layer gradient
    output = layerOutputs{end};
    error = output - target;
    delta = error; % Linear activation derivative is 1
    
    % Gradient for output layer weights
    hiddenOutput = layerOutputs{end-1}; % Output of last hidden layer
    dW_final = delta * hiddenOutput';
    db_final = sum(delta, 2);
    
    % Find position in weight vector for output layer
    idx = numel(w) - (layerSizes(end)*nOutputs + nOutputs);
    grad(idx+1:idx+numel(dW_final)) = dW_final(:);
    grad(idx+numel(dW_final)+1:idx+numel(dW_final)+numel(db_final)) = db_final(:);
    
    % Backpropagate error
    if nLayers > 1
        delta = (weights{end}' * delta) .* (1 - hiddenOutput.^2); % tanh derivative
    end
    
    % Backpropagate through hidden layers
    for l = nLayers-1:-1:1
        % Get layer input (output of previous layer)
        if l == 1
            layerInput = layerOutputs{1}; % Original input
        else
            layerInput = layerOutputs{l};
        end
        
        % Compute gradients
        dW = delta * layerInput';
        db = sum(delta, 2);
        
        % Find position in weight vector for this layer
        if l == 1
            inputSize = size(inp, 1);
        else
            inputSize = layerSizes(l-1);
        end
        nWeights = layerSizes(l) * inputSize;
        nBiases = layerSizes(l);
        
        if l == 1
            idx = 0;
        else
            % Calculate starting index for this layer's weights
            idx = 0;
            for i = 1:l-1
                if i == 1
                    prevSize = size(inp, 1);
                else
                    prevSize = layerSizes(i-1);
                end
                idx = idx + layerSizes(i) * prevSize + layerSizes(i);
            end
        end
        
        % Store gradients
        grad(idx+1:idx+nWeights) = dW(:);
        grad(idx+nWeights+1:idx+nWeights+nBiases) = db(:);
        
        % Propagate error backward if not first layer
        if l > 1
            delta = (weights{l}' * delta) .* (1 - layerInput.^2); % tanh derivative
        end
    end
end

function [weights, biases] = unpackWeights(w, layerSizes, nOutputs, inputSize)
    % UNPACKWEIGHTS Extracts all weights and biases from the weight vector
    idx = 0;
    nLayers = length(layerSizes);
    weights = cell(nLayers, 1);
    biases = cell(nLayers, 1);
    
    % Hidden layers
    for i = 1:nLayers - 1
        currentLayerSize = layerSizes(i);
        
        % Extract weights
        nWeights = currentLayerSize * inputSize;
        weights{i} = reshape(w(idx + 1 : idx + nWeights), currentLayerSize, inputSize);
        idx = idx + nWeights;
        
        % Extract bias
        biases{i} = reshape(w(idx + 1 : idx + currentLayerSize), currentLayerSize, 1);
        idx = idx + currentLayerSize;
        
        inputSize = currentLayerSize;
    end
    
    % Output layer
    weights{end} = reshape(w(idx + 1 : idx + layerSizes(end)*nOutputs), nOutputs, layerSizes(end));
    idx = idx + layerSizes(end)*nOutputs;
    biases{end} = reshape(w(idx + 1 : idx + nOutputs), nOutputs, 1);
end

function w = initializeWeights(layerSizes, nOutputs, inputSize)
    % INITIALIZEWEIGHTS Initializes weights with Xavier/Glorot initialization
    totalWeights = 0;
    nLayers = length(layerSizes);
    
    % Calculate total number of weights needed
    % Input to first hidden layer
    totalWeights = totalWeights + layerSizes(1) * inputSize + layerSizes(1);
    
    % Hidden to hidden layers
    for i = 2:nLayers
        totalWeights = totalWeights + layerSizes(i) * layerSizes(i-1) + layerSizes(i);
    end
    
    % Last hidden to output layer
    totalWeights = totalWeights + nOutputs * layerSizes(end) + nOutputs;
    
    % Initialize all weights
    w = zeros(totalWeights, 1);
    idx = 0;
    
    % Initialize input to first hidden layer
    range = sqrt(6 / (inputSize + layerSizes(1)));
    nWeights = layerSizes(1) * inputSize;
    w(idx+1:idx+nWeights) = -range + 2*range * rand(nWeights, 1);
    idx = idx + nWeights;
    
    % Biases for first hidden layer
    w(idx+1:idx+layerSizes(1)) = zeros(layerSizes(1), 1);
    idx = idx + layerSizes(1);
    
    % Initialize hidden to hidden layers
    for i = 2:nLayers
        range = sqrt(6 / (layerSizes(i-1) + layerSizes(i)));
        nWeights = layerSizes(i) * layerSizes(i-1);
        w(idx+1:idx+nWeights) = -range + 2*range * rand(nWeights, 1);
        idx = idx + nWeights;
        
        % Biases
        w(idx+1:idx+layerSizes(i)) = zeros(layerSizes(i), 1);
        idx = idx + layerSizes(i);
    end
    
    % Initialize last hidden to output layer
    range = sqrt(6 / (layerSizes(end) + nOutputs));
    nWeights = nOutputs * layerSizes(end);
    w(idx+1:idx+nWeights) = -range + 2*range * rand(nWeights, 1);
    idx = idx + nWeights;
    
    % Output biases
    w(idx+1:idx+nOutputs) = zeros(nOutputs, 1);
end