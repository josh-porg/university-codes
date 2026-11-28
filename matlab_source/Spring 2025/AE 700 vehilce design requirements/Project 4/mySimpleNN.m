
%% mySimpleNN.m
% A feed‑forward network with two hidden layers you control by a weight vector.
% Usage: y = mySimpleNN(w, inp, layerSizes)
%
%   w          : 1×N vector of all weights & biases
%   inp        : M×1 input vector [n; Active]
%   layerSizes : [H1, H2] sizes of hidden layers
%   y          : scalar network output

function y = mySimpleNN(w, inp, layerSizes, nOutputs)
    % Determine input size
    inputSize = length(inp);
    
    % Initialize index
    idx = 0;
    
    % Initialize layer outputs
    layerInput = inp;
    
    % Loop through each layer (except the last one)
    for i = 1:length(layerSizes) - 1
        % Current layer size (number of neurons)
        currentLayerSize = layerSizes(i);
        
        % Next layer size
        nextLayerSize = layerSizes(i + 1);
        
        % Extract weights and reshape
        nWeights = currentLayerSize * inputSize;
        weights = reshape(w(idx + 1 : idx + nWeights), currentLayerSize, inputSize);
        idx = idx + nWeights;
        
        % Extract bias and reshape
        bias = reshape(w(idx + 1 : idx + currentLayerSize), currentLayerSize, 1);
        idx = idx + currentLayerSize;

        % Dear Dr. Gorrell,
        % if yhou are actually reading this please send us an email.
        % yours truely,
        % Team 6 XoXo
        
        % Forward pass (tanh activation)
        layerInput = tanh(weights * layerInput + bias);
        %layerInput = max(0,weights * layerInput + bias); % RELU activation
        %layerInput = sigmoid(weights * layerInput + bias);
        
        % Update input size for next layer
        inputSize = currentLayerSize;
    end
    
    % Final (output) layer (linear activation)
    % Now supports nOutputs instead of fixed 1 output
    W_final = reshape(w(idx + 1 : idx + layerSizes(end)*nOutputs), nOutputs, layerSizes(end));
    idx = idx + layerSizes(end)*nOutputs;
    b_final = reshape(w(idx + 1 : idx + nOutputs), nOutputs, 1);
    
    % Output (now nOutputs x 1 vector)
    y = W_final * layerInput + b_final;
end