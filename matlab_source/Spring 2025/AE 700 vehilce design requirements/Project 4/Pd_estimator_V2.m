%% Load and prepare data
load('simulation_data.mat');

% Combine features and target
X = [all_n_history, all_Active_history, all_Phi_history];
y = all_Pd_history;

% Remove NaN values
valid_samples = ~any(isnan(X), 2) & ~isnan(y);
X = X(valid_samples, :);
y = y(valid_samples);

% Split data (60% train, 20% val, 20% test)
rng(42); % For reproducibility
n = size(X, 1);
idx = randperm(n);
train_idx = idx(1:floor(0.6*n));
val_idx = idx(floor(0.6*n)+1:floor(0.8*n));
test_idx = idx(floor(0.8*n)+1:end);

X_train = X(train_idx, :)';
y_train = y(train_idx)';
X_val = X(val_idx, :)';
y_val = y(val_idx)';
X_test = X(test_idx, :)';
y_test = y(test_idx)';

%% Create and train network
net = fitnet([16 16 16]); % 3 hidden layers with 16 neurons each

% Configure training parameters
net.divideFcn = 'divideind'; % Manual data division
net.divideParam.trainInd = 1:length(train_idx);
net.divideParam.valInd = length(train_idx)+1:length(train_idx)+length(val_idx);
net.divideParam.testInd = [];

net.trainFcn = 'trainlm'; % Levenberg-Marquardt (fast for small/medium nets)
% Alternatives: 'trainbr' (Bayesian regularization), 'trainscg' (scaled conjugate gradient)

net.performFcn = 'mse'; % Mean squared error
net.trainParam.epochs = 1000;
net.trainParam.max_fail = 20; % Early stopping if val error increases for 20 epochs
net.trainParam.min_grad = 1e-7;
net.trainParam.showWindow = true; % Show training GUI

% Train the network
[net, tr] = train(net, [X_train X_val], [y_train y_val]);

%% Evaluate performance
% Predictions
y_train_pred = net(X_train);
y_val_pred = net(X_val);
y_test_pred = net(X_test);

% Calculate RMSE
train_rmse = sqrt(mean((y_train_pred - y_train).^2));
val_rmse = sqrt(mean((y_val_pred - y_val).^2));
test_rmse = sqrt(mean((y_test_pred - y_test).^2));

fprintf('Training RMSE: %.4f\n', train_rmse);
fprintf('Validation RMSE: %.4f\n', val_rmse);
fprintf('Test RMSE: %.4f\n', test_rmse);

%% Plot results
figure;
subplot(3,1,1);
plot(y_train, y_train_pred, 'o');
title('Training Set Performance');
xlabel('Actual Pd');
ylabel('Predicted Pd');
grid on;

subplot(3,1,2);
plot(y_val, y_val_pred, 'o');
title('Validation Set Performance');
xlabel('Actual Pd');
ylabel('Predicted Pd');
grid on;

subplot(3,1,3);
plot(y_test, y_test_pred, 'o');
title('Test Set Performance');
xlabel('Actual Pd');
ylabel('Predicted Pd');
grid on;

%% Save the trained network
save('trained_network.mat', 'net');

% Data preprocessing
net.input.processFcns = {'removeconstantrows', 'mapminmax'};
net.output.processFcns = {'removeconstantrows', 'mapminmax'};

% Regularization
net.performParam.regularization = 0.1; % L2 regularization

% Bayesian Regularization (often works better)
net.trainFcn = 'trainbr'; 
net.trainParam.epochs = 500;

