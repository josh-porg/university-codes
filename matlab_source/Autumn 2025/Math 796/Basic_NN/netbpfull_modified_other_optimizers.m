function netbp_full
%NETBP_FULL
%   Extended version of netbp, with more graphics
%
%   Set up data for neural net test
%   Use backpropagation to train 
%   Visualize results
%
% C F Higham and D J Higham, Aug 2017
%
%%%%%%% DATA %%%%%%%%%%%
% xcoords, ycoords, targets
x1 = [0.1,0.3,0.1,0.6,0.4,0.6,0.5,0.9,0.4,0.7];
x2 = [0.1,0.4,0.5,0.9,0.2,0.3,0.6,0.2,0.4,0.6];
y = [ones(1,5) zeros(1,5); zeros(1,5) ones(1,5)];
% new definition of y that allows easier modification
y_pre_prelim = [0,1,1,1,1];
y_prelim = [y_pre_prelim, ~y_pre_prelim];
y = [y_prelim; ~y_prelim]

figure(1)
clf
a1 = subplot(1,1,1);
plot(x1(1:5),x2(1:5),'ro','MarkerSize',12,'LineWidth',4)
hold on
plot(x1(6:10),x2(6:10),'bx','MarkerSize',12,'LineWidth',4)
a1.XTick = [0 1];
a1.YTick = [0 1];
a1.FontWeight = 'Bold';
a1.FontSize = 16;
xlim([0,1])
ylim([0,1])

%print -dpng pic_xy.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize weights and biases 
%rng(5000);
W2 = 0.5*randn(2,2);
W3 = 0.5*randn(3,2);
W4 = 0.5*randn(2,3);
b2 = 0.5*randn(2,1);
b3 = 0.5*randn(3,1);
b4 = 0.5*randn(2,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% store sizes of weight and biases matrixies for reshaping to/from
% vetor/matrix formatts
sizes = {[size(W2,1),size(W2,2)], [size(W3,1),size(W3,2)], [size(W4,1),size(W4,2)], [size(b2,1),size(b2,2)], [size(b3,1),size(b3,2)], [size(b4,1),size(b4,2)]};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward and Back propagate 
% Pick a training point at random
eta = 0.05;
new_eta = eta;
Niter = 1e6;
savecost = zeros(Niter,1);

design_point = [W2(:); W3(:); W4(:); b2; b3; b4];

%dp = optimvar("dp", length(design_point));
% prob = optimproblem("Objective", vector_cost(dp, sizes));
% Set initial point
%x0.dp = design_point;  % Fixed: need to specify field name

% Use direct GA call (this approach works with nested functions)
objective_function = @(x) vector_cost(x, sizes);
lb = -1e5 * ones(length(design_point), 1);  % Lower bounds
ub = 1e5 * ones(length(design_point), 1);   % Upper bounds

% Configure GA options
% ga_options = optimoptions('ga', ...
%     'Display', 'iter', ...
%     'MaxGenerations', 1000, ...
%     'PopulationSize', 5000, ...
%     'InitialPopulationMatrix', design_point');  % Use current weights as starting point

% these work to keep it from converging to local mins but dont stop running
% 
% ga_options = optimoptions('ga', ...
%     'Display', 'iter', ...
%     'MaxGenerations', 5000, ...  % Increased max generations
%     'PopulationSize', 1000, ...   % Larger population for better exploration
%     'InitialPopulationMatrix', design_point', ...  % Use current weights as starting point
%     'MaxStallGenerations', 200, ...  % Allow more stall generations before stopping
%     'FunctionTolerance', 1e-8, ...   % Tighter function tolerance
%     'ConstraintTolerance', 1e-6, ...
%     'CreationFcn', @gacreationuniform, ...  % Ensure good initial diversity
%     'CrossoverFraction', 0.8, ...   % Higher crossover rate
%     'MutationFcn', {@mutationadaptfeasible, 0.1}, ...  % Adaptive mutation
%     'SelectionFcn', @selectiontournament, ...
%     'EliteCount', 50, ...  % Keep best individuals
%     'UseParallel', false);  % Set to true if you have Parallel Computing Toolbox


ga_options = optimoptions('ga', ...
    'Display', 'iter', ...
    'MaxGenerations', 5000, ...  % Increased max generations
    'PopulationSize', 1000, ...   % Larger population for better exploration
    'InitialPopulationMatrix', design_point', ...  % Use current weights as starting point
    'MaxStallGenerations', 200, ...  % Allow more stall generations before stopping
    'FunctionTolerance', 1e-8, ...   % Tighter function tolerance
    'ConstraintTolerance', 1e-6, ...
    'CreationFcn', @gacreationuniform, ...  % Ensure good initial diversity
    'CrossoverFraction', 0.8, ...   % Higher crossover rate
    'MutationFcn', {@mutationadaptfeasible, 0.1}, ...  % Adaptive mutation
    'SelectionFcn', @selectiontournament, ...
    'EliteCount', 50, ...  % Keep best individuals
    'UseParallel', false);  % Set to true if you have Parallel Computing Toolbox

% run the solver
[solution_vector, fval, eflag, output] = ga(objective_function, length(design_point), ...
    [], [], [], [], lb, ub, [], ga_options);

% Turn the vector back into individual weight and bias matrices
strtIdx = 1;
for j = 1:length(sizes)
    nElem = prod(sizes{j}); % Fixed variable name consistency
    % Extract segment from linear vector
    seg = solution_vector(strtIdx : strtIdx + nElem - 1);  % Fixed: use solution_vector
    % Reshape into original matrix size
    ogMat{j} = reshape(seg, sizes{j}(1), sizes{j}(2));
    % Update index for next matrix
    strtIdx = strtIdx + nElem;
end

% Assign them to the weights and bias matrices
W2 = ogMat{1};
W3 = ogMat{2};
W4 = ogMat{3};
b2 = ogMat{4};
b3 = ogMat{5};
b4 = ogMat{6};

% Display final cost
final_cost = cost(W2,W3,W4,b2,b3,b4);
fprintf('GA optimization completed. Final cost: %g\n', final_cost);

% for counter = 1:Niter
%     k = randi(10);
%     x = [x1(k); x2(k)];
% 
% 
%     % Forward pass
%     a2 = activate(x,W2,b2);
%     a3 = activate(a2,W3,b3);
%     a4 = activate(a3,W4,b4);
% 
%     % Backward pass
%     delta4 = a4.*(1-a4).*(a4-y(:,k));
%     delta3 = a3.*(1-a3).*(W4'*delta4);
%     delta2 = a2.*(1-a2).*(W3'*delta3);
% 
%     % for clarity lets compute the gradient
%     d_W2 = delta2*x';
%     d_W3 = delta3*a2';
%     d_W4 = delta4*a3';
%     d_b2 = delta2;
%     d_b3 = delta3;
%     d_b4 = delta4;
% 
% 
%     % for this to work with my inexact step size code I need the descint
%     % function (because this is unconstrained that is simply the cost
%     % function) the format must be takes a vector x and produces a scalar
%     % cost. thus i must write a wrapper to make a design point into a
%     % vector and a wrapper to take said vecotr and pass it into cost
%     design_point = [W2(:); W3(:); W4(:); b2; b3; b4];
%     descent_direction = -[d_W2(:); d_W3(:); d_W4(:); d_b2; d_b3; d_b4];
% 
%     % compute the step size (learning rate) using inexact line search
%     gamma = .1;
%     mu = .5;
%     % TODO: modify the code to only compute the exat step size every n
%     % itterations? or after n itterations?
%     % if(counter > 5e5)
%     %     new_eta = inexactStepSize(design_point, descent_direction, @(x) vector_cost(x, sizes), gamma, mu);
%     % end
%         %new_eta = 0.05;
% 
%     if new_eta < .05
%         eta = 1e-5;
%     else
%         %fprintf("new eta: %g", eta)
%         eta = new_eta;
%     end
%     % it currently seems as if the descent direction is not always an
%     % actual descent direction
% 
%     % Gradient step
%     W2 = W2 - eta*delta2*x';
%     W3 = W3 - eta*delta3*a2';
%     W4 = W4 - eta*delta4*a3';
%     b2 = b2 - eta*delta2;
%     b3 = b3 - eta*delta3;
%     b4 = b4 - eta*delta4;
% 
%     % Monitor progress
%     newcost = cost(W2,W3,W4,b2,b3,b4);   
%     if mod(counter,100000) == 0
%         % display cost to screen
%         fprintf("iteratiion %06d; cost: %g; eta: %g\n", counter, newcost, eta);
%     end
%     savecost(counter) = newcost;
% end

figure(2)
clf
semilogy([1:1e4:Niter],savecost(1:1e4:Niter),'b-','LineWidth',2)
xlabel('Iteration Number')
ylabel('Value of cost function')
set(gca,'FontWeight','Bold','FontSize',18)
print -dpng pic_cost.png

%%%%%%%%%%% Display shaded and unshaded regions 
N = 500;
Dx = 1/N;
Dy = 1/N;
xvals = [0:Dx:1];
yvals = [0:Dy:1];
for k1 = 1:N+1
    xk = xvals(k1);
    for k2 = 1:N+1
        yk = yvals(k2);
        xy = [xk;yk];
        a2 = activate(xy,W2,b2);
        a3 = activate(a2,W3,b3);
        a4 = activate(a3,W4,b4);
        Aval(k2,k1) = a4(1);
        Bval(k2,k1) = a4(2);
     end
end
[X,Y] = meshgrid(xvals,yvals);

figure(3)
clf
a2 = subplot(1,1,1);
Mval = Aval>Bval;
contourf(X,Y,Mval,[0.5 0.5])
hold on
colormap([1 1 1; 0.8 0.8 0.8])
% Automatic plotting based on binary classification
% Assumes y is a 2xN matrix where each column is one-hot encoded
% Find indices for each class
class_a_indices = find(y(1,:) == 1);  % Where first row is 1
class_b_indices = find(y(2,:) == 1);  % Where second row is 1
% Plot class A (category 1)
plot(x1(class_a_indices), x2(class_a_indices), 'ro', 'MarkerSize', 12, 'LineWidth', 4)
hold on
% Plot class B (category 2)  
plot(x1(class_b_indices), x2(class_b_indices), 'bx', 'MarkerSize', 12, 'LineWidth', 4)
hold off
%legend('Category A', 'Category B', 'Location', 'best') % Optional: Add legend
a2.XTick = [0 1];
a2.YTick = [0 1];
a2.FontWeight = 'Bold';
a2.FontSize = 16;
xlim([0,1])
ylim([0,1])

print -dpng pic_bdy_bp.png

    function costval = cost(W2,W3,W4,b2,b3,b4)

        costvec = zeros(10,1);
        for i = 1:10 % Loop through ALL 10 datapoints
            x =[x1(i);x2(i)];
            a2 = activate(x,W2,b2);
            a3 = activate(a2,W3,b3);
            a4 = activate(a3,W4,b4);
            costvec(i) = norm(y(:,i) - a4,2); % L2 norm of error for point i
        end
        costval = norm(costvec,2)^2;  % Sum of squares of all individual errors
    end % of nested function

    % a wrapper for the cost function to allow a vecotr of parameters as an
    % input.
    function costval = vector_cost(design_point, sizes)
        % turn the vecotr into the individual weight and bias matricies
        startIdx = 1;
        for k = 1:length(sizes)
            nElements = prod(sizes{k}); % number of elements in kth matrix
            % Extract segment from linear vector
            segment = design_point(startIdx : startIdx + nElements - 1);
            % Reshape into original matrix size
            originalMatrix{k} = reshape(segment, sizes{k}(1), sizes{k}(2));
            % Update index for next matrix
            startIdx = startIdx + nElements;
        end

        % assign them to the wieghts and bias matricies
        W2 = originalMatrix{1};
        W3 = originalMatrix{2};
        W4 = originalMatrix{3};
        b2 = originalMatrix{4};
        b3 = originalMatrix{5};
        b4 = originalMatrix{6};

        %call the cost function
        costval = cost(W2,W3,W4,b2,b3,b4);
    end


end
