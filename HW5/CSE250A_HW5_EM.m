% CSE 250A HW-5
% Expectation-Maximization (EM) Algorithm for Noisy-OR.

clear; clc; close all;

% Initialization. ---------------------------------------------------------
% Load the data files
x = load('spectX-1.txt'); % Load input data
y = load('spectY-1.txt'); % Load output data
NumOfSamples = size(x,1); % Number of samples/tests = Row number of the input
NumOfVariables = size(x,2); % Number of variables = Column number of the input

% Define the maximum number of iterations
MaxIterations = 256;

% Initialize all the p_i with 0.05.
p = 0.05 * ones(NumOfVariables,1)'; 

% Initialize arrays for the count of mistakes and log-likelihood.
Mistakes = zeros(MaxIterations+1,1);
LogLikelihood = zeros(MaxIterations+1,1);

% EM iterations. ----------------------------------------------------------
% Notations: i --> Number of iterations
%            t --> Number of samples/tests
%            n --> Number of input variables
for i = 0:MaxIterations
    % Compute count of mistakes and log-likelihood.
    for t = 1:NumOfSamples
        PY_X = 1 - prod(1 - p.*x(t,:)); % Calculate P(Y=1|X).
        LogLikelihood(i+1) = LogLikelihood(i+1) + y(t)*log(PY_X) + (1-y(t))*log(1 - PY_X); % Calculate log-likelihood.
        
        % Mistake count.
        if (y(t) == 0 && PY_X >= 0.5) || (y(t) == 1 && PY_X < 0.5)
            Mistakes(i+1) = Mistakes(i+1) + 1;
        end
    end
    
    % Normalize the log-likelihood.
    LogLikelihood(i+1) = LogLikelihood(i+1) / NumOfSamples;

    % Calculate P(Z=1|X,Y).
    PZ_XY = zeros(NumOfSamples, NumOfVariables);

    for t = 1:NumOfSamples
        PY_X = 1 - prod(1 - p.*x(t,:)); % Compute P(Y=1|X).
        for n = 1:NumOfVariables
            if x(t,n) == 1
                PZ_XY(t,n) = (y(t)*p(n)) / PY_X;
            end
        end
    end
    
    % Update parameter p_i.
    for n = 1:NumOfVariables
        Ti = sum(x(:,n) == 1); % Count samples with X=1.
        p(n) = sum(PZ_XY(x(:,n) == 1, n)) / Ti;
    end
end

% Display the results. ----------------------------------------------------
Results = array2table([(0:1:MaxIterations)' Mistakes LogLikelihood], 'VariableNames', {'Iteration', 'Number of Mistakes', 'Log-likelihood'});
DispIterations = [0 1 2 4 8 16 32 64 128 256];
disp(Results(DispIterations+1,:));
