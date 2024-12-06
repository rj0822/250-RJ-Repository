% CSE 250A HW-6
% Viterbi Algorithm for Hidden Markov Model (HMM).

clear; clc; close all;

% Initialization. ---------------------------------------------------------
% Load the data files.
Emission = load('emissionMatrix.txt'); % Load emission matrix.
Transition = load('transitionMatrix.txt'); % Load transition matrix.
Observation = load('observations.txt'); % Load observations.
InitialState = load('initialStateDistribution.txt'); % Load initial state.

% Define parameters.
n = size(Transition, 1); % Number of states.
T = size(Observation,2); % Number of observations (number of time steps).

% Initialize arrays for probability and backtracing indices.
Probability = zeros(n, T); % Log probability for each state.
Indices = zeros(n, T); % Backtracking indices.

% Viterbi algorithm. ------------------------------------------------------
% Calculate the first column.
for i = 1:n
    Probability(i, 1) = log(InitialState(i)) + log(Emission(i, Observation(1)+1));
end

% Process each time step.
for t = 2:T % Loop over each time step.
    for j = 1:n % Loop over each state.
        CurrentProbability_max = -inf; % Initialize the current log probability as -inf.
        CurrentIndex = -1; % Initialize the current backtracking index as -1.

        for i = 1:n
            CurrentProbability = Probability(i, t-1) + log(Transition(i,j));
            if CurrentProbability > CurrentProbability_max
                CurrentProbability_max = CurrentProbability;
                CurrentIndex = i;
            end
        end

        Indices(j, t) = CurrentIndex;
        Probability(j, t) = CurrentProbability_max + log(Emission(j, Observation(t)+1));
    end
end

% Backtracking to find the most likely sequence of hidden states. ---------
Sequence = zeros(1, T); % Initialize the most likely sequence.
[~, Sequence(T)] = max(Probability(:, T)); % Find the maximum probability in the last column.

for i = T-1:-1:1
    Sequence(i) = Indices(Sequence(i+1), i+1);
end

% Plot and answer check. --------------------------------------------------
% Plot the most likely sequence of hidden states vs time step.
set(0, 'DefaultAxesFontSize', 15);
set(0, 'DefaultTextFontSize', 15);

set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');

figure('Position', [0, 0, 800, 600]);
plot(1:T, Sequence, 'b', 'LineWidth', 2);
grid on;
box on;
xlim([0 450000]);
ylim([0 30]);
xticks(0:50000:450000);
yticks(0:5:30);
xlabel('Time Step');
ylabel('Hidden State Index');
title('Most Likely Sequence of Hidden States');

% Answer check.
OutputIndices = Sequence(1); % Initialize the output with the first hidden state.

% If the state is not repeated, then add it to the output.
for i = 2:T
    if Sequence(i) ~= Sequence(i-1)
        OutputIndices = [OutputIndices, Sequence(i)];
    end
end

% Map the indices to the output string.
Letter = 'abcdefghijklmnopqrstuvwxyz ';
String = Letter(OutputIndices);
disp(String);
