% CSE 250A HW7
% Movie Recommendation System Based on Expectation-Maximization (EM) Algorithm.

clear; clc; close all;

% Load the data file. -----------------------------------------------------
StudentID = readlines('hw7_ids.txt'); % Student ID.
StudentID(end) = [];
MovieList = readlines('hw7_movies.txt'); % Movie list.
MovieList(end) = [];

% Initial probabilities.
InitialProbR = load('hw7_probR_init.txt');
InitialProbZ = load('hw7_probZ_init.txt');

% Movie rating inputs.
% Raw ratings: 1 --> good, 0 --> bad, ? --> unseen. The initial rating is string type.
% Updated ratings: 1 --> good, 0 --> bad, NaN --> unseen. The updated rating is double type.
RawRatings = readlines('hw7_ratings.txt');
RawRatings(end) = [];

Ratings = zeros(size(RawRatings,1), size(MovieList,1));

for i = 1:size(Ratings,1)
    Ratings(i,:) = strsplit(RawRatings(i,:));
end

% Sanity check. -----------------------------------------------------------
% Calculate the mean ratings of each movie.
MeanRatings = zeros(size(MovieList,1),1);

for i = 1:size(MovieList,1)
    MeanRatings(i) = mean(Ratings(~isnan(Ratings(:,i)),i));
end

% Print out the movie titles from least to most popular along with mean ratings.
MovieMeanRatings = array2table(MeanRatings, 'VariableNames', {'Mean Rating'});
MovieMeanRatings = addvars(MovieMeanRatings, MovieList, 'Before', 1, 'NewVariableNames', 'Movie');
MovieMeanRatingsSorted = sortrows(MovieMeanRatings, 2);
disp('(a) Mean popularity ratings of all the movies:')
disp('The mean popularity ratings does not match my own preference from the observation.');
disp(MovieMeanRatingsSorted);

% EM algorithm implementation. --------------------------------------------
% Initialization.
k = 4; % Number of movie-goer types.
N = size(MovieList,1); % Number of movies.
T = size(Ratings,1); % Number of students.
MaxIter = 256; % Maximum number of itertions.
Posteriors = zeros(k, T);
ProbZ = InitialProbZ;
ProbR = InitialProbR;
L = zeros(MaxIter+1, 1);

print_itr = 1;

for iter = 0:MaxIter
    % E-step.
    Posteriors = zeros(k, T);

    for t = 1:T % Loop over each student.
        Recommended = find(Ratings(t, :) == 1); % Movies recommended by student t
        NotRecommended = find(Ratings(t, :) == 0); % Movies not recommended by student t
        
        % Calculate log-likelihood.
        P = 0;

        for i = 1:k
            P = P + ProbZ(i) * prod(ProbR(Recommended, i)) * prod(1 - ProbR(NotRecommended, i));
        end

        L(iter+1) = L(iter+1) + log(P);
        
        % Calculate denominator for E-Step.
        Estep_denominator = 0;

        for i = 1:k
            Estep_denominator = Estep_denominator + ProbZ(i) * prod(ProbR(Recommended, i)) * prod(1 - ProbR(NotRecommended, i));
        end
        
        % Update posterior probabilities.
        for i = 1:k
            Estep_numeritor = ProbZ(i) * prod(ProbR(Recommended, i)) * prod(1 - ProbR(NotRecommended, i));
            Posteriors(i, t) = Estep_numeritor / Estep_denominator;
        end
    end

    L(iter+1) = L(iter+1) / T; % Normalization;
    
    % M-step.
    ProbZ_updated = zeros(1, k);
    ProbR_updated = zeros(N, k);

    for i = 1:k
        temp = sum(Posteriors(i, :));
        ProbZ_updated(i) = temp/T;

        for j = 1:N
            % Calculate numerator for M-Step.
            t_seen = find(Ratings(:, j) == 1);
            numer_seen = sum(Posteriors(i, t_seen));
            
            t_unseen = find(isnan(Ratings(:, j)));
            numer_unseen = ProbR(j, i) * sum(Posteriors(i, t_unseen));
            
            ProbR_updated(j, i) = (numer_seen + numer_unseen) / temp;
        end
    end
    
    % Update probabilities.
    ProbZ = ProbZ_updated;
    ProbR = ProbR_updated;
    
end

% Display the results. 
Results = array2table([(0:1:MaxIter)' L], 'VariableNames', {'Iteration', 'Log-likelihood'});
DispIterations = [0 1 2 4 8 16 32 64 128 256];
disp('(e) Normalized log-likelihood (note the log-likelihood increases at each iteration):')
disp(Results(DispIterations+1,:));

% Personal movie recommendations. -----------------------------------------
StudentIndex = find(StudentID == 'A59026267');

% Compute expected ratings for unseen movies.
UnseenMovies = find(isnan(Ratings(StudentIndex, :)));
ExpectedRatings = zeros(length(UnseenMovies), 1);

for i = 1:length(UnseenMovies) % Loop over all unseen movies.
    MovieIndex = UnseenMovies(i); % Index of an unseen movie in the movie list.
    PosteriorDenominator = 0;
    Posterior = zeros(k, 1);

    for j = 1:k
        RatedGood = find(Ratings(StudentIndex, :) == 1);
        RatedBad = find(Ratings(StudentIndex, :) == 0);
        
        P_good = prod(ProbR(RatedGood, j));
        P_bad = prod(1-ProbR(RatedBad, j));

        % Calculate posterior probability.
        Posterior(j) = ProbZ(j) * P_good * P_bad;
        PosteriorDenominator = PosteriorDenominator + Posterior(j);
    end

    % Normalize the posterior probability.
    Posterior = Posterior/PosteriorDenominator;

    % Calculate the expected rating for the unseen movie.
    for j = 1:k
        ExpectedRatings(i) = ExpectedRatings(i) + Posterior(j) * ProbR(MovieIndex, j);
    end
end

% Display the results.
UnseenMovieRatings = array2table(ExpectedRatings, 'VariableNames', {'Expected Rating'});
UnseenMovieRatings = addvars(UnseenMovieRatings, MovieList(UnseenMovies), 'Before', 1, 'NewVariableNames', 'Unseen Movie');
UnseenMovieRatingsSorted = sortrows(UnseenMovieRatings, 'Expected Rating', 'descend');
disp('(f) Unseen movies with expected ratings:')
disp('Compare to the mean popularity ratings in Part(a), this is much better.')
disp(UnseenMovieRatingsSorted);
