%% *CSE 250A HW-4*
% *Problem 4.2: Statistical language modeling* 
% *(a) Compute the maximum likelihood estimate of the unigram distribution _Pu_(_w_) 
% over words _w_. Print out a table of all the tokens (i.e., words) that start 
% with the letter “M”, along with their numerical unigram probabilities (not counts).*

clear; clc; close all;

% Read the vocabulary data file.
vocab = readlines('hw4_vocab.txt'); 
vocab(end) = [];

% Read the unigram counts data file and convert it into double precision numbers.
unigram = str2double(readlines('hw4_unigram.txt'));
unigram(end) = [];

% Calculate unigram probability.
unigram_prob = unigram/sum(unigram);

% Save the unigram probability result in a table.
vocab_unigram = table(vocab, unigram, unigram_prob, 'VariableNames', {'Word', 'Unigram Count', 'Unigram Probability'});

% Show the results for words starting with letter M.
vocab_unigram_M = vocab_unigram(startsWith(vocab_unigram.Word, 'M'), :);
vocab_unigram_M.Word = char(vocab_unigram_M.Word);
disp(vocab_unigram_M);
% *(b) Compute the maximum likelihood estimate of the bigram distribution _Pb_(_w_′|_w_). Print out a table of the 10 most likely words to follow the word “THE”, along with their numerical bigram probabilities.*

% Input the word that other words will follow (Word1).
Word1 = 'THE';

% Read the bigram counts data file into a table.
bigram = readtable('hw4_bigram.txt', 'Delimiter', '\t', 'ReadVariableNames', false, 'TextType', 'string');
bigram.Properties.VariableNames = {'Index(w1)', 'Index(w2)', 'Count(w1, w2)'};

% Find the index of Word1.
Index_Word1 = find(strcmp(vocab, Word1));

% Extract the bigram data for Word1 and calculate the bigram probability.
bigram_Word1 = bigram(bigram{:,1} == Index_Word1, :);
bigram_Word1.('Probability(w1, w2)') = bigram_Word1{:,3} / sum(bigram_Word1{:,3});

% Sort the bigram data based on the bigram probability.
bigram_Word1 = sortrows(bigram_Word1, 'Probability(w1, w2)', 'descend');

% Find the corresponding Word2 based on the index.
EmptyColumns = table(strings(height(bigram_Word1),1), strings(height(bigram_Word1),1), 'VariableNames', {'Word 1', 'Word 2'});
bigram_Word1 = [EmptyColumns bigram_Word1];

for i = 1:size(bigram_Word1,1)
    bigram_Word1{i,1} = vocab(bigram_Word1{i,3});
    bigram_Word1{i,2} = vocab(bigram_Word1{i,4});
end

% Display the result.
bigram_Word1.('Word 1') = char(bigram_Word1.('Word 1'));
bigram_Word1.('Word 2') = char(bigram_Word1.('Word 2'));
disp(bigram_Word1(1:10,:));
% *(c) Consider the sentence “The stock market fell by one hundred points last week.” Ignoring punctuation,compute and compare the log-likelihoods of this sentence under the unigram and bigram models. In the equation for the bigram log-likelihood, the token ⟨s⟩ is used to mark the beginning of a sentence. Which model yields the highest log-likelihood?*
% _Unigram model:_

% Input the sentence.
Input = 'THE STOCK MARKET FELL BY ONE HUNDRED POINTS LAST WEEK';

% Split the input sentence into multiple words.
Input = string(split(Input));

% Initialize the log-likelihood with zero.
Lu = 0;

% Initialize the error index with zero.
Error = 0;

for i = 1:size(Input,1)
    WordIndex = find(strcmp(vocab, Input(i))); % Find the index in the vocabulary list for the i-th word.

    if isempty(WordIndex) % If the case is not find, then jump out of the loop and return to -inf.
        Error = 1;
        break;
    end

    Lu = Lu + log(vocab_unigram{WordIndex,3}); % Update the log-likelihood.
end

% Display the result.
if Error == 0
    fprintf('Log-likelihood under the unigram model =  %f \n', Lu);
else
    fprintf('Log-likelihood under the unigram model =  -inf \n');
end
%% 
% _Bigram model:_

% Input the sentence.
Input = '<s> THE STOCK MARKET FELL BY ONE HUNDRED POINTS LAST WEEK';

% Split the input sentence into multiple words.
Input = string(split(Input));

% Initialize the log-likelihood with zero.
Lb = 0;

for i = 1:size(Input,1)-1
    Word1Index = find(strcmp(vocab, Input(i))); % Find the index in the vocabulary list for the i-th word.
    Word2Index = find(strcmp(vocab, Input(i+1))); % Find the index in the vocabulary list for the word following the i-th word.

    if isempty(Word1Index) || isempty(Word2Index)  % If the case is not find, then jump out of the loop and return to -inf.
        Error = 1;
        break;
    end

    P_Word1 = find(bigram{:,1} == Word1Index); 
    P_Word2_Word1 = find(bigram{:,1} == Word1Index & bigram{:,2} == Word2Index);

    if isempty(P_Word2_Word1)  % If the case is not find, then jump out of the loop and return to -inf.
        Error = 1;
        break;
    end

    Lb = Lb + log(bigram{P_Word2_Word1,3}/sum(bigram{P_Word1,3})); % Update the log-likelihood.
end

% Display the result.
if Error == 0
    fprintf('Log-likelihood under the unigram model =  %f \n', Lb);
else
    fprintf('Log-likelihood under the unigram model =  -inf \n');
end
%% 
% The bigram model yields a higher log-likelihood.
% *(d) Consider the sentence “The sixteen officials sold fire insurance.” Ignoring punctuation, compute and compare the log-likelihoods of this sentence under the unigram and bigram models. Which pairs of adjacent words in this sentence are not observed in the training corpus? What effect does this have on the log-likelihood from the bigram model?*
% _Unigram model:_

% Input the sentence.
Input = 'THE SIXTEEN OFFICIALS SOLD FIRE INSURANCE';

% Split the input sentence into multiple words.
Input = string(split(Input));

% Initialize the log-likelihood with zero.
Lu = 0;

% Initialize the error index with zero.
Error = 0;

for i = 1:size(Input,1)
    WordIndex = find(strcmp(vocab, Input(i))); % Find the index in the vocabulary list for the i-th word.

    if isempty(WordIndex) % If the case is not find, then jump out of the loop and return to -inf.
        Error = 1;
        break;
    end

    Lu = Lu + log(vocab_unigram{WordIndex,3}); % Update the log-likelihood.
end

% Display the result.
if Error == 0
    fprintf('Log-likelihood under the unigram model =  %f \n', Lu);
else
    fprintf('Log-likelihood under the unigram model =  -inf \n');
end
%% 
% _Bigram model:_

% Input the sentence.
Input = '<s> THE SIXTEEN OFFICIALS SOLD FIRE INSURANCE';

% Split the input sentence into multiple words.
Input = string(split(Input));

% Initialize the log-likelihood with zero.
Lb = 0;

for i = 1:size(Input,1)-1
    Word1Index = find(strcmp(vocab, Input(i))); % Find the index in the vocabulary list for the i-th word.
    Word2Index = find(strcmp(vocab, Input(i+1))); % Find the index in the vocabulary list for the word following the i-th word.

    if isempty(Word1Index) || isempty(Word2Index)  % If the case is not find, then jump out of the loop and return to -inf.
        Error = 1;
        break;
    end

    P_Word1 = find(bigram{:,1} == Word1Index); 
    P_Word2_Word1 = find(bigram{:,1} == Word1Index & bigram{:,2} == Word2Index);

    if isempty(P_Word2_Word1)  % If the case is not find, then jump out of the loop and return to -inf.
        Error = 2;
        break;
    end

    Lb = Lb + log(bigram{P_Word2_Word1,3}/sum(bigram{P_Word1,3})); % Update the log-likelihood.
end

% Display the result.
if Error == 0
    fprintf('Log-likelihood under the unigram model =  %f \n', Lb);
else
    fprintf('Log-likelihood under the unigram model =  -inf \n');
end

if Error == 1
    disp(['The word does not exist: ', vocab(Word1Index)]);
end

if Error == 2
   fprintf('The first word sequence does not exist: %s %s\n', vocab(Word1Index), vocab(Word2Index));
end
%% 
% The first pair of adjacent words in this sentence not observed in the training 
% set is "SIXTEEN OFFICIALS". Also, it is observed that the "SOLD FIRE" is not 
% observed by checking the conditional probability "P_Word2_Word1" in the program.
% 
% The log-likelihood of the bigram model returns to minus infinity if the case 
% is not observed in the training set.
% (e) Mixture model.

% Input the sentence.
Input = '<s> THE SIXTEEN OFFICIALS SOLD FIRE INSURANCE';

% Split the input sentence into multiple words.
Input = string(split(Input));

% Define lambda values.
lambda = 0:0.01:1;

% Initialize the log-likelihood with zero.
Lm = zeros(length(lambda),1);

for j = 1:size(lambda,2)
    for i = 1:size(Input,1)-1
        Word1Index = find(strcmp(vocab, Input(i))); % Find the index in the vocabulary list for the i-th word.
        Word2Index = find(strcmp(vocab, Input(i+1))); % Find the index in the vocabulary list for the word following the i-th word.
    
        P_Word1 = find(bigram{:,1} == Word1Index);
        P_Word2_Word1 = find(bigram{:,1} == Word1Index & bigram{:,2} == Word2Index);
    
        Pu = vocab_unigram{WordIndex,3};
        Pb = bigram{P_Word2_Word1,3}/sum(bigram{P_Word1,3});

        Pm = lambda(j)*Pu + (1-lambda(j))*Pb;
        Lm(j) = Lm(j) + log(Pb); % Update the log-likelihood.
    end
end