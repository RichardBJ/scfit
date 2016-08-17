%% Getting started fitting an exponential mixture distribution
% Illustrate the use of expmixdist functions

%% Simulate an exponential mixture distribution
% Three exponential components with means of 0.173, 2.84, and 12.5
% in the proportion 0.22, 0.61, 0.17
taus = [0.173, 2.84, 12.5];
weights = [0.22, 0.61, 0.17];
n = 5000;
wait_times = [exprnd(taus(1), floor(n*weights(1)), 1);
              exprnd(taus(2), floor(n*weights(2)), 1);
              exprnd(taus(3), floor(n*weights(3)), 1)];

%% Fit waiting times to an exponential mixture distribution

tau_guess = [0.1, 10, 1000];
weight_guess = ones(1, 3) ./ 3;
[ taushat, weightshat ] = emdistfit(wait_times, tau_guess, weight_guess);

%% Plot results

emhist(wait_times, taushat, weightshat);
