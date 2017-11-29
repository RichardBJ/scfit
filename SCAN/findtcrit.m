function [ tcrit ] = findtcrit( taus, weights, varargin )
%FINDTCRIT Find tcrit for the longest component of an exp mix distribution
%   Input
%      taus
%          Means of the individual exponential components
%      weights
%          Weights of the exponential components (must sum to 1)
%      'method', value (optional)
%          The method used to define the critical time (see below). the
%          default is 2.
%
%   There is no unique criterion for the optimum way to divide an
%   experimental record into bursts but at least three methods have been
%   proposed
%   1. Minimize the _total number_ of misclassified intervals. Jackson et
%      al (1983)
%   2. Choose tcrit so that _equal numbers_ of short and long intervals are
%      misclassified. Magleby and Pallotta (1983) and Clapham and Neher
%      (1984)
%   3. Choose tcrit so that _equal proportions_ of short and long intervals
%      misclassfied. Colquhouh and Sakmann, 1985
%   This function uses method 2, equalizing the number of short and long
%   intervals misclassified.
%
% Reference
%   Colquhoun and Sigworth 1995. Fitting and Statistical Analysis of
%   Single-Channel Records. Ch. 19 pp. 535-536 in Single-Channel Recording
%   2nd Edition Editors Sakmann and Neher

p = inputParser;
addRequired(p, 'taus', @(x) isnumeric(x));
addRequired(p, 'weights', @(x) isnumeric(x) && sum(x) == 1);
addParameter(p, 'Method', 2, @(x) isnumeric(x) && isscalar(x) && any(x == [1, 2, 3]));
parse(p, taus, weights, varargin{:});

% Sort components so that the slowest component is last (has a larger tau)
[t,ix] = sort(taus);
w = weights(ix);

% The faster component (smaller tau)
t1 = t(end-1);
w1 = w(end-1);

% The slower component (larger tau)
t2 = t(end);
w2 = w(end);

if p.Results.Method == 1
    fxn = @(x) w1*(1/t1)*exp(-x/t1) - w2*(1/t2)*exp(-x/t2);
elseif p.Results.Method == 2
    fxn = @(x) w1*exp(-x/t1) + w2*exp(-x/t2) - w2;
elseif p.Results.Method == 3
    fxn = @(x) exp(-x/t1) + exp(-x/t2) - 1;
end

% Not sure where this guess came from
guess = (t1*t2/(t1-t2))*log((t1*w2)/(w1*t2));

tcrit = fzero(fxn,guess);

end

