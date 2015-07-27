function [ tcrit ] = findtcrit( taus, weights )
%FINDTCRIT Find tcrit for the longest component of an exp mix distribution
%   Detailed explanation goes here

[t,ix] = sort(taus);
w = weights(ix);

t1 = t(end-1);
w1 = w(end-1);

t2 = t(end);
w2 = w(end);

fxn = @(x) w1/t1*exp(-x/t1) + w2/t2*exp(-x/t2) - w2/t2;

guess = (t1*t2/(t1-t2))*log((t1*w2)/(w1*t2));

tcrit = fzero(fxn,guess);

end

