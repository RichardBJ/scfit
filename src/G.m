function [ y ] = G( q, a, f, s)
%G Laplace transform of the matrix describing the probability densities
%for intervals of the semi-Markov process that occur when the Markov
%process with the transition rate matrix Q transitions from set l to set m
%   See equation 2.9 and 2.10 of Hawkes, Jalali, and Colquhoun (1990) 
%   Phil Trans R Soc Lond A

qaa = q(a,a);
id = eye(size(qaa));

y = (s*id - qaa)\q(a,f);

end

