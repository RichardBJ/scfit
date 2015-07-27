function [ y ] = eG( Q, l, m, tau, s )
%EG Laplace transform of the matrix describing the probability densities
%for intervals of the semi-Markov process that occur when the Markov
%process with the transition rate matrix Q transitions from set l to set m
%and all events of duration less than tau are missed
%   See equations 2.12, 2.19 and 2.20 of Hawkes, Jalali and Colquhoun
%   (1990) Phil Trans R Soc Lond A

id = eye(size(Q(l,l)));

y = (id - G(Q,l,m,s)*S(Q,m,tau,s)*G(Q,m,l,s))\(G(Q,l,m,s)*L(Q,m,tau,s));

end

