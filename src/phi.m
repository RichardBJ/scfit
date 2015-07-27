function [ p ] = phi( q, l, m, td )
%PHI Equilibrium vector for the Markov chain embedded by the semi-Markov
%process that occurs when the Markov process with transition rate matrix Q
%goes from set l to set m
%   See equations 2.14 and 2.14 of Hawkes, Jalali, and Colquhoun (1990)
%   Phil Trans R Soc Lond A
%   Also see Hawkes and Sykes (1990) IEEE Trans on Reliability

p = [zeros(1,length(l)), 1] / [eye(length(l))-eG(q,l,m,td,0)*eG(q,m,l,td,0), ones(length(l),1)];
end

