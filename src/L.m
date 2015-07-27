function [ y ] = L( q, l, tau, s )
%L Used to calculate eG
%   See equation 2.17 and 2.18 in Hawkes, Jalali, and Colquhoun (1990)
%   Phil Trans R Soc Lond A

id = eye(size(q(l,l)));

y = expm(-(s*id-q(l,l))*tau);

end

