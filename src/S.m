function [ y ] = S( q, l, tau, s )
%S Used to calculate eG
%   See equation 2.16 and 2.18 in Hawkes, Jalali, and Colquhoun (1990)
%   Phil Trans R Soc Lond A

id = eye(size(q(l,l)));

y = id - expm(-(s*id-q(l,l))*tau);

end

