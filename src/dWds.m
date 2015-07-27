function [ dWds ] = dWds( s, q, A, F, tau )
%DWDS Calculates the derivative of the matrix function W(s)
%   The derivative of W(s) is used in the asymptotic approximation to the
%   matrix function R, which governs the pdf of dwell times when events are
%   missed (see Hawkes et al 1990)
%
%   This definition of dWds is from Hawkes et al (1992) pp. 394, eq. (56)

qAA = q(A,A);
qAF = q(A,F);
qFF = q(F,F);
qFA = q(F,A);
idA = eye(size(qAA));
idF = eye(size(qFF));
M = s*idF - qFF;
%From eq (2.16) in Hawkes et al (1990), SFF*(s) = I-expm(-(s*I-Q(F,F))*tau)
SFF = idF - expm(-M*tau);
%From eq (4) in Hawkes et al (1992), GFA*(s) = inv(s*I-Q(F,F))*Q(F,A)
GFA = M\qFA;

dWds = idA + qAF*(SFF/M - tau*(idF-SFF))*GFA;

end

