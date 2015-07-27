function [ W ] = W( s, q, A, F, tau ) %#codegen
%W Returns the matrix W from Hawkes et al (1992) Phil Trans R Soc Lond B p.393
%   The asymptotic behavior of the matrix function, R, depends on values of
%   s that render the matrix W(s) singular.  W(s) is defined as
%       $$W(s) = sI - H(s)$$
%   where
%       $$H(s) = Q_{AA}+Q_{AF}\left (\int_0^{\tau} e^{-st}e^{Qt}dt \right )
%       Q_{FA}$$
%   or, if s is not an eigenvalue of Q(F,F), then inv(s*I - Q(F,F)) exits
%   and
%       H(s) = Q(A,A) + Q(A,F)*inv(s*I - Q(F,F))*(I - exp(-(s*I-Q(F,F))*tau))*Q(F,A)

tol = 1e-12;
qAA = q(A,A);
idA = eye(size(qAA));
qAF = q(A,F);
qFF = q(F,F);
qFA = q(F,A);
idF = eye(size(qFF));
eigFF = eig(qFF);

M = s*idF - qFF;
%check whether (sI-qFF)^-1 exists, which will not occur if sI-qFF is singular, 
%if sI-qFF is singular, the det(sI-qFF)=0 or equivalently s is an eigenvalue of qFF
if any(abs(eigFF-s)<tol)
    % since (sI-qFF)^-1 does not exist, we cannot shortcut calculation of
    % the integral e^(-(sI-qFF)*t), so let's estimate it numerically
    x = linspace(0,tau);
    y = zeros([size(M),numel(x)]);
    for ii=1:numel(x)
        y(:,:,ii) = real(expm(-x(ii)*M));
    end
    % Numerically estimate the integral using trapz()
    H1 = trapz(x,y,3);
    H = qAA + qAF*H1*qFA;
    fprintf('In W(s)\ns is an eigenvalue of qFF, so (sI-qFF)^-1 does not exist.');
else
    H = qAA + qAF/M*(idF - expm(-M*tau))*qFA;
end

W = s*idA - H;

end

