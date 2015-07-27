function [ pdf ] = hjcdist( q, A, F, td, t, varargin ) %#codegen
%HJCDIST Returns the probability density function for the distribution
%of open (or shut) times for the gating mechanism given by the Q matrix, q
%   HJCDIST uses the exact correction for missed events given by Hawkes,
%   Jalali and Colquhoun (1990) using the asymptotic approximation given by
%   Hawkes, Jalali, and Colquhoun (1992)
%   INPUT
%       q - Q matrix
%       A - vector of the open (or shut) states. To get the distribution of
%       open times, pass the open states. To get the shut time distribution
%       pass the shut states.  Whichever one you don't pass here, pass to F
%       F - vector of the shut (or open) states
%       td - resolution imposed on the dwell times (aka dead time)
%       t - vector of times at which to return the pdf
%       islog (optional) - Boolean value indicating whether t is the
%       logarithm of the dwell times

pdf = zeros(length(t),1);

if nargin < 6
    islog = true;
else
    islog = varargin{1};
end

% Calculation of the exact distribution may be numerically unstable for
% mMax > 5. At least for mechanisms with many states. -KKO 140923
mMax = 3;
[C, lambda] = Cvals(q,A,F,td,mMax);
[s, areaR] = asymptoticRvals(q,A,F,td);
if any(isinf(s))
    error('Not all the roots of det(W) were found. The asymptotic approximation of R is unreliable.');
end
eqFFt = expm(q(F,F)*td);
uF = ones(length(F),1);
phiA = phi(q,A,F,td);

if islog
    t=10.^t;
end

for ii=1:length(t)
%     fprintf('In HJCDIST, t = %.4f\n', t(ii));
    if islog
        pdf(ii) = log(10).*t(ii).*phiA*R(C,lambda,td,s,areaR,mMax,t(ii)-td)*q(A,F)*eqFFt*uF;
    else
        pdf(ii) = phiA*R(C,lambda,td,s,areaR,mMax,t(ii)-td)*q(A,F)*eqFFt*uF;
    end
    if pdf(ii) <= 0
        ct = ii;
    end
end    

end

