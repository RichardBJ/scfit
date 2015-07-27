function ll = Qmatrixloglik (params, Q, idxtheta, M, b, A, F, td, dwells) %#codegen
%QMATRIXLOGLIK Log likelihood of rates in a gating mechanism given the
%exact sequence of open and closed durations
%   params - log10(rates to vary)
%   Q - Q matrix
%   idxtheta - vector of indices into the Q matrix that give all of the rate constants
%   M - coefficient matrix of theta such that M*parameters + b = theta
%   b - constant matrix
%   A - vector of the open states
%   F - vector of the shut states
%   td - resolution (or dead time) imposed on the data
%   dwells - vector of durations in each state. Must alternate starting 
%       from set A to set F and end on a dwell in set F

if any(isnan(params))
    ll=nan;
    return;
end

ndwells = length(dwells)./2;
nstates = length(Q);
qnew = zeros(nstates);

theta = M*params + b;
qnew(idxtheta) = 10.^theta;
qnew = qnew - diag(sum(qnew,2));

mMax = 2; % Number of multiples of tau to use for exact correction for missed events
[Co, lambda] = Cvals(qnew,A,F,td,mMax);
[so, areaRo] = asymptoticRvals(qnew,A,F,td);
if any(isinf(so))
    ll=nan;
    return;
end
Cs = Cvals(qnew,F,A,td,mMax);
[s_s, areaRs] = asymptoticRvals(qnew,F,A,td);
if any(isinf(s_s))
    ll=nan;
    return;
end
eqAAt = expm(qnew(A,A)*td);
eqFFt = expm(qnew(F,F)*td);
phiA = phi(qnew,A,F,td);

% Calculation of the forward recursions, which are used to calculate the
% likelihood, will run into numerical problems because the probabilities
% will almost always be less than one and the likelihood will therefore
% decay with increasing number of dwell times.  To avoid this, we
% continually scale the inital probability vector (see Rabiner 1989
% and Qin et al (1997) Proc R Soc Lond B

scalefactor = zeros(ndwells,1);
p = phiA;
for ii=1:ndwells
    time = dwells([2*ii-1, 2*ii]);
    p = p * R(Co,lambda,td,so,areaRo,mMax,abs(time(1)-td))*qnew(A,F)*eqFFt * R(Cs,lambda,td,s_s,areaRs,mMax,abs(time(2)-td))*qnew(F,A)*eqAAt;
    scalefactor(ii) = 1./sum(p);
    p = p .* scalefactor(ii);
end
ll = -(-sum(log10(scalefactor)));

end