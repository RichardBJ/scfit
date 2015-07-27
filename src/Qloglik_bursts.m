function ll = Qloglik_bursts (params, Q, idxtheta, M, b, A, F, ...
    td, tcrit, dwells) %#codegen
%QLOGLIK_Bursts Log likelihood of rates in a gating mechanism given the
%exact sequence of open and closed durations
%   params - log10(rates to vary)
%   Q - Q matrix
%   idxAll - vector giving the index of all the rate constants in Q
%   idxvarytoall - vector of indices such that idxAll(idxvarytoall) = idxVary
%   idxconstraintoall - similar to idxvarytoall
%   R1 - square matrix from qr factorization of gamma, the coefficient
%   matrix giving the linear constraints on log10(rates)
%   R2 - matrix corresponding to free parameters; given from qr
%   factorization of gamma
%   U - square unitary matix from qr factorization of gamma
%   xi - constants for linear constraints such that gamma*theta = xi
%   A - vector of the open states
%   F - vector of the shut states
%   td - resolution (or dead time) imposed on the data
%   tcrit - critical gap duration
%   dwells - vector of durations in each state. Must alternate starting 
%       from set A to set F and end on a dwell in set F

if any(isnan(params))
    ll=nan;
    return;
end

nBursts = size(dwells,2);
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

% Calculation of the forward recursions, which are used to calculate the
% likelihood, will run into numerical problems because the probabilities
% will almost always be less than one and the likelihood will therefore
% decay with increasing number of dwell times.  To avoid this, we
% continually scale the inital probability vector (see Rabiner 1989
% and Qin et al (1997) Proc R Soc Lond B

[phib, ef] = CHSvectors(qnew,A,F,areaRs,-1./s_s,td,tcrit);
LL = zeros(nBursts,1);
for ii=1:nBursts
    tmpdwells = dwells(:,ii);
    tmpdwells(isnan(tmpdwells))=[];
    numdwells = numel(tmpdwells);
    pA = phib;
    pF = zeros(1, length(F));
    scalefactor = zeros(numdwells,1);
    for jj=1:numdwells
        if mod(jj,2)==1
            pF = pA * R(Co,lambda,td,so,areaRo,mMax,abs(tmpdwells(jj)-td))*qnew(A,F)*eqFFt;
            scalefactor(jj) = 1./sum(pF);
            pF = pF .* scalefactor(jj);
        else
            pA = pF * R(Cs,lambda,td,s_s,areaRs,mMax,abs(tmpdwells(jj)-td))*qnew(F,A)*eqAAt;
            scalefactor(jj) = 1./sum(pA);
            pA = pA .* scalefactor(jj);
        end
    end
    LL(ii) = -sum(log10(scalefactor)) + log10(pF*ef);
end
ll = -sum(LL);

end