function [ phib, ef ] = CHSvectors( Q, A, F, areaR, mu, res, tcrit ) %#codegen
%CHSVECTORS Gives the initial and final vectors, phib and ef, respectively
%used for maximum likelihood estimation of rate constants from bursts of
%channel activity
%   Colquhoun, Hawkes and Srodzinski (1996) Phil Trans R Soc Lond A
%   equations 5.8 and 5.11
%   
%   Q - Q matrix
%   A - vector of open (or shut) states
%   F - vector of shut (or open) states
%   areaR - matrix returned from asymptoticRvals that is used to determine
%       the asymptotic approximation of the distribution of e-opens and
%       e-shuts. If A is the open states, then areaR should correspond to
%       the shut states and vice versa
%   mu - time constants for the asymptotic shut time distribution - this 
%       matrix is returned from asymptoticRvals
%   res - resolution (or dead time) imposed on the data
%   tcrit - critical gap length separating bursts of activity arising from
%       a single ion channel

kF = numel(F);
kA = numel(A);
uA = ones(kA,1);
const = Q(F,A)*expm(Q(A,A)*res);
time = tcrit-res;

phiF = phi(Q,F,A,res);

Hfa = zeros(kF,kA);
for ii=1:kF
    Hfa = Hfa + areaR(:,:,ii)*const*mu(ii)*exp(-time/mu(ii));
end

ef = Hfa*uA;
phib = (phiF*Hfa) / (phiF*Hfa*uA);

end

