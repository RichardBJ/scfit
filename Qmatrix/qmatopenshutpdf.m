function [opentaus, openareas, shuttaus, shutareas] = qmatopenshutpdf(q,ostates,cstates)
%QMATOPENSHUTPDF Calculates the taus and corresponding areas
%for the mixture of exponential distributions decribing the equilibrium 
%open times and shut times for the gating mechanism described by the Q matrix
%   q - the full Q matrix
%   ostates - vector containing the indices of the open states
%   cstates - vector containing the indices of the shut states

qopens = q(ostates,ostates);
qshuts = q(cstates,cstates);
qos = q(ostates,cstates);
qso = q(cstates,ostates);
peq = eqoccupy(q);
po = peq(ostates);
ps = peq(cstates);
uo = ones(length(ostates),1);
us = ones(length(cstates),1);

%open time distribution
phio = (ps*qso)./(ps*qso*uo);
[opentaus, ~, A] = qmatvals(qopens);
opentaus=opentaus';
openareas = zeros(1,length(ostates));
for ii=1:length(ostates)
    openareas(ii) = -opentaus(ii)*phio*A(:,:,ii)*qopens*uo;
end

%shut time distribution
phis = (po*qos)./(po*qos*us);
[shuttaus, ~, A] = qmatvals(qshuts);
shuttaus=shuttaus';
shutareas = zeros(1,length(cstates));
for ii=1:length(cstates)
    shutareas(ii) = -shuttaus(ii)*phis*A(:,:,ii)*qshuts*us;
end

end

