function [U,R1,R2,idxtheta,Gamma,R,idxconstrain,nConstrain] = constrainqr(gamma,idxall,idxvary)
%CONSTRAINQR qr factorization for constraints

nRates = length(idxall);
nVary = length(idxvary);
nConstrain = nRates-nVary;

% Need to generate theta such that theta = log10([q(idxconstrain); q(idxvary)]);
% Generating theta doesn't actually need to be done here, but we need to
% arrange gamma correctly so that when we construct theta in this way then
% gamma * theta == xi
% Currently, gamma * log10(q(idxall)) == xi
% So we need to know how to go from idxall to idxtheta

idxconstrain = setdiff(idxall,idxvary);
assert(nConstrain==numel(idxconstrain),'Number of constraints does not match.');

idxtheta = [idxconstrain;idxvary]; % Note that theta == log10(q(idxtheta));
[~,iAllToTheta] = ismember(idxtheta,idxall);
Gamma = gamma(:,iAllToTheta);

% Now do the qr factorization of Gamma
[U,R] = qr(Gamma);
R1 = R(:,1:nConstrain);
R2 = R(:,nConstrain+1:end);

end