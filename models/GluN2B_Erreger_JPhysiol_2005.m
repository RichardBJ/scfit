%% Model of GluN2B Activation from Erreger et al 2005 J Physiol

%% Pack the Q matrix
% GluN2B Rates in s^-1 or uM^-1*s^-1
ksp = 48;
ksm = 230;
kfp = 2836;
kfm = 175;
kd1p = 550;
kd1m = 81.4;
kd2p = 112;
kd2m = 0.91;
kon = 2.83;
koff = 38.1;

% Q matrix
q=zeros(8);
q(1,2)=2*kon;
q(2,1)=koff;
q(2,3)=kon;
q(3,2)=2*koff;
q(3,4)=kfp;
q(4,3)=kfm;
q(3,5)=ksp;
q(5,3)=ksm;
q(4,6)=ksp;
q(6,4)=ksm;
q(5,6)=kfp;
q(6,5)=kfm;
q(3,7)=kd1p;
q(7,3)=kd1m;
q(3,8)=kd2p;
q(8,3)=kd2m;

% open state
A = 6;

% shut states
F = [1:5, 7:8];

% scale rates to ms^-1
q = q*1e-3;

% Concentration-dependence of rates
c=zeros(8);
c(1,2)=1;
c(2,3)=1;

%% Constraints on the rate constants
% Identify the indices of the rates in Q matrix
idxall = sub2ind(size(q), [1 2 2 3 3 3 3 3 4 4 5 5 6 6 7 8], [2 1 3 2 4 5 7 8 3 6 3 6 4 5 3 3])';
idxconstrain = sub2ind(size(q), [2 3 4 6 5 6], [3 2 6 4 6 5]);
idxvary = setdiff(idxall,idxconstrain);

% Constrain binding and unbinding rates
src = sub2ind(size(q), [1 2], [2 1]);
tgt = sub2ind(size(q), [2 3], [3 2]);
[gamma,xi] = constrainrate (q, idxall, 'constrain', src, tgt, [0.5 2]);

% Constrain gating rates to have a fast step, kf, and a slow step, ks, in
% any order
src = sub2ind(size(q), [3 3 4 5], [4 5 3 3]);
tgt = sub2ind(size(q), [5 4 6 6], [6 6 5 4]);
[tmpgamma, tmpxi] = constrainrate(q,idxall,'constrain',src,tgt,[1 1 1 1]);
gamma=[gamma;tmpgamma];
xi=[xi;tmpxi];
