%% Model of GluN2C Activation from Dravid et al 2008 J Physiol
% From Scheme 3a

%% Pack the Q matrix
% GluN2C Rates in s^-1 or uM^-1*s^-1
%
%    D1    RA2b
%     \   /   \
% R-RA-RA2a   RA2d-O1-O2
%     /   \   /
%    D2    RA2c

ksp = 65;
ksm = 43;
kfp = 570;
kfm = 3300;
kd1p = 97;
kd1m = 22;
kd2p = 0.12;
kd2m = 1.7;
kon = 1.4;
koff = 8.8;
k3p = 590;
k3m = 3400;
k4p = 1300;
k4m = 2500;
 
q=zeros(10);
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
q(6,7)=k3p;
q(7,6)=k3m;
q(7,8)=k4p;
q(8,7)=k4m;
q(3,9)=kd1p;
q(9,3)=kd1m;
q(3,10)=kd2p;
q(10,3)=kd2m;

% open state
A = 7:8;

% shut states
F = [1:6, 9:10];

% scale rates to ms^-1
q = q*1e-3;

% Concentration-dependence of rates
c=zeros(10);
c(1,2)=1;
c(2,3)=1;

%% Constraints on the rate constants
% Identify the indices of the rates in Q matrix
idxall = find(q)';
idxconstrain = sub2ind(size(q), [2 3 5 6 4 6], [3 2 6 5 6 4]);
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
