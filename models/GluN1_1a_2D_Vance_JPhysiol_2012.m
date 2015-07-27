%% Model of GluN1-1a/GluN2D Activation from Vance et al 2005 J Physiol

%% Pack the Q matrix
% GluN2D Rates in s^-1 or uM^-1*s^-1 from Vance et al (2012) J Physiol for GluN1-1a/GluN2D
bp = 0.25; % in uM^-1*s^-1
bm = 0.21;
k1p = 130;
k1m = 990;
k2p = 880;
k2m = 11000;
k3p = 6800;
k3m = 6200;
k4p = 7400;
k4m = 3600;
k5p = 4.6;
k5m = 3.4;
k6p = 140;
k6m = 480;
k7p = 410;
k7m = 9100;
k8p = 820;
k8m = 12000;
k9p = 5600;
k9m = 4500;

% Q matrix
q=zeros(14);
q(1,2)=2*bp;
q(1,8)=k5p;
q(2,1)=bm;
q(2,3)=bp;
q(2,9)=k5p;
q(3,2)=2*bm;
q(3,4)=k1p;
q(3,10)=k5p;
q(4,3)=k1m;
q(4,5)=k2p;
q(4,11)=k5p;
q(5,4)=k2m;
q(5,6)=k3p;
q(5,12)=k5p;
q(6,5)=k3m;
q(6,7)=k4p;
q(7,6)=k4m;
q(8,1)=k5m;
q(8,9)=2*bp;
q(9,2)=k5m;
q(9,8)=bm;
q(9,10)=bp;
q(10,3)=k5m;
q(10,9)=2*bm;
q(10,11)=k6p;
q(11,4)=k5m;
q(11,10)=k6m;
q(11,12)=k7p;
q(12,5)=k5m;
q(12,11)=k7m;
q(12,13)=k8p;
q(13,12)=k8m;
q(13,14)=k9p;
q(14,13)=k9m;

% open state
A = [6 7 13 14];

% shut states
F = [1:5, 8:12];

% scale rates to ms^-1
q = q*1e-3;

% Concentration-dependence of the rates
c = zeros(size(q));
c(1,2)=1;c(2,3)=1;c(8,9)=1;c(9,10)=1;

%% Constraints on the rate constants
% Identify the indices of the rates in q matrix
idxall=find(q);
idxvary = sub2ind(size(q), [1 1 2 3 4 4 5 5 8 10 11 11 12 12 6 6 7 13 13 14], [2 8 1 4 3 5 4 6 1 11 10 12 11 13 5 7 6 12 14 13])';

% Constrain binding and unbinding rates
src = sub2ind(size(q), [1 1 1 2 2 2], [2 2 2 1 1 1]);
tgt = sub2ind(size(q), [2 8 9 3 9 10], [3 9 10 2 8 9]);
[gamma,xi] = constrainrate(q,idxall,'constrain',src,tgt,[0.5 1 0.5 2 1 2]);

% Constrain all rates going from the high Po arm to the low Po arm to be
% equal
src = sub2ind(size(q), [1 1 1 1 8 8 8 8], [8 8 8 8 1 1 1 1]);
tgt = sub2ind(size(q), [2 3 4 5 9 10 11 12], [9 10 11 12 2 3 4 5]);
[tmpg,tmpxi] = constrainrate(q,idxall,'constrain',src,tgt,[1 1 1 1 1 1 1 1]);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];
