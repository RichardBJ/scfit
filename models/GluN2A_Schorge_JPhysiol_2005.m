%% Model of GluN2A Activation from Schorge et al 2005 J Physiol

%% Pack the Q matrix
% GluN2A Rates in s^-1 or M^-1*s^-1
alpha1 = 734;
beta1 = 35200;
alpha2 = 11100;
beta2 = 788;
don = 1130;
doff = 925;
gma = 22500;
gpa = 2013;
gmb = 150;
gpb = 66.1;
kmg = 29.4;
kpg = 1.03e7;
kme = 25.4;
kpe = 8.04e6;

% Q matrix
%     c     n
%    / \    |
%   b   f   j
%  / \ / \ / \
% a   e   i   l-m
%  \ / \ / \ /
%   d   h   k
%    \ /    |
%     g     z

q = zeros(15);
q(1,2) = 2*kpe;
q(1,4) = 2*kpg;
q(2,1) = kme;
q(2,3) = kpe;
q(2,5) = 2*kpg;
q(3,2) = 2*kme;
q(3,6) = 2*kpg;
q(4,1) = kmg;
q(4,5) = 2*kpe;
q(4,7) = kpg;
q(5,2) = kmg;
q(5,4) = kme;
q(5,6) = kpe;
q(5,8) = kpg;
q(6,3) = kmg;
q(6,5) = 2*kme;
q(6,9) = kpg;
q(7,4) = 2*kmg;
q(7,8) = 2*kpe;
q(8,5) = 2*kmg;
q(8,7) = kme;
q(8,9) = kpe;
q(9,8) = 2*kme;
q(9,6) = 2*kmg;
q(9,10) = gpb;
q(9,11) = gpa;
q(10,9) = gmb;
q(10,14) = don;
q(10,12) = gpa;
q(11,9) = gma;
q(11,12) = gpb;
q(11,15) = beta2;
q(12,10) = gma;
q(12,11) = gmb;
q(12,13) = beta1;
q(13,12) = alpha1;
q(14,10) = doff;
q(15,11) = alpha2;

% open states
A = [1 2];

% shut states
F = 3:15;

% scale rates to ms^-1
q = q*1e-3;

% Concentration-dependence of the rates
c = zeros(15);
idx = sub2ind(size(c),[1 2 4 5 7 8],[2 3 5 6 8 9]);
c(idx)=1;
idx2 = sub2ind(size(q),[1 4 2 5 3 6],[4 7 5 8 6 9]);
c(idx2) = 2;

%% Constraints on the rate constants
% Identify the indices of the rates in q matrix
idxall = find(q);
idxvary = sub2ind(size(q), [15 15 14 12 7 7 6 5 5 2 6 3 1], [14 12 15 15 6 5 7 7 2 5 3 6 4]);

% Constrain glutamate binding rates
src = sub2ind(size(q), [15 15 15 15 15], [14 14 14 14 14]);
tgt = sub2ind(size(q), [12 9 14 11 8], [11 8 13 10 7]);
[gamma,xi] = constrainrate (q, idxall, 'constrain', src, tgt, [1 1 0.5 0.5 0.5]);

% Constrain glutamate unbinding rates
src = sub2ind(size(q), [14 14 14 14 14], [15 15 15 15 15]);
tgt = sub2ind(size(q), [11 8 13 10 7], [12 9 14 11 8]);
[tmpgamma, tmpxi] = constrainrate(q,idxall,'constrain',src,tgt,[1 1 2 2 2]);
gamma=[gamma;tmpgamma];
xi=[xi;tmpxi];

% Constrain glycine binding rates
src = sub2ind(size(q), [15 15 15 15 15], [12 12 12 12 12]);
tgt = sub2ind(size(q), [14 13 12 11 10], [11 10 9 8 7]);
[tmpgamma, tmpxi] = constrainrate(q,idxall,'constrain',src,tgt,[1 1 0.5 0.5 0.5]);
gamma=[gamma;tmpgamma];
xi=[xi;tmpxi];

% Constrain glycine unbinding rates
src = sub2ind(size(q), [12 12 12 12 12], [15 15 15 15 15]);
tgt = sub2ind(size(q), [11 10 9 8 7], [14 13 12 11 10]);
[tmpgamma, tmpxi] = constrainrate(q,idxall,'constrain',src,tgt,[1 1 2 2 2]);
gamma=[gamma;tmpgamma];
xi=[xi;tmpxi];

% Constrain gating rates
src = sub2ind(size(q), [7 7 6 5], [6 5 7 7]);
tgt = sub2ind(size(q), [5 6 4 4], [4 4 5 6]);
[tmpgamma, tmpxi] = constrainrate(q,idxall,'constrain',src,tgt,[1 1 1 1]);
gamma=[gamma;tmpgamma];
xi=[xi;tmpxi];

% Fix the opening rate
src = sub2ind(size(q), 4, 1);
[tmpgamma, tmpxi] = constrainrate(q,idxall,'fix',src,[],35);
gamma=[gamma;tmpgamma];
xi=[xi;tmpxi];
