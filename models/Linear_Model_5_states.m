%% 5 state linear model with 2 open states (4 and 5)
% 1 <-> 2 <-> 3 <-> 4 <-> 5

%% Pack the Q matrix
% Rates in s^-1 or uM^-1*s^-1
k12 = 1e3;
k21 = 0.5e3;
k23 = 1.1e3;
k32 = 0.7e3;
k34 = 550;
k43 = 81.4;
k45 = 1.2e3;
k54 = 0.91e3;

% Q matrix
q=zeros(5);
q(1,2)=k12;
q(2,1)=k21;
q(2,3)=k23;
q(3,2)=k32;
q(3,4)=k34;
q(4,3)=k43;
q(4,5)=k45;
q(5,4)=k54;

% open state
A = [4,5];

% shut states
F = [1:3];

% scale rates to ms^-1
q = q*1e-3;

% Concentration-dependence of rates
c=zeros(size(q));

%% Constraints on the rate constants
% Identify the indices of the rates in Q matrix
idxall = find(q);
% idxconstrain = sub2ind(size(q), [2 3 4 6 5 6], [3 2 6 4 6 5]);
idxconstrain = [];
idxvary = setdiff(idxall,idxconstrain);

% Constrain binding and unbinding rates
% src = sub2ind(size(q), [1 2], [2 1]);
% tgt = sub2ind(size(q), [2 3], [3 2]);
% [gamma,xi] = constrainrate (q, idxall, 'constrain', src, tgt, [0.5 2]);

% Constrain gating rates to have a fast step, kf, and a slow step, ks, in
% any order
% src = sub2ind(size(q), [3 3 4 5], [4 5 3 3]);
% tgt = sub2ind(size(q), [5 4 6 6], [6 6 5 4]);
% [tmpgamma, tmpxi] = constrainrate(q,idxall,'constrain',src,tgt,[1 1 1 1]);
% gamma=[gamma;tmpgamma];
% xi=[xi;tmpxi];
