%% CH82 Ion Channel Mechanism

q = zeros(5);
tmp = [1 4; 4 1; 1 2; 2 1; 2 3; 3 2; 3 4; 4 3;4 5; 5 4];
idxall = sub2ind(size(q),tmp(:,1),tmp(:,2));
idxconstrain = idxall([3,4]);
idxvary = setdiff(idxall,idxconstrain);
A = 1:2;
F = 3:5;

% Rates in s^-1 or M^-1*s^-1
alpha1 = 3000;
beta1 = 15;
kp2star = 5e8;
km2star = 0.666667;
alpha2 = 500;
beta2 = 15000;
km2 = 4000;
kp2 = 5e8;
km1 = 2000;
kp1 = 1e8;
xA = 0.1e-6; % 0.1 uM

q(idxall) = [alpha1, beta1, kp2star*xA, km2star, alpha2, beta2, km2 kp2*xA, km1, kp1*xA];
q = q*1e-3;
q = q-diag(sum(q,2));

[gamma,xi] = constrainrate(q,idxall,'constrain',idxall(8),idxall(3),1);

src = sub2ind(size(q),[1 2 3 4],[2 3 4 1]);
tgt = sub2ind(size(q),[1 4 3 2],[4 3 2 1]);
[tmpg,tmpxi] = constrainrate(q,idxall,'loop',src,tgt);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];

clear tmp src tgt tmpg tmpxi;

guess1 = [10e3; 1e3; 5e7*xA; 0; 100; 30e3; 1e3; 5e7*xA; 1e3; 1e7*xA]*1e-3;
guess2 = [500; 3; 15e8*xA; 0; 100; 3e3; 10e3; 15e8*xA; 10e3; 5e8*xA]*1e-3;
guess3 = [1e3; 3; 5e7*xA; 0; 5e3; 50e3; 1e3; 5e7*xA; 1e3; 1e7*xA]*1e-3;