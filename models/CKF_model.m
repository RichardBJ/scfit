%% CKF Ion Channel Mechanism
% R <--> AR <--> AR*

q = zeros(3);
q = [-1, 1, 0;
    19, -29, 10;
    0, 0.026, -0.026]; %rates in ms^-1
tmp = [1 2; 2 1; 2 3; 3 2];
idxall = sub2ind(size(q),tmp(:,1),tmp(:,2));
idxconstrain = [];
idxvary = setdiff(idxall,idxconstrain);
A = 1;
F = 2:3;

q = q-diag(sum(q,2));