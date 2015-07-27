%% NMDA Receptor Model
% Build a model of NMDA receptor activation
q(1,2) = 200;
q(2,1) = 5000;
q(2,3) = 5000;
q(3,2) = 2000;
q(3,4) = 200;
q(4,3) = 100;
q(3,5) = 2000;
q(5,3) = 1000;
q(4,6) = 100;
q(6,4) = 200;
q(4,7) = 2000;
q(7,4) = 1000;
q(5,7) = 200;
q(7,5) = 100;
q(5,8) = 1000;
q(8,5) = 2000;
q(6,9) = 2000;
q(9,6) = 1000;
q(7,9) = 100;
q(9,7) = 200;
q(7,10) = 1000;
q(10,7) = 2000;
q(8,10) = 200;
q(10,8) = 100;
q(9,11) = 1000;
q(11,9) = 2000;
q(10,11) = 100;
q(11,10) = 200;

% openstates
A = 1:2;

% shutstates
F = 3:11;

% scale rates to ms^-1
q = q*1e-3;

% q = q-diag(sum(q,2));

%% Constraints
% Description of first code block
idxall = find(q);
idxvary = sub2ind(size(q), [1 2 2 3 9 10 9 10], [2 1 3 2 6 8 11 11])';

% Constrain rates k+s
src = sub2ind(size(q), 10*ones(1,5), 8*ones(1,5));
tgt = sub2ind(size(q), [4 6 7 9 11], [3 4 5 7 10]);
[gamma,xi] = constrainrate (q, idxall, 'constrain', src, tgt, [1 2 1 2 2]);

% Constrain rates k+f
src = sub2ind(size(q), 9*ones(1,5), 6*ones(1,5));
tgt = sub2ind(size(q), [5 7 8 10 11], [3 4 5 7 9]);
[tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, [1 1 2 2 2]);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];

% Constrain rates k-s
src = sub2ind(size(q), 10*ones(1,5), 11*ones(1,5));
tgt = sub2ind(size(q), [3 4 5 7 8], [4 6 7 9 10]);
[tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, [2 1 2 1 2]);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];

% Constrain rates k-f
src = sub2ind(size(q), 9*ones(1,5), 11*ones(1,5));
tgt = sub2ind(size(q), [3 4 5 6 7], [5 7 8 9 10]);
[tmpg,tmpxi] = constrainrate (q, idxall, 'constrain', src, tgt, [2 2 1 2 1]);
gamma=[gamma;tmpg];
xi=[xi;tmpxi];