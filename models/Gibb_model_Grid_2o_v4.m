%% NMDA Receptor Model - Gibb_model_Grid_2o_v4
% Build a model of NMDA receptor activation
%
alpha1 =	2100	; %	q(1 ,2)
beta1 =	3400	; %	q(2 ,1 )
alpha2 =	500	; %	q(2 ,3 )
beta2 =	2000	; %	q(3 ,2 )
kfp =	1000	; %	q(9,6 )
kfm =	400	; %	q(9,11 )
ksp =	25	; %	q(10,8)
ksm =	5	; %	q(10,11)
%
q=zeros(11);
%
q(1,2) = alpha1;
q(2,1) = beta1;
q(2,3) = alpha2;
q(3,2) = beta2;
q(3,4) = 2*ksm;
q(4,3) = ksp;
q(3,5) = 2*kfm;
q(5,3) = kfp;
q(4,6) = ksm;
q(6,4) = 2*ksp;
q(4,7) = 2*kfm;
q(7,4) = kfp;
q(5,7) = 2*ksm;
q(7,5) = ksp;
q(5,8) = kfm;
q(8,5) = 2*kfp;
q(6,9) = 2*kfm;
q(9,6) = kfp;
q(7,9) = ksm;
q(9,7) = 2*ksp;
q(7,10) = kfm;
q(10,7) = 2*kfp;
q(8,10) = 2*ksm;
q(10,8) = ksp;
q(9,11) = kfm;
q(11,9) = 2*kfp;
q(10,11) = ksm;
q(11,10) = 2*ksp;

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
idxvary = sub2ind(size(q), [1 2 2 3 9 10 9 10], [2 1 3 2 6 8 11 11]);

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